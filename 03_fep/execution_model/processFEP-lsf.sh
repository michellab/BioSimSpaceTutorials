#!/usr/bin/env bash
# = = = 
# Preamble
# We assume interactive execution of notebook generates
# ligands.dat
# network.dat
# protocol.dat
# a paramterised protein input stored in inputs/protein/protein.prm7[rst7]
# a list of ligand mol2 files stored in inputs/ligands/
# The main output will be stored in outputs/summary.csv
# = = =
# Stage 0 set cluster specific parameters
CPUQUEUE="phase8_slough"
GPUQUEUE="phase9_48_braine"
export BSSHOME="/home/model/MD-SOFTWARE/miniconda3-BSS"
export TMPDIR="./tmp/"
export OPENMM_PLUGIN_DIR="$BSSHOME/lib/plugins/"
export OMP_NUM_THREADS=1 # avoid oversubscribing CPUs with SOMD
source /home/model/MD-SOFTWARE/BSSenv.bashrc  # imports amber, gromacs


# Max run time for various categories of jobs. Should work for most use cases. 
PARAMTIME="01:00"
PREPTIME="01:00"
RUNTIME="04:00"
ANALYSISTIME="00:30"



# = = = 
# Jobs executed on the cluster
# for LSF, we need to parse the bsub output to get unique job ID. Use a function:
function nk_jobid {
	output=$($*)
	echo $output | head -n1 | cut -d'<' -f2 | cut -d'>' -f1
}

# Stage 1 - A GPU job array executed for each ligand listed in ligands.dat
num=$(< "ligands.dat" wc -l)
tasks=$(( $num ))
echo "@@@ Parameterising dataset @@@"
#echo $tasks
ID1=$(nk_jobid bsub -J "ligprep[1-$tasks]" -q $GPUQUEUE -gpu "num=1:mode=exclusive_process" -n 1 -W $PARAMTIME -o logs/ligprep_%J_%I.out scripts/BSSligprep.sh)

echo "bsub -J "ligprep[1-$tasks]" -q $GPUQUEUE -n 1 -W $PARAMTIME -o logs/ligprep_%J_%I.out scripts/BSSligprep.sh"
echo $ID1

# Stage 2 - A series of jobs are executed once all ligands have been prepared
mapfile PERTS < network.dat
for pert in "${PERTS[@]}"
do
    IFS=' ' read -a ligpair <<< "$pert"
    win=$(( ${ligpair[2]} ))
    echo "@@@ Processing ${ligpair[0]} ${ligpair[1]} $win @@@" 
    
    # Stage 2a. Generate FEP inputs with a single job over a single CPU
    ID2=$(nk_jobid bsub -J prepFEP -w "done($ID1)" -ti -q $CPUQUEUE -n 1 -W $PREPTIME -o logs/prepFEP_%J.out scripts/BSSprepFEP.sh ${ligpair[0]} ${ligpair[1]})

    echo "bsub -J prepFEP -w done($ID1) -ti -q $CPUQUEUE -n 1 -W $PREPTIME -o logs/prepFEP_%J.out scripts/BSSprepFEP.sh ${ligpair[0]} ${ligpair[1]}"
    echo $ID2

    # Stage 2b. Job array with 2*nwindows to run over free and bound legs over GPUs
    tasks=$(( 2*${ligpair[2]} )) # not -1 (as with the slurm script) accounts for LSF arrays starting at 1.
    ID3=$(nk_jobid bsub -J "runFEP[1-$tasks]" -w "done($ID2)" -ti -q $GPUQUEUE -gpu "num=1:mode=exclusive_process" -n 1 -W $RUNTIME -o logs/runFEP_%J.out scripts/BSSrunFEP.sh ${ligpair[0]} ${ligpair[1]} ${ligpair[3]} ${ligpair[4]})

    echo " bsub -J runFEP[1-$tasks] -w done($ID2) -ti -q $GPUQUEUE -gpu num=1:mode=exclusive_process -n 1 -W $RUNTIME -o logs/runFEP_%J.out scripts/BSSrunFEP.sh ${ligpair[0]} ${ligpair[1]} ${ligpair[3]} ${ligpair[4]}"
    echo $ID3

    # Stage 2c. Single job to process output of free and bound legs
    ID4=$(nk_jobid bsub -J analyseFEP -w "done($ID3)" -ti -q $CPUQUEUE -n 1 -W $ANALYSISTIME -o logs/analyseFEP_%J.out scripts/BSSanalyseFEP.sh ${ligpair[0]} ${ligpair[1]} ${ligpair[4]})

    echo "bsub -J analyseFEP -w done($ID3) -ti -q $CPUQUEUE -n 1 -W $ANALYSISTIME -o logs/analyseFEP_%J.out scripts/BSSanalyseFEP.sh ${ligpair[0]} ${ligpair[1]} ${ligpair[4]}"
    echo $ID4

    # sleep 1 second to give lag between process IDs
    sleep 1
done
