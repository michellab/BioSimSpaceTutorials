#!/usr/bin/env bash
# = = = 
# Preamble
# We assume interactive execution of notebook generates
# ligands.dat
# network.dat
# protocol.dat
# a parameterised protein input stored in inputs/protein/protein.prm7[rst7]
# a list of ligand mol2 files stored in inputs/ligands/
# The main output will be stored in outputs/free_energies.dat 
# = = =
# Stage 0 set cluster specific parameters
CPUQUEUE="serial"
GPUQUEUE="GTX980"
export BSSHOME="/export/users/julien/miniconda3/"
export TMPDIR="./tmp/"
export OPENMM_PLUGIN_DIR="$BSSHOME/lib/plugins/"
export OMP_NUM_THREADS=1 # avoid oversubscribing CPUs with SOMD
module load cuda/10.1 # for execution on OpenMM's CUDA platform
module load amber/20 # for preparing ligands and protein
module load gromacs/20.4 # for solvating

# Max run time for various categories of jobs. Should work for most use cases. 
PARAMTIME="01:00:00"
PREPTIME="01:00:00"
RUNTIME="04:00:00"
ANALYSISTIME="00:30:00"


# = = = 
# Jobs executed on the cluster

# Stage 1 - A GPU job array executed for each ligand listed in ligands.dat
num=$(< "ligands.dat" wc -l)
tasks=$(( $num -1 ))
echo "@@@ Parameterising dataset @@@"
#echo $tasks
ID1=$(sbatch --parsable --array [0-$tasks] --partition=$GPUQUEUE --ntasks=1 --gres=gpu:1 --time=$PARAMTIME --job-name=ligprep --output=logs/ligprep_%A_%a.out scripts/BSSligprep.sh)
echo "sbatch --parsable --array [0-$tasks] --partition=$GPUQUEUE --ntasks=1 --gres=gpu:1 --time=$PARAMTIME --job-name=ligprep --output=logs/ligprep_%A_%a.out scripts/BSSligprep.sh"
#echo $ID1

# Stage 2 - A series of jobs are executed once all ligands have been prepared
mapfile PERTS < network.dat
for pert in "${PERTS[@]}"
do
    IFS=' ' read -a ligpair <<< "$pert"
    win=$(( ${ligpair[2]} ))
    echo "@@@ Processing ${ligpair[0]} ${ligpair[1]} $win @@@" 

    # Stage 2a. Generate FEP inputs with a single job over a single CPU
    ID2=$(sbatch --parsable --dependency=afterany:${ID1} --partition=$CPUQUEUE --ntasks=1 --time=$PREPTIME --job-name=prepFEP --output=logs/prepFEP_%A.out scripts/BSSprepFEP.sh ${ligpair[0]} ${ligpair[1]})
    echo "sbatch --parsable --dependency=afterany:${ID1}  --partition=$CPUQUEUE --ntasks=1 --time=$PREPTIME --job-name=prepFEP --output=logs/prepFEP_%A.out scripts/BSSprepFEP.sh ${ligpair[0]} ${ligpair[1]}"
    echo $ID2

    # Stage 2b. Job array with 2*nwindows to run over free and bound legs over GPUs
    tasks=$(( 2*${ligpair[2]} - 1 ))
    ID3=$(sbatch --parsable --dependency=afterany:${ID2} --partition=$GPUQUEUE --ntasks=1 --gres=gpu:1 --time=$RUNTIME --job-name=runFEP --output=logs/runFEP_%A_%a.out --array=[0-$tasks] scripts/BSSrunFEP.sh ${ligpair[0]} ${ligpair[1]} ${ligpair[3]} ${ligpair[4]})
    echo "sbatch --parsable --dependency=afterany:${ID2} --partition=$GPUQUEUE --ntasks=1 --gres=gpu:1 --time=$RUNTIME --job-name=runFEP --output=logs/runFEP_%A_%a.out --array=[0-$tasks] scripts/BSSrunFEP.sh ${ligpair[0]} ${ligpair[1]} ${ligpair[3]} ${ligpair[4]}"
    echo $ID3

    # Stage 2c. Single job to process output of free and bound legs 
    ID4=$(sbatch --parsable --dependency=afterany:${ID3} --partition=$CPUQUEUE --ntasks=1 --time=$ANALYSISTIME --job-name=analyseFEP --output=logs/analyseFEP_%A.out scripts/BSSanalyseFEP.sh ${ligpair[0]} ${ligpair[1]} ${ligpair[4]})
    echo "sbatch --parsable --dependency=afterany:${ID3} --partition=$CPUQUEUE --ntasks=1 --time=$ANALYSISTIME --job-name=analyseFEP --output=logs/analyseFEP_%A.out scripts/BSSanalyseFEP.sh ${ligpair[0]} ${ligpair[1]} ${ligpair[4]}"
    echo $ID4

    # sleep 1 second to give lag between process IOs
    sleep 1
done
