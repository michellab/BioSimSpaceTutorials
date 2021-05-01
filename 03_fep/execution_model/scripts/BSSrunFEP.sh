#!/usr/bin/env bash
date

lig0=$1
lig1=$2
lambdastring=$3
engine=$4

IFS=',' read -r -a lambdas <<< "$lambdastring"

nwin=${#lambdas[@]}


if [ "$SLURM_ARRAY_TASK_ID" -ge "$nwin" ]; then
    stage="free"
    idx=$(( $SLURM_ARRAY_TASK_ID - $nwin ))
else
    stage="bound"
    idx=$SLURM_ARRAY_TASK_ID 
fi


echo "idx is " $idx
lambda=${lambdas[$idx]}
echo "lambda is: "$lambda
echo "stage is: "$stage
echo "engine is: "$engine


if [[ $engine == *"SOMD"* ]]; then
    echo "USING SOMD"
    echo "cd outputs/SOMD/$lig0~$lig1/$stage/lambda_$lambda/"
    cd outputs/SOMD/$lig0~$lig1/$stage/lambda_$lambda/
    
    # for SOMD, append a few settings to the config file. These are up-to-date recommended settings.
    echo "hydrogen mass repartitioning factor = 1.0" >> somd.cfg
    echo "cutoff type = cutoffperiodic" >> somd.cfg
    echo "cutoff distance = 10*angstrom" >> somd.cfg
    echo "barostat = True" >> somd.cfg
    echo "andersen = True" >> somd.cfg
    echo "precision = mixed" >> somd.cfg
    echo "center solute = True" >> somd.cfg
    echo "reaction field dielectric = 82.0" >> somd.cfg
    echo "minimal coordinate saving = False" >> somd.cfg    

    # run SOMD simulation.    
    echo "$BSSHOME/bin/somd-freenrg -C ./somd.cfg -l $lambda -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert -p CUDA 1> somd.log 2> somd.err"
    $BSSHOME/bin/somd-freenrg -C ./somd.cfg -l $lambda -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert -p CUDA 1> somd.log 2> somd.err

elif [[ $engine == *"GROMACS"* ]]; then
    echo "USING GROMACS"
    # here run a gromacs exe on a lambda folder as above.
else
    echo "The FEP engine $engine is not supported"
fi

exit 0
