#!/usr/bin/env bash
date

export OMP_NUM_THREADS=1 # avoid oversubscribing CPUs when using SOMD
export OPENMM_PLUGIN_DIR=/export/users/julien/sire.app/lib/plugins/

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
    
    echo "$BSSHOME/biosimspace.app/bin/somd-freenrg -C ./somd.cfg -l $lambda -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert -p CUDA 1> somd.log 2> somd.err"
    $BSSHOME/biosimspace.app/bin/somd-freenrg -C ./somd.cfg -l $lambda -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert -p CUDA 1> somd.log 2> somd.err
elif [[ $engine == *"GROMACS"* ]]; then
    echo "USING GROMACS"
    # here run a gromacs exe on a lambda folder as above.
else
    echo "The FEP engine $engine is not supported"
fi

wait
exit 0
