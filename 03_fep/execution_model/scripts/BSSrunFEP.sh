#!/usr/bin/env bash
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
date
#$BSSHOME/sire.app/bin/python scripts/BSSrunFEP.py $1 $2 $3 $SLURM_ARRAY_TASK_ID

lig0=$1
lig1=$2
lambdastring=$3
engine=$4

IFS=',' read -r -a lambdas <<< "$lambdastring"

echo "${lambdas[@]}"
nwin=${#lambdas[@]}

SLURM_ARRAY_TASK_ID=$5

echo $SLURM_ARRAY_TASK_ID
echo $nwin


if [ $SLURM_ARRAY_TASK_ID -ge $nwin ]
then
    stage="free"
    idx=$(( $SLURM_ARRAY_TASK_ID - $nwin ))
else
    stage="bound"
    idx=$SLURM_ARRAY_TASK_ID 
fi

echo "idx is " $idx
lambda=${lambdas[$idx]}
echo "lambda is " $lambda

if [ $engine = "SOMD" ]
then
    echo "USING SOMD"
    #cd run/$lig0~$lig1/$stage/$lambda/
    echo "cd run/$lig0~$lig1/$stage/$lambda/"
    #$BSSHOME/sire.app/bin/somd-freenrg -C sim.cfg -l $lambda 1> somd.log 2> somd.err
    echo "$BSSHOME/sire.app/bin/somd-freenrg -C sim.cfg -l $lambda 1> somd.log 2> somd.err"
elif [ $engine = "GROMACS" ]
then
    echo "USING GROMACS"
else
    echo "The FEP engine $engine is not supported"
fi

sleep 5
exit 0
