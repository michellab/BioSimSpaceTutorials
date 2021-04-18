#!/usr/bin/env bash
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
date
$BSSHOME/sire.app/bin/python scripts/BSSrunFEP.py $1 $2 $3 $SLURM_ARRAY_TASK_ID
sleep 5
exit 0
