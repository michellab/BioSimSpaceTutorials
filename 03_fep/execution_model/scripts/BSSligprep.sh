#!/usr/bin/env bash

## tmp for local testing
SLURM_JOB_NAME=("ligprep")
SLURM_ARRAY_TASK_ID=0
BSSHOME="/home/jscheen"
# change biosimspace.app back to sire.app?
### /

echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
date
$BSSHOME/biosimspace.app/bin/ipython scripts/BSSligprep.py $SLURM_ARRAY_TASK_ID

sleep 5
exit 0
