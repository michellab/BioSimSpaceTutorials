#!/usr/bin/env bash

## tmp for local testing
#SLURM_JOB_NAME=("ligprep")

echo $SLURM_JOB_NAME
date
$BSSHOME/bin/ipython scripts/BSSligprep.py $SLURM_ARRAY_TASK_ID

date

exit 0
