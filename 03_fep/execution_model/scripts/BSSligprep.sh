#!/usr/bin/env bash
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
date
/home/users/common/BioSimSpace/sire.app/bin/python scripts/BSSligprep.py $SLURM_ARRAY_TASK_ID

sleep 5
exit 0
