#!/usr/bin/env bash

## tmp for local testing
SLURM_JOB_NAME=("ligprep")

TMPDIR="./tmp/"

export OMP_NUM_THREADS=1 # avoid oversubscribing CPUs when using SOMD
source /export/users/marie/gromacs-2020.4/install/bin/GMXRC # for solvating

echo $SLURM_JOB_NAME
date
$BSSHOME/biosimspace.app/bin/ipython scripts/BSSligprep.py $SLURM_ARRAY_TASK_ID

date

exit 0
