#!/usr/bin/env bash
echo $SLURM_JOB_NAME
date
/home/users/common/BioSimSpace/sire.app/bin/python scripts/BSSprepFEP.py $1 $2
sleep 5
exit 0
