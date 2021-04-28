#!/usr/bin/env bash
echo $SLURM_JOB_NAME
date

$BSSHOME/biosimspace.app/bin/ipython scripts/BSSprepFEP.py $1 $2
date
sleep 5
exit 0
