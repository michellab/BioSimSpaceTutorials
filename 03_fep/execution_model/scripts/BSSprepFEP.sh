#!/usr/bin/env bash
echo $SLURM_JOB_NAME
date

#tmp for local testing.
BSSHOME="/home/jscheen"
TMPDIR="./tmp/"
$BSSHOME/biosimspace.app/bin/ipython scripts/BSSprepFEP.py ejm_50 ejm_31

#$BSSHOME/sire.app/bin/python scripts/BSSprepFEP.py $1 $2
sleep 5
exit 0
