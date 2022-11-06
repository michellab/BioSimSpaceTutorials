#!/bin/bash
#SBATCH -o somd-array-gpu-%A.%a.out
#SBATCH -p serial
#SBATCH -n 1
#SBATCH --time 24:00:00
#SBATCH --array=0-7

lam=${lamvals[SLURM_ARRAY_TASK_ID]}

echo $((1*SLURM_ARRAY_TASK_ID))
sleep $((1*SLURM_ARRAY_TASK_ID))


if [ "$SLURM_ARRAY_TASK_ID" -eq "7" ]
then
  echo ID is 7
fi

