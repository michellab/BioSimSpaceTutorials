#!/bin/bash
#SBATCH -o somd-array-gpu-%A.%a.out
#SBATCH -p GTX1080
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --time 24:00:00
#SBATCH --array=0-17

sleep $((30*SLURM_ARRAY_TASK_ID))

module load cuda

echo "CUDA DEVICES:" $CUDA_VISIBLE_DEVICES

lamvals=( 0.000 0.028 0.056 0.111 0.167 0.222 0.278 0.333 0.389 0.444 0.500 0.556 0.611 0.667 0.722 0.778 0.889 1.000 )
lam=${lamvals[SLURM_ARRAY_TASK_ID]}

echo "lambda is: " $lam

mkdir lambda-$lam
cd lambda-$lam

export OPENMM_PLUGIN_DIR=/export/users/finlayclark/biosimspace.app/lib/plugins/

srun ~/biosimspace.app/bin/somd-freenrg -C ../../input/sim.cfg -l $lam -p CUDA
cd ..

if [ "$SLURM_ARRAY_TASK_ID" -eq "17" ]
then
  wait
  sleep 30
  sbatch ../ljcor.sh
  sleep 30 
  sbatch ../mbar.sh
fi

