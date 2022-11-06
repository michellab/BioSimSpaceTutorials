#!/bin/bash
#SBATCH -o somd-array-gpu-%A.%a.out --exclude=node02
#SBATCH -p GTX980
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --time 24:00:00
#SBATCH --array=0-35

sleep $((30*SLURM_ARRAY_TASK_ID))

module load cuda

echo "CUDA DEVICES:" $CUDA_VISIBLE_DEVICES

lamvals=( 0.000 0.025 0.050 0.075 0.100 0.125 0.150 0.175 0.200 0.225 0.250 0.275 0.300 0.325 0.350 0.375 0.400 0.425 0.450 0.475 0.500 0.525 0.550 0.575 0.600 0.625 0.650 0.675 0.700 0.725 0.750 0.800 0.850 0.900 0.950 1.000 )
lam=${lamvals[SLURM_ARRAY_TASK_ID]}

echo "lambda is: " $lam

mkdir lambda-$lam
cd lambda-$lam

#export OPENMM_PLUGIN_DIR=/export/users/finlayclark/biosimspace.app/lib/plugins/

srun /export/users/finlayclark/sire.app/bin/somd-freenrg -C ../../input/sim.cfg -l $lam -p CUDA
cd ..

if [ "$SLURM_ARRAY_TASK_ID" -eq "35" ]
then
  wait
  sleep 30
  sbatch ../ljcor.sh
  sleep 30
  sbatch ../mbar.sh
  sleep 30
  sbatch ../boresch_ana_corr.sh
fi

