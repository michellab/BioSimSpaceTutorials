#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --job-name=prod
#SBATCH -o prod_%A_%a.out
#SBATCH -e prod_%A_%a.err

date

# this was run via command line using (remove spaces in the first word):
# s b a t c h --array=0-8 run_production_slurm.sh

# ensure that the path to SOMD can be found. This could be via the conda install for example,
# which could be in the following format depending on the path:
export BSS="/home/anna/anaconda3/bin/activate biosimspace-dev"
source $BSS

lamvals=( 0.0000 0.1250 0.2500 0.3750 0.5000 0.6250 0.7500 0.8750 1.0000 )

# define lambda based on slurm array
lam=${lamvals[SLURM_ARRAY_TASK_ID]}

trans_dir=$(pwd)

# iterate over dir (for each leg) based on no of repeats
for dir in 'bound' 'free'; do

cd $dir

echo "Lambda is $lam"

cd lambda_$lam
somd-freenrg -c somd.rst7 -t somd.prm7 -m somd.pert -C somd.cfg -p CUDA

cd $trans_dir

done

echo "done."
