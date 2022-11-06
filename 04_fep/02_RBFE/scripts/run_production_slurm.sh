#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --job-name=prod
#SBATCH -o ../slurm_logs/prod_%A_%a.out
#SBATCH -e ../slurm_logs/prod_%A_%a.err

# sourcing
module load cuda/11.6
somdfreenrg="/export/users/XXX/sire.app/bin/somd-freenrg"
source $BSS
source extract_execution_model_bash.sh

# keep_traj="None"

date
echo "Folder for these runs is : $MAINDIRECTORY"
echo "The transformation is $1 using $2 windows and $3 as the MD engine"

# define no of windows based on sys arg 2
if [ $3 = 11 ]; then
lamvals=( 0.0000 0.1000 0.2000 0.3000 0.4000 0.5000 0.6000 0.7000 0.8000 0.9000 1.0000 )
fi
if [ $3 = 17 ]; then
lamvals=( 0.0000 0.0625 0.1250 0.1875 0.2500 0.3125 0.3750 0.4375 0.5000 0.5625 0.6250 0.6875 0.7500 0.8125 0.8750 0.9375 1.0000 )
fi

# define lambda based on slurm array
lam=${lamvals[SLURM_ARRAY_TASK_ID]}

# change to the trans dir, abort and message if not there
cd $MAINDIRECTORY/outputs/$2/$1
if [[ ! -d $MAINDIRECTORY/outputs/$2/$1 ]]; then
    echo "$MAINDIRECTORY/outputs/$2/$1 does not exist. Production run aborted..."
    exit
fi
trans_dir=$(pwd)

for dir in 'bound' 'free'; do
cd $dir
cd lambda_$lam
$somdfreenrg -c somd.rst7 -t somd.prm7 -m somd.pert -C somd.cfg -p CUDA
cd $repeat_dir

# delete simulation data 
# if [[ $keep_traj == "None" ]]; then
#     echo "Removing all trajectories..."
#     rm eq/lambda_$lam/*.dcd eq/lambda_$lam/*.s3* eq/lambda_$lam/latest*
#     rm lambda_$lam/*.dcd lambda_$lam/*.s3* lambda_$lam/latest*
# fi

cd $trans_dir

done

echo "done."
