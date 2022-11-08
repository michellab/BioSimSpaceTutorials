#!/bin/bash
#SBATCH -o somd-array-gpu-%A.%a.out
#SBATCH -p main
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --time 24:00:00

echo "CUDA DEVICES:" $CUDA_VISIBLE_DEVICES


srun ~/sire_boresch_2022.app/bin/somd-freenrg -C ../../input/sim.cfg -l 0.250 -p CUDA

