#!/bin/bash
#SBATCH -o boresch_analytical_correction-%A.%a.out
#SBATCH -p GTX1080
#SBATCH -n 1
#SBATCH --time 24:00:00

sleep 30 

srun /export/users/finlayclark/sire.app/pkgs/sire-2021.1.0/bin/boresch_analytical_correction -C ../input/sim.cfg --verbose > boresch_analytical_correction.dat

cd ..

wait


