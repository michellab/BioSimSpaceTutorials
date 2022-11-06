#!/bin/bash
#SBATCH -o analyse-free-nrg-%A.%a.out
#SBATCH -p GTX1080
#SBATCH -n 1
#SBATCH --time 01:00:00

sleep 30 

export OPENMM_PLUGIN_DIR=/home/users/common/BioSimSpace/sire.app/lib/plugins

srun /home/users/common/BioSimSpace/sire.app/bin/analyse_freenrg mbar -i lambda*/simfile.dat -p 83 --overlap --temperature 298.0 > freenrg-MBAR-5ns-overlap.dat

#bzip2 lambda-*/*




