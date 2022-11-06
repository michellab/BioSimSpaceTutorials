#!/bin/bash
#SBATCH -o analyse-free-nrg-%A.%a.out
#SBATCH -p GTX1080
#SBATCH -n 1
#SBATCH --time 01:00:00

sleep 30 

export OPENMM_PLUGIN_DIR=/export/users/finlayclark/biosimspace.app/lib/plugins/

srun /export/users/finlayclark/biosimspace.app/bin/analyse_freenrg mbar -i lambda*/simfile.dat -p 83 --overlap --temperature 298.0 > freenrg-MBAR-5ns-overlap.dat

#bzip2 lambda-*/*




