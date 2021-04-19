#!/bin/bash
#SBATCH --job-name=steeredMD
#SBATCH --gres=gpu:1
#SBATCH --time=80:00:00
#SBATCH -o steered_MD.out

#load modules
#ensure AMBERHOME is set
source /home/adele/software/amber20_20.04/amber20/amber.sh

#get the right python environment

#run python script
#output will be saved in current directory
~/sire.app/bin/python ~/Documents/BioSimSpaceTutorials/02_steered_md/01_run_sMD.py --topology ~/Documents/BioSimSpaceTutorials/02_steered_md/data/system.prm7 --coordinates ~/Documents/BioSimSpaceTutorials/02_steered_md/data/system.rst7 --reference ~/Documents/BioSimSpaceTutorials/02_steered_md/data/reference.pdb --steering_runtime 1 --total_runtime 1.2 --force 3500

