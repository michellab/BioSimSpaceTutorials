#!/bin/bash
#SBATCH --job-name=steeredMD
#SBATCH --gres=gpu:1
#SBATCH --time=80:00:00
#SBATCH -o steered_MD.out

#load modules
#ensure AMBERHOME is set
source /home/adele/software/amber20_20.04/amber20/amber.sh

#get the right python environment
source /home/adele/anaconda3/conda.sh
conda activate BSS-env

#run python script
#output will be saved in current directory
python ~/Documents/BioSimSpaceTutorials/02_steered_md/01_run_sMD.py --topology ~/Documents/BioSimSpaceTutorials/02_steered_md/data/system.prm7 --coordinates ~/Documents/BioSimSpaceTutorials/02_steered_md/data/system.rst7 --reference ~/Documents/BioSimSpaceTutorials/02_steered_md/data/reference.pdb --residues 178-184 --steering_runtime 0.2 --total_runtime 0.25 --force 3500

