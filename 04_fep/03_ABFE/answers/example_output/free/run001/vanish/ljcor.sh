#!/bin/bash

cd lambda-0.000
/home/finlayclark/anaconda3/envs/biosimspace-dev/bin/lj-tailcorrection -C ../../input/sim.cfg -l 0.00 -r traj000000001.dcd -s 1 1> ../freenrg-LJCOR-lam-0.000.dat 2> /dev/null

cd ..
cd lambda-1.000
/home/finlayclark/anaconda3/envs/biosimspace-dev/bin/lj-tailcorrection -C ../../input/sim.cfg -l 1.00 -r traj000000001.dcd -s 1 1> ../freenrg-LJCOR-lam-1.000.dat 2> /dev/null
cd ..

wait

# utility script to get final LJ correction term
python ~/Documents/research/scripts/abfe/parselj.py freenrg-LJCOR-lam-0.000.dat freenrg-LJCOR-lam-1.000.dat > freenrg-LJCOR.dat
