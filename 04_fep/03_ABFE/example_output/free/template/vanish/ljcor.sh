#!/bin/bash

cd lambda-0.000
~/biosimspace.app/bin/lj-tailcorrection -C ../../input/sim.cfg -l 0.00 -r traj000000001.dcd -s 1 1> ../freenrg-LJCOR-lam-0.000.dat 2> /dev/null

cd ..
cd lambda-1.000
~/biosimspace.app/bin/lj-tailcorrection -C ../../input/sim.cfg -l 1.00 -r traj000000001.dcd -s 1 1> ../freenrg-LJCOR-lam-1.000.dat 2> /dev/null
cd ..

wait

# utility script to get final LJ correction term
/usr/bin/python3 /export/users/finlayclark/restraint_comparison_mif/study/scripts freenrg-LJCOR-lam-0.000.dat freenrg-LJCOR-lam-1.000.dat > freenrg-LJCOR.dat
