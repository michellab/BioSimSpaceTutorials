#!/bin/bash
-rw-rw-r--  1 finlayclark finlayclark 2.2K Sep  3 12:48 gradients.dat
-rw-rw-r--  1 finlayclark finlayclark 1.1K Sep  3 12:48 gradients.s3
-rw-rw-r--  1 finlayclark finlayclark 4.7M Sep  3 12:48 latest.pdb
-rw-rw-r--  1 finlayclark finlayclark  144 Sep  3 12:48 moves.dat
-rw-rw-r--  1 finlayclark finlayclark  21M Sep  3 12:48 simfile.dat
-rw-rw-r--  1 finlayclark finlayclark  32M Sep  3 12:49 sim_restart.s3
-rw-rw-r--  1 finlayclark finlayclark  32M Sep  3 12:48 sim_restart.s3.previous
-rw-rw-r--  1 finlayclark finlayclark 8.8M Sep  3 12:48 SYSTEM.s3

#bound
for dir in bound/run00{1..5}/; 
do
	for leg in restrain discharge vanish
		do 
			pushd $dir/$leg/output;
			rm freenrg-*;
			rm somd*;
			rm boresch*;
			rm slurm*;
			rm analyse*;
			rm traj*;
			for lam in lam*;
			do
			       	pushd $lam;
				rm traj*;
				rm gradients.s3  				
				rm latest.pdb
				rm moves.dat
				rm sim_restart.s3*
				rm SYSTEM.s3
				popd;
				mv $lam ${lam//-/_} 
			done
			popd;
		done
done

#free
for dir in free/run00{1..5}/; 
do
	for leg in discharge vanish
		do 
			pushd $dir/$leg/output;
			rm freenrg-*;
			rm somd*;
			rm boresch*;
			rm slurm*;
			rm analyse*;
			rm traj*;
			for lam in lam*;
			do
			       	pushd $lam;
				rm traj*;
				rm gradients.s3  				
				rm latest.pdb
				rm moves.dat
				rm sim_restart.s3*
				rm SYSTEM.s3
				popd;
				mv $lam ${lam//-/_} 
			done
			popd;
		done
done
