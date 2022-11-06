To process the entire FEP dataset,

if necessary edit the variables \$CPUQUEUE and ​\$GPUQUEUE near the top of the bash script. 



Then, when using a cluster configured with SLURM, type:

bash processFEP-slurm.sh

or when using a cluster configured with LSF, type:

bash processFEP-lsf.sh



On running this script some messages should be printed out during which jobs will be submitted. The order of dependent jobs is:

- ligprep
- prepFEP
- runFEP
- analyseFEP



After job completion the main results can be found in ```./outputs/summary.csv```.

