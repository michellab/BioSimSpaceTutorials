import BioSimSpace as BSS
import sys
#sys.argv[1] is ligand0
#sys.argv[2] is ligand1
#sys.argv[3] is the number of windows to use per perturbation
#sys.argv[4] is the SLURM array index.
# if index < nwindows run the i'th window of the bound leg.  
# if index > nwindows run the i'th window of the free leg. 

# the protocol and simulation engine are specified in 'protocol.dat'
# the input should be compatible with the engine specified in protocol.dat

# verify that a trajectory has not been already generated.
# if trajectory complete exit.
# if trajectory partial resume job. 
# else start from scratch.
print ("%s %s %s %s %s" % (sys.argv[0], sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]))
