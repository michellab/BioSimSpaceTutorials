import BioSimSpace as BSS
import sys
#sys.argv[1] is ligand0
#sys.argv[2] is ligand1
# Load equilibrated free input for ligand0. Complain if input not found
# Load equilibrated free input for ligand1. Complain if input not found
# Extract ligand1 from frame.
# Align on ligand0
# Generate merge molecule.
# Write output files for perturbation using engine specified in 'protocol.dat' (so prm7/rst7 and pert if SOMD)

# Repeat same procedure for equilibrated bound input
print ("%s %s %s" % (sys.argv[0],sys.argv[1],sys.argv[2]))
