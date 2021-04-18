import BioSimSpace as BSS
import sys
# sys.argv[1] is ligand0
# sys.argv[2] is ligand1
# verify that the result has not already been generated. if so exit. 
# verify that the run folder for ligand0-ligand1 exists
# verify that output files are available for the chosen engine
# if error complain and abort. 
# compute free energy change for the bound leg
# compute free energy change for the free leg
# compute the relative binding free energy
# append the binding energy to a file the output folder
print ("%s %s %s" % (sys.argv[0], sys.argv[1], sys.argv[2]))
