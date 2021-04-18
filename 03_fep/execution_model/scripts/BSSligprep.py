import BioSimSpace as BSS
import sys
# Open 'ligands.dat', find sys.argv[1] 'th entry 
# Load matching input with BSS.read.IO
# parameterise ligand
# make a copy of ligand, solvate ligand'
# combine ligand with protein
# solvate complex
# free input
# X steps energy minimise with BSS.Protocol.Minimize
# X ps NVT equilibration with restraints on all non solvent atoms
# X ps NVT equilibration with restraints on backbone atoms and ligand
# X ps NVT equilibration unrestrained
# X ps NPT equilibration restraints on non solvent atoms
# X ps NPT equilibration unrestrained
# bound input
# X steps energy minimise with BSS.Protocol.Minimize
# X ps NVT equilibration with restraints on all non solvent atoms
# X ps NVT equilibration with restraints on backbone atoms and ligand
# X ps NVT equilibration unrestrained
# X ps NPT equilibration restraints on non solvent atoms
# X ps NPT equilibration unrestrained
# Save last snapshot of free equilibration
# Save last snapshot of bound equilibration

print ("%s %s" % (sys.argv[0],sys.argv[1]))
idx = int( sys.argv[1] )
stream = open("ligands.dat","r")
lines = stream.readlines()
lig = lines[idx]
print ("Parameterising %s" % lig)

