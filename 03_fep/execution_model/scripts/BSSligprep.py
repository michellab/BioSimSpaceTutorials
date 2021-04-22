import BioSimSpace as BSS
from BioSimSpace import _Exceptions
import sys
import csv


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

### preamble.
print (f"{sys.argv[0]} {sys.argv[1]}")
idx = int( sys.argv[1] )

#################
### Open 'ligands.dat', find sys.argv[1] 'th entry 
stream = open("ligands.dat","r")
lines = stream.readlines()
lig_name = lines[idx].rstrip()

#################
### Load matching input with BSS.read.IO.
lig = BSS.IO.readMolecules(f"inputs/ligands/{lig_name}.mol2")[0]

#################
### parameterise ligand by deriving requested FF from protocol.dat.
stream = open("protocol.dat","r")
lines = stream.readlines()
ff_query = lines[0].rstrip()

# do tests on protocol, parameterise with correct FF.
if not "ligand" in ff_query:
	raise NameError("Please supply a ligand force field on the first line of protocol.dat. The line should look like (e.g.):\n"+\
		"ligand forcefield = GAFF2")
else:
	if "GAFF1" in ff_query:
		print(f"Parameterising {lig_name} using GAFF1 force field.")
		lig_p = BSS.Parameters.gaff(lig).getMolecule()

	elif "GAFF2" in ff_query:
		print(f"Parameterising {lig_name} using GAFF2 force field.")
		lig_p = BSS.Parameters.gaff2(lig).getMolecule()

	else:
		raise NameError(f"Force field not supported: {ff_query}. Please use either of [GAFF1, GAFF2], or" \
							+" edit this script to use other force fields available in BSS.Parameters.forceFields().")
# at this point BSS should have thrown an error if parameterisation failed, so no checks needed.

#################
### make a copy of ligand, solvate original ligand object.
lig_p_copy = lig_p.copy()

# need to derive some settings from protocol.dat again. 
stream = open("protocol.dat","r")
lines = stream.readlines()

# get the solvent force field.
solvent_query = lines[2].rstrip().replace(" ","").split("=")[-1]

# get the box size settings.
boxsize_query = lines[3].rstrip().replace(" ","").split("=")[-1]
box_axis_length = boxsize_query.split("*")[0]
box_axis_unit = boxsize_query.split("*")[1]
if box_axis_unit.lower() == "nm" or box_axis_unit.lower() == "nanometer":
	box_axis_unit = BSS.Units.Length.nanometer

elif box_axis_unit.lower() == "a" or box_axis_unit.lower() == "angstrom":
	box_axis_unit = BSS.Units.Length.angstrom
else:
	raise NameError("Input unit not recognised. Please use any of ['spc', 'spce', 'tip3p', 'tip4p', 'tip5p'] in " \
	+"the fourth line of protocol.dat in the shape of (e.g.):\nbox edges = 10*angstrom")

# get the box type settings.
boxtype_query = lines[4].rstrip().replace(" ","").split("=")[-1]
if boxtype_query.lower() == "orthorhombic":
	box, angles = BSS.Box.cubic(int(box_axis_length) * BSS.Units.Length.angstrom)
elif boxtype_query.lower() == "triclinic" or boxtype_query.lower() == "octahedral":
	box, angles = BSS.Box.truncatedOctahedron(int(box_axis_length) * BSS.Units.Length.angstrom)
else:
	raise NameError("Input box type not recognised. Please use any of ['orthorhombic', 'octahedral', 'triclinic']" \
	+"the fifth line of protocol.dat in the shape of (e.g.):\nbox type = orthorhombic")



lig_p_solvated = BSS.Solvent.solvate(solvent_query, molecule=lig_p_copy,
                               box=box, angles=angles)

#################
### combine ligand copy with protein

#################
### solvate complex
































