import BioSimSpace as BSS
from BioSimSpace import _Exceptions
import sys
import csv


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
at this point BSS should have thrown an error if parameterisation failed, so no checks needed.


############### TMP LOAD pre-param ligand#####
#BSS.IO.saveMolecules("inputs/tmp", lig_p, ["PRM7", "RST7"])
lig_p = BSS.IO.readMolecules(["inputs/tmp.prm7", "inputs/tmp.rst7"])[0]


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


"""
figure out the ligand dimensions, add the specified water box together with the largest ligand axis.
This is a workaround to account for adding ions in some cases. Based on:
https://github.com/michellab/BioSimSpace/issues/111
"""
box_min, box_max = lig_p.getAxisAlignedBoundingBox()
box_size = [y - x for x, y in zip(box_min,box_max)]
box_sizes = [x + int(box_axis_length) * box_axis_unit for x in box_size]

# do the same for ligand +protein system.
protein = BSS.IO.readMolecules(["inputs/protein/protein.rst7", "inputs/protein/protein.prm7"])[0]
system = lig_p + protein
box_min_s, box_max_s = system.getAxisAlignedBoundingBox()
box_size_s = [y - x for x, y in zip(box_min_s,box_max_s)]
box_sizes_s = [x + int(box_axis_length) * box_axis_unit for x in box_size_s]


# get the box type settings.
boxtype_query = lines[4].rstrip().replace(" ","").split("=")[-1]
if boxtype_query.lower() == "orthorhombic":
	print("Solvating ligand.")
	box, angles = BSS.Box.cubic(max(box_sizes))
	lig_p_solvated = BSS.Solvent.solvate(solvent_query, molecule=lig_p,
                               box=box, angles=angles)

	print("Solvating ligand + protein.")
	box, angles = BSS.Box.cubic(max(box_sizes_s))
	system_solvated = BSS.Solvent.solvate(solvent_query, molecule=system,
	                               box=box, angles=angles)



elif boxtype_query.lower() == "triclinic" or boxtype_query.lower() == "octahedral":
	print("Solvating ligand.")
	box, angles = BSS.Box.truncatedOctahedron(max(box_sizes))
	lig_p_solvated = BSS.Solvent.solvate(solvent_query, molecule=lig_p,
                           box=box, angles=angles)

	print("Solvating ligand + protein.")
	box, angles = BSS.Box.truncatedOctahedron(max(box_sizes_s))
	system_solvated = BSS.Solvent.solvate(solvent_query, molecule=system,
	                               box=box, angles=angles)

else:
	raise NameError("Input box type not recognised. Please use any of ['orthorhombic', 'octahedral', 'triclinic']" \
	+"in the fifth line of protocol.dat in the shape of (e.g.):\nbox type = orthorhombic")

#################
### Minimise and equilibrate our systems to make the ligand relax inside the pocket.

# first save and reload the systems to bypass some amber checks that block the next few steps.
BSS.IO.saveMolecules("inputs/lig_tmp", lig_p_solvated, ["PRM7", "RST7"])
lig_p_solvated = BSS.IO.readMolecules(["inputs/lig_tmp.prm7", "inputs/lig_tmp.rst7"])

def MinimEquil(system, minim_steps=250, runtime_nvt=50, runtime_npt=300):
	"""
	Given a solvated system (BSS object), run a minimisation and restrained equilibration workflow.
	Returns the equilibrated system.
	"""
	#####
	print(f"Minimising in {minim_steps} steps..")
	protocol = BSS.Protocol.Minimisation(steps=minim_steps)
	minimised = BSS.Process.Amber(system, protocol).start().getSystem(block=True)

	#####
	print(f"NVT equilibration for {runtime_nvt} ps while restraining all non-solvent atoms..")
	protocol = BSS.Protocol.Equilibration(
					runtime=runtime_nvt*BSS.Units.Time.picosecond, 
					temperature_start=0*BSS.Units.Temperature.kelvin, 
					temperature_end=300*BSS.Units.Temperature.kelvin,
					restraint="all")

	equil1 = BSS.Process.Amber(minimised, protocol).start().getSystem(block=True)

	#####
	print(f"NVT equilibration for {runtime_nvt} ps while restraining backbone atoms and ligand..")
	protocol = BSS.Protocol.Equilibration(
					runtime=runtime_nvt*BSS.Units.Time.picosecond, 
					temperature_start=0*BSS.Units.Temperature.kelvin, 
					temperature_end=300*BSS.Units.Temperature.kelvin,
					restraint="backbone")

	equil2 = BSS.Process.Amber(equil1, protocol).start().getSystem(block=True)

	#####
	print(f"NVT equilibration for {runtime_nvt} ps without restraints..")
	protocol = BSS.Protocol.Equilibration(
					runtime=runtime_nvt*BSS.Units.Time.picosecond, 
					temperature_start=0*BSS.Units.Temperature.kelvin, 
					temperature_end=300*BSS.Units.Temperature.kelvin,
					)

	equil3 = BSS.Process.Amber(equil2, protocol).start().getSystem(block=True)

	#####
	print(f"NPT equilibration for {runtime_npt} ps while restraining non-solvent heavy atoms..")
	protocol = BSS.Protocol.Equilibration(
					runtime=runtime_npt*BSS.Units.Time.picosecond, 
					pressure=1*BSS.Units.Pressure.atm,
					temperature=300*BSS.Units.Temperature.kelvin,
					restraint="heavy")
	equil4 = BSS.Process.Amber(equil3, protocol).start().getSystem(block=True)

	#####
	print(f"NPT equilibration for {runtime_npt} ps without restraints..")
	protocol = BSS.Protocol.Equilibration(
					runtime=runtime_npt*BSS.Units.Time.picosecond, 
					pressure=1*BSS.Units.Pressure.atm,
					temperature=300*BSS.Units.Temperature.kelvin,
					)
	equil_fin = BSS.Process.Amber(equil4, protocol).start().getSystem(block=True)

	return equil_fin

lig_s_equil = MinimEquil(lig_p_solvated, minim_steps=250, runtime_nvt=50, runtime_npt=300)

# free input
# X steps energy minimise with BSS.Protocol.Minimize
# X ps NVT equilibration with restraints on all non solvent atoms		# "all"
# X ps NVT equilibration with restraints on backbone atoms and ligand 	# "backbone"
# X ps NVT equilibration unrestrained
# X ps NPT equilibration restraints on non solvent atoms				# "heavy"
# X ps NPT equilibration unrestrained

# bound input
# X steps energy minimise with BSS.Protocol.Minimize					# "all"
# X ps NVT equilibration with restraints on all non solvent atoms		# "backbone"
# X ps NVT equilibration with restraints on backbone atoms and ligand
# X ps NVT equilibration unrestrained
# X ps NPT equilibration restraints on non solvent atoms				# "heavy"
# X ps NPT equilibration unrestrained


# Save last snapshot of free equilibration
# Save last snapshot of bound equilibration






























