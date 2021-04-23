import BioSimSpace as BSS
import sys
import csv


print ("%s %s %s" % (sys.argv[0],sys.argv[1],sys.argv[2]))

# Load equilibrated free inputs for both ligands. Complain if input not found
print(f"Loading ligands {sys.argv[1]} and {sys.argv[2]}.")
ligs_path = "inputs/ligands/"
ligand_1 = BSS.IO.readMolecules([f"{ligs_path}{sys.argv[1]}_lig_equil_solv.rst7", f"{ligs_path}{sys.argv[1]}_lig_equil_solv.prm7"])
ligand_2 = BSS.IO.readMolecules([f"{ligs_path}{sys.argv[2]}_lig_equil_solv.rst7", f"{ligs_path}{sys.argv[2]}_lig_equil_solv.prm7"])

# Extract ligands.
ligand_1 = ligand_1.getMolecule(0)
ligand_2 = ligand_2.getMolecule(0)

# Align ligand1 on ligand2
print("Mapping and aligning..")
mapping = BSS.Align.matchAtoms(ligand_1, ligand_2, sanitize=True)
ligand_1_a = BSS.Align.rmsdAlign(ligand_1, ligand_2, mapping)

# Generate merged molecule.
print("Merging..")
merged_ligs = BSS.Align.merge(ligand_1_a, ligand_2, mapping)


## now repeat above steps, but for the protein + ligand systems.
# Load equilibrated bound inputs for both ligands. Complain if input not found
print(f"Loading bound ligands {sys.argv[1]} and {sys.argv[2]}.")
ligs_path = "inputs/ligands/"
system_1 = BSS.IO.readMolecules([f"{ligs_path}{sys.argv[1]}_sys_equil_solv.rst7", f"{ligs_path}{sys.argv[1]}_sys_equil_solv.prm7"])
system_2 = BSS.IO.readMolecules([f"{ligs_path}{sys.argv[2]}_sys_equil_solv.rst7", f"{ligs_path}{sys.argv[2]}_sys_equil_solv.prm7"])

# Extract ligands.
system_ligand_1 = system_1.getMolecule(0)
system_ligand_2 = system_2.getMolecule(0)

# extract protein for ligand 2.
protein = system_2.getMolecule(1)

# Align ligand1 on ligand2
print("Mapping and aligning..")
mapping = BSS.Align.matchAtoms(system_ligand_1, system_ligand_2, sanitize=True)
system_ligand_1_a = BSS.Align.rmsdAlign(system_ligand_1, system_ligand_2, mapping)

# Generate merged molecule.
print("Merging..")
system_merged_ligs = BSS.Align.merge(system_ligand_1_a, system_ligand_2, mapping)

# then add the protein back onto the merged object.
merged_system = protein + system_merged_ligs 


#################################### solvate merged molecule (same approach as BSSligprep.py).
########################
def solvateSystem(system):
    """
    Given a system (merged molecule or merged molecule + protein), solvate based on settings in protocol.dat.
    Returns solvated system.
    """
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

    box_min, box_max = system.getAxisAlignedBoundingBox()
    box_size = [y - x for x, y in zip(box_min,box_max)]
    box_sizes = [x + int(box_axis_length) * box_axis_unit for x in box_size]

    boxtype_query = lines[4].rstrip().replace(" ","").split("=")[-1]
    if boxtype_query.lower() == "orthorhombic":
        box, angles = BSS.Box.cubic(max(box_sizes))
        system_solvated = BSS.Solvent.solvate(solvent_query, molecule=system,
                                   box=box, angles=angles)


    elif boxtype_query.lower() == "triclinic" or boxtype_query.lower() == "octahedral":
        box, angles = BSS.Box.truncatedOctahedron(max(box_sizes))
        system_solvated = BSS.Solvent.solvate(solvent_query, molecule=system,
                               box=box, angles=angles)

    else:
        raise NameError("Input box type not recognised. Please use any of ['orthorhombic', 'octahedral', 'triclinic']" \
        +"in the fifth line of protocol.dat in the shape of (e.g.):\nbox type = orthorhombic")

    return system_solvated

print("Solvating merged ligands.")
merged_ligs_solvated = solvateSystem(merged_ligs)
print("Solvating merged ligands + protein.")
merged_system_solvated = solvateSystem(merged_system)

########################### now set up the SOMD or GROMACS MD directories. 
#first, figure out which engine and what runtime the user has specified in protocol.
stream = open("protocol.dat","r")
lines = stream.readlines()

### get the requested engine.
engine_query = lines[7].rstrip().replace(" ","").split("=")[-1].upper()
if engine_query not in ["SOMD", "GROMACS"]:
    raise NameError("Input MD engine not recognised. Please use any of ['SOMD', 'GROMACS']" \
    +"on the eighth line of protocol.dat in the shape of (e.g.):\nengine = SOMD")

### get the requested runtime.
runtime_query = lines[6].rstrip().replace(" ","").split("=")[-1].split("*")[0]
try:
    runtime_query = int(runtime_query)
except ValueError:
    raise NameError("Input runtime value not supported. Please use an integer" \
    +" on the seventh line of protocol.dat in the shape of (e.g.):\nsampling = 2*ns")

# make sure user has set ns or ps.
runtime_unit_query = lines[6].rstrip().replace(" ","").split("=")[-1].split("*")[1]
if runtime_unit_query not in ["ns", "ps"]:
    raise NameError("Input runtime unit not supported. Please use 'ns' or 'ps'" \
    +" on the seventh line of protocol.dat in the shape of (e.g.):\nsampling = 2*ns")

if runtime_unit_query == "ns":
    runtime_unit = BSS.Units.Time.nanosecond
elif runtime_unit_query == "ps":
    runtime_unit = BSS.Units.Time.picosecond

### get the number of lambda windows for this pert.
num_lambda = None
with open("lambdas_per_pert.dat", "r") as lambdas_file:
    reader = csv.reader(lambdas_file)
    for row in reader:

        if (row[0] == sys.argv[1] and row[1] == sys.argv[2]) or \
        (row[1] == sys.argv[1] and row[0] == sys.argv[2]):
            num_lambda = int(row[2])
if not num_lambda:
    raise NameError(f"The perturbation {sys.argv[1]}~{sys.argv[2]} (or the reverse) was not found in network.dat.")

# define the free energy protocol with all this information. User could customise settings further here, see docs.
freenrg_protocol = BSS.Protocol.FreeEnergy(num_lam=num_lambda, runtime=runtime_query*runtime_unit)


############# Set up the directory environment.
# testing is already done by BSS.
print(f"Setting up {engine_query} directory environment in outputs/{engine_query}/{sys.argv[1]}~{sys.argv[2]}.")
# set up a SOMD bound+free folder with standard settings.
print("Bound..")
BSS.FreeEnergy.Binding(
                    merged_system_solvated, 
                    freenrg_protocol, 
                    engine=f"{engine_query}",
                    work_dir=f"outputs/{engine_query}/{sys.argv[1]}~{sys.argv[2]}"
)

# set up a SOMD free + vacuum folder. Note that the free folder is overwritten, which is what we want.
print("Free..")
BSS.FreeEnergy.Solvation(
                    merged_system_solvated, 
                    freenrg_protocol, 
                    engine=f"{engine_query}",
                    work_dir=f"outputs/{engine_query}/{sys.argv[1]}~{sys.argv[2]}"
)













