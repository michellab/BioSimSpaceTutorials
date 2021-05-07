import BioSimSpace as BSS
from BioSimSpace import _Exceptions
import os,sys
import csv
print(BSS.__version__)

print ("%s %s %s" % (sys.argv[0],sys.argv[1],sys.argv[2]))

# Load equilibrated free inputs for both ligands. Complain if input not found. These systems already contain equil. waters.
print(f"Loading ligands {sys.argv[1]} and {sys.argv[2]}.")
ligs_path = "prep/ligands/"
ligand_1_sys = BSS.IO.readMolecules([f"{ligs_path}{sys.argv[1]}_lig_equil_solv.rst7", f"{ligs_path}{sys.argv[1]}_lig_equil_solv.prm7"])
ligand_2_sys = BSS.IO.readMolecules([f"{ligs_path}{sys.argv[2]}_lig_equil_solv.rst7", f"{ligs_path}{sys.argv[2]}_lig_equil_solv.prm7"])

# Extract ligands.
ligand_1 = ligand_1_sys.getMolecule(0)
ligand_2 = ligand_2_sys.getMolecule(0)

# Align ligand2 on ligand1
print("Mapping and aligning..")
print(ligand_1, ligand_2)
mapping = BSS.Align.matchAtoms(ligand_1, ligand_2, sanitize=True, complete_rings_only=True)
ligand_1_a = BSS.Align.rmsdAlign(ligand_1, ligand_2, mapping)

# Generate merged molecule.
print("Merging..")
merged_ligs = BSS.Align.merge(ligand_1_a, ligand_2, mapping)

#### Get equilibrated waters and waterbox information for both bound and free. Get all information from lambda==0
# Following is work around because setBox() doesn't validate correctly boxes with lengths and angles

ligand_1_sys.removeMolecules(ligand_1)
ligand_1_sys.addMolecules(merged_ligs)
system_free = ligand_1_sys



################ now repeat above steps, but for the protein + ligand systems.
# Load equilibrated bound inputs for both ligands. Complain if input not found
print(f"Loading bound ligands {sys.argv[1]} and {sys.argv[2]}.")
ligs_path = "prep/protein/"
system_1 = BSS.IO.readMolecules([f"{ligs_path}{sys.argv[1]}_sys_equil_solv.rst7", f"{ligs_path}{sys.argv[1]}_sys_equil_solv.prm7"])
system_2 = BSS.IO.readMolecules([f"{ligs_path}{sys.argv[2]}_sys_equil_solv.rst7", f"{ligs_path}{sys.argv[2]}_sys_equil_solv.prm7"])

# Extract ligands and protein. Do this based on nAtoms and nResidues, as sometimes
# the order of molecules is switched, so we can't use index alone.
# JM - bugfix in BSS making the below redundant ?
system_ligand_1 = None
protein = None
n_residues = [mol.nResidues() for mol in system_1]
n_atoms = [mol.nAtoms() for mol in system_1]
for i, (n_resi, n_at) in enumerate(zip(n_residues[:20], n_atoms[:20])):
    if n_resi == 1 and n_at > 5:
        system_ligand_1 = system_1.getMolecule(i)
    elif n_resi > 1:
        protein = system_1.getMolecule(i)
    else:
        pass

# loop over molecules in system to extract the ligand  
system_ligand_2 = None

n_residues = [mol.nResidues() for mol in system_2]
n_atoms = [mol.nAtoms() for mol in system_2]
for i, (n_resi, n_at) in enumerate(zip(n_residues, n_atoms)):
    # grab the system's ligand and the protein. ignore the waters.
    if n_resi == 1 and n_at > 5:
        system_ligand_2 = system_2.getMolecule(i)
    else:
        pass

# extract ions.
#ions_bound = system_2.search("not mols with atomidx 2")

if system_ligand_1 and system_ligand_2 and protein:
    print("Using molecules ligand_1, ligand_2, protein:")
    print(system_ligand_1, system_ligand_2, protein)
else:
    raise _Exceptions.AlignmentError("Could not extract ligands or protein from input systems. Check that your ligands/proteins are properly prepared by BSSligprep.sh!")

# Align ligand2 on ligand1
print("Mapping..")
mapping = BSS.Align.matchAtoms(system_ligand_1, system_ligand_2, sanitize=True, complete_rings_only=True)

print("Aligning..")
system_ligand_1_a = BSS.Align.rmsdAlign(system_ligand_1, system_ligand_2, mapping)

# Generate merged molecule.
print("Merging..")
system_merged_ligs = BSS.Align.merge(system_ligand_1_a, system_ligand_2, mapping)



#### Get equilibrated waters and waterbox information for both bound and free. Get all information from lambda==0
#waters_bound = system_1.getWaterMolecules()
#waterbox_bound = system_1.getBox()
# now make final systems with merged, the equil. protein of lambda==0 and equil. waters of lambda==0.
#system_bound = system_merged_ligs + protein + ions_bound + waters_bound
# restore box information.
#system_bound.setBox(waterbox_bound)

system_1.removeMolecules(system_ligand_1)
system_1.addMolecules(system_merged_ligs)
system_bound = system_1

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
with open("network.dat", "r") as lambdas_file:
    reader = csv.reader(lambdas_file, delimiter=" ")
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
# set up a bound+free folder with standard settings.
# Use solvation to prep the bound leg
# Is it necessary to duplicate pert/prm7/rst7 inputs for each lambda folder ? 


print("Bound..")
workdir=f"outputs/{engine_query}/{sys.argv[1]}~{sys.argv[2]}"
BSS.FreeEnergy.Solvation(
                    system_bound, 
                    freenrg_protocol, 
                    engine=f"{engine_query}",
                    work_dir=workdir
)

cmd = "mv %s/free %s/bound" % (workdir,workdir)
os.system(cmd)

# set up a free + vacuum folder. Note that the free folder is overwritten, which is what we want because
# we've equilibrated the ligand and ligand+protein separately.
print("Free..")
BSS.FreeEnergy.Solvation(
                    system_free, 
                    freenrg_protocol, 
                    engine=f"{engine_query}",
                    work_dir=workdir
)

cmd = "rm -rf %s/vacuum" % (workdir)
os.system(cmd)


# Until we figure out what is going on overwrite BSS generated cfg files with parameter set that works 
if engine_query == 'SOMD':
    cfg_template = """platform = CUDA
nmoves = 50000
ncycles = 10
buffered coordinates frequency = 1250
save coordinates = True
timestep = 2 * femtosecond
constraint = hbonds-notperturbed
hydrogen mass repartitioning factor = 1.0
cutoff type = cutoffperiodic
cutoff distance = 10*angstrom
barostat = True
thermostat = True
energy frequency = 250
precision = mixed
minimise = True
equilibrate = False
center solute = True
reaction field dielectric = 82.0
minimal coordinate saving = True
"""
    lams = freenrg_protocol.getLambdaValues()

    lamarraystring = "lambda array = %.4f" % lams[0]
    for lam in lams[1:]:
        lamarraystring += ", %.4f" % lam
    #print (lamarraystring)

    cfg_template += lamarraystring + "\n"
    import glob
    config_files = glob.glob(f"{workdir}/*/lambda*/somd.cfg")
    for cfg in config_files:
        f = open(cfg, 'w'); f.write(cfg_template); f.close()





