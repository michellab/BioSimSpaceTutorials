import BioSimSpace as BSS
import sys



print ("%s %s %s" % (sys.argv[0],sys.argv[1],sys.argv[2]))

# Load equilibrated free inputs for both ligands. Complain if input not found
ligs_path = "inputs/ligands/"
ligand_1 = BSS.IO.readMolecules([f"{ligs_path}{sys.argv[1]}.rst7", f"{ligs_path}{sys.argv[1]}.prm7"])
ligand_2 = BSS.IO.readMolecules([f"{ligs_path}{sys.argv[2]}.rst7", f"{ligs_path}{sys.argv[2]}.prm7"])

# Extract ligand2 from frame.
ligand_2 = ligand_2.getMolecule(0)

# Align on ligand1
mapping = BSS.Align.matchAtoms(ligand_1, ligand_2, sanitize=True)

# Generate merge molecule.
# Write output files for perturbation using engine specified in 'protocol.dat' (so prm7/rst7 and pert if SOMD)

# Repeat same procedure for equilibrated bound input








# # get the molecule objects from our dictionary.
# lig_1 = ligands_p[pert[0]]     # 
# lig_2 = ligands_p[pert[1]]


# # derive the perturbation name (to name our simulation folder).
# pert_name = pert[0]+"~"+pert[1]

# # generate a mapping between the two molecules.
# mapping = BSS.Align.matchAtoms(lig_1, lig_2)

# # align ligand A to ligand B.
# lig_1_a = BSS.Align.rmsdAlign(lig_1, lig_2, mapping)

# # merge the aligned molecules into a single object.
# merged = BSS.Align.merge(lig_1_a, lig_2, mapping)

# # insert the merged molecules into the protein and solvate the system.
# system = protein + merged

# # no GROMACS install on current system; uncomment in tutorial.
# system = BSS.Solvent.tip3p(molecule=system, box=3*[10*BSS.Units.Length.nanometer])

# # set up a SOMD folder with standard settings.
# BSS.FreeEnergy.Binding(
#                     system, 
#                     protocol, 
#                     engine="SOMD",
#                     work_dir="outputs/SOMD/"+pert_name
# )