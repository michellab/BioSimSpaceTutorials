import BioSimSpace as BSS
import os

input_dir = '/work/braine_md101_scratch/dlukauskis/BSS_tutorial/input'

if not os.path.isfile(f'{input_dir}/solvated.prm7'):
    # Parameterise the protein and the ligand
    protein = BSS.IO.readMolecules(f'{input_dir}/protein.pdb')
    protein_water = protein.getWaterMolecules()

    protein_search = protein.search('not water')
    protein = protein_search.getResult(0)

    protein = BSS.Parameters.ff14SB(protein, water_model="tip3p").getMolecule()

    # parameterise the ligand
    ligand = BSS.IO.readMolecules(f'{input_dir}/ligand.mol2').getMolecule(0)
    ligand = BSS.Parameters.gaff2(ligand,
                              net_charge=BSS.Parameters.formalCharge(ligand)).getMolecule()

    # Here the system is your protein plus the merged molecule.
    box_min, box_max = protein.getAxisAlignedBoundingBox()

    # Work out the box size from the difference in the coordinates.
    box_size = [y - x for x, y in zip(box_min, box_max)]

    # how much to pad each side of the protein (nonbonded cutoff = 10 A)
    padding = 15 * BSS.Units.Length.angstrom

    box_length = max(box_size) + 2*padding

    cmplx = protein + ligand

    solvated = BSS.Solvent.tip3p(molecule = cmplx.toSystem(),
                            box=3*[box_length])

    BSS.IO.saveMolecules(f'{input_dir}/solvated',solvated, ['PDB','RST7','PRM7'])

else:
    solvated = BSS.IO.readMolecules([f'{input_dir}/solvated.rst7',f'{input_dir}/solvated.prm7'])

if not os.path.isfile('minimised.pdb'):
    # minimize
    protocol = BSS.Protocol.Minimisation()
    process = BSS.Process.OpenMM(solvated, protocol, platform='CUDA')

    # Start the process in the background.
    process.start()

    # Wait for the process to finish.
    process.wait()

    # Get the minimised molecular system.
    minimised = process.getSystem()

    BSS.IO.saveMolecules('minimised',minimised, ['PDB','RST7','PRM7'])

else:
    minimised = BSS.IO.readMolecules(['minimised.rst7','minimised.prm7'])

if not os.path.isfile('equilibrated.pdb'):
    # equilibrate
    protocol = BSS.Protocol.Equilibration(runtime=2*BSS.Units.Time.nanosecond)
    process = BSS.Process.OpenMM(minimised, protocol, platform='CUDA')

    # Start the process in the background.
    process.start()

    # Wait for the process to finish.
    process.wait()

    # Get the equilibrated system
    equilibrated = process.getSystem()

    BSS.IO.saveMolecules('equilibrated',equilibrated, ['PDB','RST7','PRM7'])
else:
    equilibrated = BSS.IO.readMolecules(['equilibrated.rst7','equilibrated.prm7'])

p1, p2 = BSS.Metadynamics.CollectiveVariable.makeFunnel(equilibrated)

new_upper_bound = BSS.Metadynamics.Bound(value=3.5*BSS.Units.Length.nanometer)
funnel_cv = BSS.Metadynamics.CollectiveVariable.Funnel(p1, p2, upper_bound = new_upper_bound)

protocol = BSS.Protocol.Metadynamics(funnel_cv, 
                                     runtime = 100*BSS.Units.Time.nanosecond,
                                     hill_height = 1.5*BSS.Units.Energy.kj_per_mol,
                                     restart_interval = 50000,
                                     bias_factor = 10)

BSS.Process.OpenMM(equilibrated, protocol, platform='CUDA',work_dir='fun-metaD').run()
