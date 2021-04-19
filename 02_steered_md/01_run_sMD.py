import BioSimSpace as BSS
import os
from argparse import ArgumentParser
from shutil import copyfile

def create_CV(reference, system):
    """Create a collective variable to use in sMD

    Parameters
    ----------
    reference : BioSimSpace._SireWrappers._molecule.Molecule
        reference structure for RMSD calculation

    system : BioSimSpace._SireWrappers._system.System
        steered system

    Returns
    -------
    rmsd_cv : BioSimSpace.Metadynamics.CollectiveVariable._rmsd.RMSD
        CV for use in steering
    """
    rmsd_indices = []
    for residue in reference.getResidues():
        if 178<=residue.index()<=184:
            for atom in residue.getAtoms():
                if atom.element()!='Hydrogen (H, 1)':
                    rmsd_indices.append(atom.index())

    rmsd_cv = BSS.Metadynamics.CollectiveVariable.RMSD(system, reference, 0, rmsd_indices)

    return rmsd_cv

def setup_sMD_protocol(steering_runtime, total_runtime, force_constant, cv):
    """Create steered MD protocol

    Parameters
    ----------
    steering_runtime : BioSimSpace.Units.Time
        steering time. First 4 ps will be used to apply force

    total_runtime : BioSimSpace.Units.Time
        total simulation time. Usually steering_runtime + some short duration for relaxation of system

    force_constant : int
        force constant to use for steering

    Returns
    -------
    protocol : BioSimSpace.Protocol._steering.Steering
        steering protocol
    """
    start = 0* BSS.Units.Time.nanosecond
    apply_force = 4 * BSS.Units.Time.picosecond

    nm = BSS.Units.Length.nanometer
    restraint_1 = BSS.Metadynamics.Restraint(cv.getInitialValue(), 0)
    restraint_2 = BSS.Metadynamics.Restraint(cv.getInitialValue(), force_constant)
    restraint_3 = BSS.Metadynamics.Restraint(0*nm, force_constant)
    restraint_4 = BSS.Metadynamics.Restraint(0*nm, 0)

    protocol = BSS.Protocol.Steering(cv, [start, apply_force, steering_runtime, total_runtime], [restraint_1, restraint_2, restraint_3, restraint_4], runtime=total_runtime)

    return protocol

def setup_in_AMBER(system, protocol, cv):
    """setup an sMD process to run with AMBER

    Parameters
    ----------
    system : BioSimSpace._SireWrappers._system.System
        steered system

    protocol : BioSimSpace.Protocol._steering.Steering
        steering protocol

    cv : BioSimSpace.Metadynamics.CollectiveVariable._rmsd.RMSD
        CV for use in steering

    Returns
    -------
    process : BioSimSpace.Process._amber.Amber
        steering process
    """
    #get plumed input
    plumed = BSS.Process._plumed.Plumed('.')
    conf, files = plumed.createConfig(system, protocol)

    #create new protocol and process in Amber
    protocol = BSS.Protocol.Production(runtime=protocol.getRunTime())
    process = BSS.Process.Amber(system, protocol, exe=f'{os.environ["AMBERHOME"]}/bin/pmemd.cuda')

    #write input files
    plumed_file = open(f'{process.workDir()}/plumed.in', 'w')
    plumed_file.writelines(map(lambda x: x + '\n', conf))
    plumed_file.close()
    ref_file = open(f'{process.workDir()}/reference.pdb', 'w')
    ref_file.writelines(map(lambda x: x + '\n', cv.getReferencePDB()))
    ref_file.close()

    #adjust config
    process.setConfig(process.getConfig()[:-1]+['  plumed=1,', '  plumedfile="plumed.in",', ' /'])

    return process

def run_sMD(topology, coordinates, reference, steering_runtime, total_runtime, force_constant):
    """Run steered MD with PLUMED and BioSimSpace.

    Parameters
    ----------
    topology : str
        AMBER topology file

    coordinates : str
        equilibrated system coordinates

    reference : str
        steering target PDB structure

    steering_runtime : BioSimSpace.Units.Time
        steering time

    total_runtime : BioSimSpace.Units.Time
        total simulation runtime

    force_constant : int
        force to use for steering

    Returns
    -------
    None
    """
    system = BSS.IO.readMolecules([topology, coordinates])
    reference = BSS.IO.readMolecules(reference).getMolecule(0)

    rmsd_cv = create_CV(reference, system)

    protocol = setup_sMD_protocol(steering_runtime, total_runtime, force_constant, rmsd_cv)

    process = setup_in_AMBER(system, protocol, rmsd_cv)
    copyfile(topology, f'{process.workDir()}/amber.prm7')

    process.start()
    process.wait()

    files = ['COLVAR', 'stdout', 'amber.out', 'amber.nc']
    for file in files:
        copyfile(f'{process.workDir()}/{file}', file)

    return None

def __main__():
    parser = ArgumentParser(description='Set up and run sMD with PLUMED and BioSimSpace')
    parser.add_argument('--topology', type=str, required=True, help='AMBER topology file')
    parser.add_argument('--coordinates', type=str, required=True, help='equilibrated system coordinates')
    parser.add_argument('--reference', type=str, required=True, help='target molecule PDB')
    parser.add_argument('--steering_runtime', type=float, required=True, help='steering time in ns')
    parser.add_argument('--total_runtime', type=float, required=True, help='total simulation runtime')
    parser.add_argument('--force', type=int, required=True, help='steering force_constant')
    args = parser.parse_args()

    run_sMD(args.topology, args.coordinates, args.reference, args.steering_runtime*BSS.Units.Time.nanosecond, args.total_runtime*BSS.Units.Time.nanosecond, args.force)

    return None

if __name__=='__main__':
    __main__()
