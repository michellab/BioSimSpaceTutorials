import BioSimSpace as BSS
import os
from shutil import copyfile
from argparse import ArgumentParser

def load_system(snapshot, phosphate_parameters, output_dir=''):
    """Load a snapshot, resolvate it and parameterise.
    
    Parameters
    ----------
    snapshot : str
        snapshot PDB path
    
    phosphate_parameters : str
        path to additional phospho residue parameters
        
    output_dir : str
        output directory
    
    Returns
    -------
    system : BioSimSpace._SireWrappers._system.System
        parameterised and solvated system
    """
    system = BSS.IO.readMolecules(snapshot)
    
    #pull out protein and substrate
    for i in range(system.nMolecules()):
        if system.getMolecule(i).nResidues()==11:
            peptide = system.getMolecule(i)
        elif system.getMolecule(i).nResidues()>200:
            protein = system.getMolecule(i)
        
    #parameterise
    protein_parm = BSS.Parameters.ff14SB(protein).getMolecule()
    peptide_parm = BSS.Parameters.ff14SB(peptide, leap_commands=[f'addPath {phosphate_parameters}', 'source leaprc.phosaa10']).getMolecule()
    new_system = protein_parm.toSystem()
    new_system.addMolecules(peptide_parm)
    BSS.IO.saveMolecules(f'{output_dir}system_dry', new_system, 'prm7')
    
    #solvate
    solvated = BSS.Solvent.tip3p(new_system, shell=10*BSS.Units.Length.angstrom, ion_conc=0.15, is_neutral=True)
    BSS.IO.saveMolecules(f'{output_dir}system', solvated, ['prm7', 'rst7'])
    solvated = BSS.IO.readMolecules([f'{output_dir}system.prm7', f'{output_dir}system.rst7'])
    
    return solvated

def minimise(system, steps=2000, output_dir=''):
    """Minimise system with AMBER.
    
    Parameters
    ----------
    system : BioSimSpace._SireWrappers._system.System
        system to minimise
    
    steps : int
        minimisation steps
    
    output_dir : str
        output directory
    
    Returns
    -------
    minimised : BioSimSpace._SireWrappers._system.System
        minimised system
    """
    protocol = BSS.Protocol.Minimisation(steps)
    process = BSS.Process.Amber(system, protocol)
    minimised = process.start().getSystem(block=True)

    BSS.IO.saveMolecules(f'{output_dir}system_minimised', minimised, 'rst7')
    
    return minimised

def heat(system, runtime=20*BSS.Units.Time.picosecond, output_dir=''):
    """Heat system to 300 K with AMBER in the NVT ensemble.
    
    Parameters
    ----------
    system : BioSimSpace._SireWrappers._system.System
        system to heat
    
    runtime : BioSimSpace.Units.Time
        heating time
    
    output_dir : str
        output directory
    
    Returns
    -------
    heated : BioSimSpace._SireWrappers._system.System
        heated system
    """
    protocol = BSS.Protocol.Equilibration(runtime=runtime, temperature_start=0*BSS.Units.Temperature.kelvin, temperature_end=300*BSS.Units.Temperature.kelvin)
    process = BSS.Process.Amber(system, protocol)
    heated = process.start().getSystem(block=True)
    
    BSS.IO.saveMolecules(f'{output_dir}system_heated', heated, 'rst7')
    
    return heated

def equilibrate(system, runtime=100*BSS.Units.Time.picosecond, output_dir=''):
    """Equilibrate a system with AMBER in the NPT ensemble. Requires $AMBERHOME to be set.
    
    Parameters
    ----------
    system : BioSimSpace._SireWrappers._system.System
        system to equilibrate
    
    runtime : BioSimSpace.Units.Time
        equilibration time. The first 10 ps will be run with sander and the rest with pmemd.cuda
    
    output_dir : str
        output directory
    
    Returns
    -------
    equilibrated : BioSimSpace._SireWrappers._system.System
        equilibrated system
    """
    protocol = BSS.Protocol.Equilibration(runtime=runtime, temperature=300*BSS.Units.Temperature.kelvin, pressure=1*BSS.Units.Pressure.atm)
    process = BSS.Process.Amber(system, protocol)        
    equilibrated = process.start().getSystem(block=True)
    
    BSS.IO.saveMolecules(f'{output_dir}system_equilibrated.rst7', equilibrated, 'rst7')
    
    return equilibrated

def run_MD(system, runtime=100*BSS.Units.Time.nanosecond, output_dir=''):
    """Run equilibrium MD with AMBER (pmemd.cuda). Requires $AMBERHOME to be set.
    
    Parameters
    ----------
    system : BioSimSpace._SireWrappers._system.System
        prepared system
    
    runtime : BioSimSpace.Units.Time
        simulation duration
    
    output_dir : str
        output directory
    
    Returns
    -------
    None
    """
    protocol = BSS.Protocol.Production(runtime=runtime)
    process = BSS.Process.Amber(system, protocol, exe=f'{os.environ["AMBERHOME"]}/bin/pmemd.cuda')
    
    process.start()
    process.wait()
    
    #copy over results
    [copyfile(f'{process.workDir()}/{file}', f'{output_dir}{file}') for file in ['amber.nc', 'stdout']]
    
    return None

def seededMD(snapshot, phosphate_parameters, output_dir=''):
    """Set up a system from a snapshot PDB and run an equilibrium MD simulation using AMBER (pmemd.cuda). Requires $AMBERHOME to be set to a correct AMBER installation.
    
    Parameters
    ----------
    snapshot : str
        snapshot PDB path
    
    phosphate_parameters : str
        path to additional phospho residue parameters
        
    output_dir : str
        output directory
    
    Returns
    -------
    None
    """
    system = load_system(snapshot, phosphate_parameters, output_dir)
    minimised = minimise(system, output_dir=output_dir)
    heated = heat(minimised, output_dir=output_dir)
    equilibrated = equilibrate(heated, output_dir=output_dir)
    run_MD(equilibrated, output_dir=output_dir)
    
    return None

def __main__():
    parser = ArgumentParser(description='Set up and run seeded MD from snapshots. For PTP1B with peptide substrate')
    parser.add_argument('--snapshot', type=str, required=True, help='path to snapshot PDB')
    parser.add_argument('--phosphate_params', type=str, required=True, help='path to folder that contains additional phospho residue parameters')
    parser.add_argument('--output', type=str, default='', help='output directory. If not specified, output will be saved in current directory')
    args = parser.parse_args()
    
    seededMD(args.snapshot, args.phosphate_params, args.output)
    
    return None

if __name__=='__main__':
    __main__()
