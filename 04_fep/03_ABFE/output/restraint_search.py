# This script was used to generate the trajectory used for restraint selection.
import os
# Make GPU visible
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
import BioSimSpace.Sandpit.Exscientia as BSS
system = BSS.IO.readMolecules(["../input/complex/mif_mif180.prm7", "../input/complex/mif_mif180.rst7"])
lig = BSS.Align.decouple(system[0])
system.updateMolecule(0,lig)
protocol = BSS.Protocol.Production(runtime=0.5*BSS.Units.Time.nanosecond)
restraint_search = BSS.FreeEnergy.RestraintSearch(system, protocol=protocol, engine='gromacs', work_dir='restraint_search', gpu_support=True)
restraint_search.start()
restraint_search.wait()
