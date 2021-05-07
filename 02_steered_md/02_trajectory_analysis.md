# Steered MD in BioSimSpace

Author: Adele Hardie

Email: adele.hardie@ed.ac.uk

#### Requirements:
* BioSimSpace
* numpy
* pandas
* matplotlib
* A steered MD trajectory

The purpose of the steered MD simulation is to access conformational space that would take a very long time (or be inaccessible altogether) at equilibrium. To generate the data for the Markov State Model, we need to see how the system behaves given some starting conformation. To do this, we will be running seeded MD simulations, where a snapshot from the sMD trajectory is used as a starting point for an equilibrium MD simulation.

Start by importing required libraries:


```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import BioSimSpace as BSS
import os
```

## Plot steering output

PLUMED outputs a file with the CV value that was used for steering, so we can see the progression during the simulation. The file is automatically named `COLVAR` by BioSimSpace. Here it is loaded as a pandas dataframe.


```python
steering_output_file = 'data/COLVAR'
df = pd.read_csv(steering_output_file, sep=' ', skipinitialspace=1, comment='#', header=None)
```

PLUMED outputs time in picoseconds and RMSD in nanometers. For easier plotting, we change time to nanoseconds and RMSD to Angstrom.


```python
df[0] = df[0]/1000
df[1] = df[1]*10
df.set_index(0, inplace=True)
```

Now the RMSD change can be plotted:


```python
ax = df[1].plot(figsize=(12,5))
ax.set_ylabel('RMSD/$\AA$')
ax.set_xlabel('time/ns')
ax.set_xlim(0, 152)
```
 
<img src="figures/COLVAR_all.png">
    
Here the loop RMSD went down to almost 1 A. This indicates that the loop conformation was very similar to the crystal structure of PTP1B with the loop closed (which was used as the target) and so we can proceed with extracting snapshots to use as seeds.

However, sMD might not work on the first try - the steering duration and force constant used is highly dependent on each individual system. Below is an example of a failed steering attempt:

<img src="figures/COLVAR_failed.png">

Here steering was carried out for 80 ns only, and the force constant used was 2500 kJ mol<sup>-1</sup>. The RMSD was decreasing as expected, but only reached around 2 A. This was deemed insufficient and a longer steering protocol with a larger force constant was decided upon in the end. Ultimately this will depend on the system you are working with and what the goal of the steering is.

## Extract snapshots

In this case we will be extracting 100 evenly spaced snapshots to be used as starting points for the seeded MD simulations.


```python
snapshot_dir = 'data'
if not os.path.exists(snapshot_dir):
    os.mkdir(snapshot_dir)
```

Get frame indices for snapshots. Note that the end point selected is not the end of the simulation, but the end of the steering part.


```python
number_of_snapshots = 100
end = 150 / 0.002
frames = np.linspace(0, end, number_of_snapshots, dtype=int)
```

Check that the snapshots roughly even sample the RMSD range:


```python
ax = df[1].iloc[frames].plot(figsize=(12,5))
ax.set_ylabel('RMSD/$\AA$')
ax.set_xlabel('time/ns')
ax.set_xlim(0, 150)
```

<img src="figures/COLVAR_snapshots.png">

Save each snapshot as a PDB:


```python
for i, index in enumerate(frames):
    frame = BSS.Trajectory.getFrame(trajectory='/home/adele/Documents/PTP1B/steering.nc', topology = '/home/adele/Documents/PTP1B/system.prm7', index=int(index))
    BSS.IO.saveMolecules(f'{snapshot_dir}/snapshot_{i+1}', frame, 'pdb')
```

These PDB files are to be used as starting points for 100 individual 100 ns simulations, starting with resolvation, minimisation and equilibration. This is very time consuming and best done on an HPC cluster. An [example script](02_run_seededMD.py) that can be used with an array submission is provided. Note that due to the additional phospho residue parameters required it is specific to PTP1B with a peptide substrate, but the `load_system()` function can be easily modified to work with other systems, while the rest of the script is transferable.
