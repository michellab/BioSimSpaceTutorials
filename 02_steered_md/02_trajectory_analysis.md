Author: Adele Hardie

Email: adele.hardie@ed.ac.uk

#### Requirements:
* BioSimSpace
* numpy
* pandas
* matplotlib
* A steered MD trajectory

## Analysing sMD trajectories
The purpose of the steered MD simulation is to access conformational space that would take a very long time (or be inaccessible altogether) at equilibrium. To generate the data for the Markov State Model, we need to see how the system behaves given some starting conformation. To do this, we will be running seeded MD simulations, where a snapshot from the sMD trajectory is used as a starting point for an equilibrium MD simulation.

Start by importing required libraries:
```python
import pandas as pd
import numpy as np
import BioSimSpace as BSS
import os
```

#### Visualize the CV progression
PLUMED outputs a file with the CV value that was used for steering.
```python
steering_output_file = 'data/COLVAR'
```
This can be read into a pandas dataframe, and the plotted. Note that the time column was originally in ps and the rmsd column was originally in nm (default PLUMED output) and have been converted to ns and &#8491;.
```python
df = pd.read_csv(steering_output_file, sep=' ', comment='#', skipinitialspace=1, names=['time/ns', 'rmsd', 'bias'])
df['time/ns'] = df['time/ns']/1000
df['rmsd'] = df['rmsd']*10
df.set_index('time/ns', inplace=True)
```
```python
ax = df['rmsd'].plot(figsize=(12,5))
ax.set_ylabel('RMSD/$\AA$')
ax.set_xlim(0, 152)
```
<img src="figures/COLVAR_all.png">
We can see how the WPD loop in PTP1B conformation moves closer and closer to the closed loop crystal structure.

#### Select snapshots
In this case we will be extracting 100 snapshots and saving them as PDB files.
```python
snapshot_dir = 'snapshots'
if not os.path.exists(snapshot_dir):
    os.mkdir(snapshot_dir)
```
numpy is used to get evenly spaced indices. The end point used is the end of the steering stage, rather than the end of the entire simulation.
```python
number_of_snapshots = 100
frames = np.linspace(0, 30000, number_of_snapshots, dtype=int)
```
In this case, the RMSD was changing roughly evenly throughout the simulation. However, just to be safe, let's check the RMSD of the chosen 100 frames only.
```python
ax = df['rmsd'].iloc[frames].plot(figsize=(12,5))
ax.set_ylabel('RMSD/$\AA$')
ax.set_xlim(0, 152)
```
<img src="figures/COLVAR_snapshots.png">
The snapshots sample a large range of RMSD values.

#### Save snapshots
Each snapshot is then saved as a PDB. We can use BioSimSpace to load a specific frame of a trajectory only.
```python
for i, index in enumerate(frames):
    frame = BSS.Trajectory.getFrame(trajectory='data/sMD.nc', topology = 'data/system.prm7', index=int(index))
    BSS.IO.saveMolecules(f'{snapshot_dir}/snapshot_{i+1}', frame, 'pdb')
```
These PDB files are to be used as starting points for 100 individual 100 ns simulations, starting with resolvation, minimisation and equilibration. This is very time consuming and best done on an HPC cluster.
