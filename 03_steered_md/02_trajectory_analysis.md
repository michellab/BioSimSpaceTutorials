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
df = pd.read_csv(steering_output_file, sep=' ')
cols = list(df.columns[2:])
df = pd.read_csv(steering_output_file, sep=' ', comment='#', names=cols)
df.head()
```

The COLVAR file contains all of the CV values (r1,t1,d1), as well as more information on the force applied and work done.

PLUMED outputs time in picoseconds and RMSD in nanometers. For easier plotting, we change time to nanoseconds and distances to Angstrom.


```python
df['time'] = df['time']/1000
df['r1'] = df['r1']*10
df['d1'] = df['d1']*10
df.set_index('time', inplace=True)
```

Now the CV changes can be plotted:


```python
fig, ax = plt.subplots(3, figsize=(12,12))

columns = ['r1', 't1', 'd1']
ylabels = ['RMSD/$\AA$', 'Dihedral/radians', 'Distance/$\AA$']

for i in range(len(columns)):
    ax[i].plot(df.index, df[columns[i]], alpha=0.7)
    ax[i].set_ylabel(ylabels[i])
    ax[i].set_xlabel('time/ns')
    ax[i].set_xlim(0, 152)
    
fig.tight_layout()
```
 
<img src="figures/COLVAR_all.png">
    
Here the loop RMSD went down to below 2 A (around 1.8 A). This indicates that the loop conformation was very similar to the crystal structure of PTP1B with the loop closed (which was used as the target) and so we can proceed with extracting snapshots to use as seeds. Additionally, the Tyr152 $\chi$ 1 angle was kept in the "down" rotamer consistently (around 1.1 radians) and the Phe196(C $\gamma$ )-Phe280(C $\gamma$ ) distance was decreased to 4 A, which corresponds to the two residues $\pi$ -stacking.

However, sMD might not work on the first try - the steering duration and force constant used is highly dependent on each individual system. Below is an example of a failed steering attempt:

<img src="figures/COLVAR_failed.png">

Here steering was carried out for 80 ns only, and the force constant used was 2500 kJ mol<sup>-1</sup>. The RMSD was decreasing as expected, but didn't go below 2 A. This was deemed insufficient and a longer steering protocol with a larger force constant was decided upon in the end. Ultimately this will depend on the system you are working with and what the goal of the steering is.

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
end = 150 / 0.005
frames = np.linspace(0, end, number_of_snapshots, dtype=int)
```

Check that the snapshots roughly evenly sample the CVs range:


```python
fig, ax = plt.subplots(3, figsize=(12,12))

columns = ['r1', 't1', 'd1']
ylabels = ['RMSD/$\AA$', 'Dihedral/radians', 'Distance/$\AA$']

for i in range(len(columns)):
    ax[i].plot(df.index[frames], df.iloc[frames][columns[i]], alpha=0.7)
    ax[i].set_ylabel(ylabels[i])
    ax[i].set_xlabel('time/ns')
    ax[i].set_xlim(0, 150)
    
fig.tight_layout()
```

<img src="figures/COLVAR_snapshots.png">

Save each snapshot as a PDB:


```python
for i, index in enumerate(frames):
    frame = BSS.Trajectory.getFrame(trajectory='/home/adele/Documents/PTP1B/steering.nc', topology = '/home/adele/Documents/PTP1B/system.prm7', index=int(index))
    BSS.IO.saveMolecules(f'{snapshot_dir}/snapshot_{i+1}', frame, 'pdb')
```

These PDB files are to be used as starting points for 100 individual 100 ns simulations, starting with resolvation, minimisation and equilibration. This is very time consuming and best done on an HPC cluster. An [example script](02_run_seededMD.py) that can be used with an array submission is provided. Note that due to the additional phospho residue parameters required it is specific to PTP1B with a peptide substrate, but the `load_system()` function can be easily modified to work with other systems, while the rest of the script is transferable.
