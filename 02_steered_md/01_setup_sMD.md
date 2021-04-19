Author: Adele Hardie

Email: adele.hardie@ed.ac.uk



## Introduction

Allosteric inhibition can be a useful alternative to conventional protein target-ing when the nature of the active site makes it difficult to design binders. This requires  assessment  of  whether  an  allosteric  binder  actually  has  an  effect  onprotein function, such as whether its presence shifts protein conformational ensemble to favour the inactive state. This can be modelled as a Markov chain by Markov State Models(MSMs).  Since the system is treated as memoryless, model building only requires local equilibrium. Therefore, it can make use of shorter MD simulations, allowing them to be run in parallel.

In order to have a more complete view of the protein ensemble, enhanced sampling methods are used, among them steered MD (sMD) (1). It introduces a bias potential that is added to the Hamiltonian, thus biasing the simulation towards a specified value of a chosen collective variable. Once the system has reached a certain conformation, those coordinates (2) can be used as starting points for equilibrium MD simulations (4) that can subsequently be used as data for constructing an MSM (4). An example summary of this is shown below:
<img src="figures/ensemble-md-protocol.png" width=300>

PLUMED is a library that, among other things, has enhanced sampling algorithms. It works with multiple MD engines, including GROMACS and AMBER. PLUMED uses a [moving restraint](https://www.plumed.org/doc-v2.5/user-doc/html/_m_o_v_i_n_g_r_e_s_t_r_a_i_n_t.html) that is calculated as follows:

$V(\vec{s},t) = \frac{1}{2} \kappa(t) ( \vec{s} - \vec{s}_0(t) )^2$     (Eq. 1)

where $\vec{s}_0$ and $\kappa$ are time dependent and specified in the PLUMED input. $\vec{s}_0$ is the target CV value and $\kappa$ is the force constant in kJ mol$^{-1}$. The values of both of them are set at specific steps, and linearly interpolated in between.

This tutorial focuses on running the prerequisite simulations using BSS.The example system used is protein tyrosine phosphatase 1B(PTP1B), which exists in two dominant conformations: WPD loop open and WPD loop closed:
<img src="figures/open-close.png" width=250>

## Set up steered MD

Running steered MD in BioSimSpace is very similar to regular simulations already covered. It only requires some additional preparation for interfacing with PLUMED, the software that takes care of biasing the Hamiltonian.

#### Setting up the system

We start by importing the required libraries.

```python
import BioSimSpace as BSS
```

Load a system with BioSimSpace. This particular system is of PTP1B with the WPD loop open (from PDB entry 2HNP) with a peptide substrate and has been minimised and equilibrated.

```python
system = BSS.IO.readMolecules(['data/system.prm7', 'data/system.rst7'])
```

#### Creating the CV

A collective variable is required to run sMD, as this is the value that will be used to calculate the biasing potential. In this case, the CV is RMSD of the heavy atoms in the WPD loop (residues 178-184) when the WPD loop is closed (i.e. steering the loop from open to closed conformation). Let's load this reference structure.

```python
reference = BSS.IO.readMolecules('data/reference.pdb').getMolecule(0)
```

Since not all of the atoms in the reference will be used to calculate the RMSD, we check all the residues and append the appropriate atom indices to the `rmsd_indices` list. Here we check all the residues instead of directly accessing the residue list in case there are some residues missing in the structure.

```python
rmsd_indices = []
for residue in reference.getResidues():
    if 178<=residue.index()<=184:
        for atom in residue.getAtoms():
            if atom.element()!='Hydrogen (H, 1)':
                rmsd_indices.append(atom.index())
```

Once we have our system and reference, and we know which atoms will be used for the RMSD calculation, we can create a `CollectiveVariable` object.

```python
rmsd_cv = BSS.Metadynamics.CollectiveVariable.RMSD(system, reference, 0, rmsd_indices)
```

One thing to note when dealing with RMSD between two different structures, is that the atoms may not be in the same order. For example, atom 1 in `system` in this case is a hydrogen, whereas in `reference` it is an oxygen. BioSimSpace takes care of this by matching up the atoms in the system to the atoms in the reference.

#### Setting up a steered MD protocol

