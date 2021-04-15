# Free Energy perturbation

Tutorial on running a standard free energy perturbation (FEP) run on a ligand series (TYK2). We start this tutorial with a set of pre-parameterised ligand files that contain geometries corresponding to their binding pose as well as a per-parameterised TYK2 ligand file. 

Using BioSimSpace we will: 

- load the collection of ligands as well as the protein
- plan which transformations we want to simulate
- setup the transformations in separate folders

After which these simulations can be performed on a series of GPUs. We have supplied the FEP results for the sake of this tutorial seeing as running this set might take several days, at least. Data analysis will be done using FreeNrgWorkflows:

- load all results
- load experimental affinity values
- statistical analyses and free energy estimations
- plotting several standard graphs





Non-BSS dependencies:

- freenrgworkflows==1.1.0

- GROMACS install for solvation



