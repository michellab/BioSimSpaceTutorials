## 1. Setting up a Free Energy Perturbation directory environment for a congeneric series of ligands

Computational chemists can support structure-activity relationship studies in medicinal chemistry by making computer models that can predict binding affinity of ligands to proteins. One of the major techniques to do this is Free Energy Perturbation (FEP). This chapter will outline the steps needed to use BioSimSpace to set up files needed for a standard FEP run.

![image-20210416115353759](./inputs/tut_imgs/tyk2_protlig.png)

​															*Figure 1: Tyrosine kinase 2 (TYK2) structure with bound ligand (ejm_48).*

#### 1.1. Importing libraries and files

First, we need to import BioSimSpace.

```python
import BioSimSpace as BSS
```

Additionally, we will need to make python navigate the folders we are working with (to make it see our input files). We will use the module [glob](https://docs.python.org/3/library/glob.html) for this purpose.

```python
import glob
```

Now that we have the required modules, we can start our work by loading our input files. We have ligands for the protein target TYK2 (a common benchmarking set in the FEP field). First, we need to tell ```glob``` in what folder the files are:

```python
path_to_ligands = "inputs/ligands"
path_to_protein = "inputs/protein/tyk2"
```

Using BioSimSpace, we can load our protein structure easily with standard syntax:

```python
protein = BSS.IO.readMolecules(["{}.rst7".format(path_to_protein), "{}.prm7".format(path_to_protein)])[0]
```

Loading our ligands is slightly more complicated as we have 16 of them. Instead of writing 16 lines of code with individual file readings, we can use ```glob``` to search through the directory we specified earlier. Together with the ligand structures we will also derive the names of the ligands to use later on in the FEP workflow.

```python
# glob returns a list of the filepaths of all PDB files in the specified directory.
ligand_files = glob.glob("{}/*.pdb".format(path_to_ligands))

# initiate empty lists to append objects to.
ligands = []
ligand_names = []

# loop over each ligand filepath,
for filepath in ligand_files:
  
    # use BSS to load the molecule object and append to the list.
    ligands.append(BSS.IO.readPDB(filepath)[0])
    
    # append the molecule name to another list so that we can use the name later on in the FEP workflow.
    ligand_names.append(filepath.split("/")[-1].replace(".pdb",""))
```

To summarise, we now have:

- a parameterised BSS protein object (TYK2)
- a list of BSS ligand objects (ligands with binding pose in TYK2)
- a list of ligand names (based on file names) of which the order corresponds to the order of ligand objects in the above list

Now that we have loaded all the input files we need, we can move on to the next step which is figuring out perturbations to simulate in our congeneric series.

#### 1.2. Using LOMAP to select molecular transformations

BioSimSpace uses software called [Lead Optimization Mapper (LOMAP)](https://github.com/MobleyLab/Lomap) to generate perturbation networks needed for FEP. LOMAP will, given our collection of ligands, automatically select which transformations to simulate. It does this based on a set of rules, where typically smaller transformations are favoured. Although typically you would run LOMAP from the commandline, we can use the BioSimSpace LOMAP functionality directly. 

```python
tranformations, lomap_scores = BSS.Align.generateNetwork(ligands)
```

LOMAP has an internal function that calculates a *LOMAP score* (between 0 and 1) that estimates whether a transformation is likely to be reliable (1), or unlikely to be reliable (0). We can check if we are happy with the edges that LOMAP made by looking at the output that BioSimSpace has generated.

```python
# Because BioSimSpace outputs transformations as a list of indices, we need to get the ligand names
# from the list of ligand names we made before.
transformations_named = [(ligand_names[transf[0]], ligand_names[transf[1]]) for transf in tranformations]

# now we can merge the named transformations and lomap scores together to see them.
for transf, score in zip(transformations_named, lomap_scores):
    print(transf, score)
```

Which will show us the *LOMAP score* per transformation:

('ejm_49', 'ejm_55') 0.2466
('ejm_49', 'ejm_31') 0.27253
('ejm_48', 'ejm_55') 0.27253
('ejm_48', 'ejm_31') 0.30119
('jmc_28', 'jmc_30') 0.90484
('jmc_28', 'jmc_23') 0.95123
('jmc_28', 'ejm_31') 0.33287
('jmc_30', 'jmc_23') 0.86071
('ejm_54', 'ejm_55') 0.86071
('ejm_54', 'ejm_42') 0.86071
('ejm_55', 'ejm_42') 0.95123
('ejm_55', 'ejm_47') 0.30119
('jmc_27', 'ejm_46') 0.90484
('jmc_27', 'jmc_23') 0.95123
('ejm_43', 'ejm_42') 0.90484
('ejm_43', 'ejm_44') 0.90484
('ejm_42', 'ejm_44') 0.81873
('ejm_42', 'ejm_50') 0.95123
('ejm_42', 'ejm_31') 0.90484
('ejm_46', 'jmc_23') 0.90484
('ejm_46', 'ejm_31') 0.36788
('ejm_47', 'ejm_31') 0.33287
('ejm_45', 'ejm_31') 0.74082
('ejm_50', 'ejm_31') 0.90484

Overall this looks good, although we have a few ligands (such as ejm_48 and ejm_49) which have poor LOMAP scores. This typically happens because they are dissimilar to the rest of the ligand series, and unfortunately unreliable estimates are difficult to prevent for these cases. One way to increase the reliability of these predictions is to increase their connections to the network; even though these new transformations are likely to also be unreliable, the added prediction will help with statistical analyses later on in the FEP workflow.

Adding a transformation is done very simply appending the desired transformation to the list of transformtions. In our case, we will add ejm_49&#8652;ejm_54.

```python
transformations_named.append(('ejm_49', 'ejm_54'))
```

Now that we have the transformations that we want to simulate, we can move on to the next step which is using BioSimSpace to set up directories to run FEP on. 

*Note that in state-of-the-art FEP software GUIs, users can view the suggested perturbation network with molecular structures. Although there is an option for BioSimSpace to output an image of this network, generating a well formatted network without overlapping molecular structures is difficult with python modules alone. We therefore recommend visualising molecular structures in external molecular viewing software, such as [PyMol](https://pymol.org/2/) or [VMD](https://www.ks.uiuc.edu/Research/vmd/).*

#### 1.3. Setting up simulation folders to run FEP on

Our ```protein``` object already contains a parameterised TYK2 structure, but our ```ligands``` were loaded from PDB files. For this tutorial, we also have pre-parameterised ligands - used force fields are FF14SB and GAFF2.

First, we load the pre-parameterised ligands into a python dictionary, which will hold both the ligand names and ligand objects:

```python
# initiate empty dictionary to save ligand names and ligand objects to.
ligands_p = {}

for ligand_name in ligand_names:
  
  # for each ligand name, we read the corresponding files using the path to the ligands we defined earlier.
  ligands_p[ligand_name] = BSS.IO.readMolecules([
                    "{}/{}.rst7".format(path_to_ligands,ligand_name),
                    "{}/{}.prm7".format(path_to_ligands,ligand_name),
                     ])[0]
```

Then, we define the FEP protocol that we would like to run. Depending on the user's goals, there are many settings that can be changed; [see the BioSimSpace documentation for the complete set of options](https://biosimspace.org/api/generated/BioSimSpace.Protocol.FreeEnergy.html#BioSimSpace.Protocol.FreeEnergy). One of the most typical settings to change in FEP is the number of &#955; windows. For demonstration purposes we will set this to seven. 

```python
protocol = BSS.Protocol.FreeEnergy(num_lam=7)
```

In FEP it is good practice to run both directions of an edge for statistical purposes. We will generate a new list of transformations by looping over the LOMAP generated transformations and appending both the forward and reverse perturbation name to a new list.

```python
transformations_named_2x = []

for pert in transformations_named:
    transformations_named_2x.append((pert[0], pert[1]))
    transformations_named_2x.append((pert[1], pert[0]))
```

Now that everything is in place, we can loop over the transformations specified by LOMAP and follow the standard BioSimSpace code for setting up a SOMD folder. The steps in this workflow are:

1. generate an atom mapping between the two ligands
2. align the two ligands using a Maximum Common Substructure (MCS) approach and the atom mapping 
3. merge the aligned ligands together into a single *merged* molecular object
4. add the *merged* object into the protein structure
5. solvate the system (i.e. the protein with *merged* object)
6. given a protocol, perturbation name and solvated system, set up a directory containing all files necessary for simulation

Which in python code would look like this:

```python
# We can loop over each transformation and do the standard BSS setup protocol.
for pert in transformations_named_2x:
    
    # get the molecule objects from our dictionary.
    lig_1 = ligands_p[pert[0]]
    lig_2 = ligands_p[pert[1]]

    # derive the perturbation name (to name our simulation folder).
    pert_name = pert[0]+"~"+pert[1]

    # 1. generate a mapping between the two molecules.
    mapping = BSS.Align.matchAtoms(lig_1, lig_2)

    # 2. align ligand A to ligand B.
    lig_1_a = BSS.Align.rmsdAlign(lig_1, lig_2, mapping)

    # 3. merge the aligned molecules into a single object.
    merged = BSS.Align.merge(lig_1_a, lig_2, mapping)

    # 4. insert the merged molecules into the protein and solvate the system.
    system = protein + merged

    # 5. solvate the system into a 10x10x10 nm water box.
    system = BSS.Solvent.tip3p(molecule=system, box=3*[10*BSS.Units.Length.nanometer])

    # 6. set up a SOMD folder with the specified protocol settings.
    BSS.FreeEnergy.Binding(
                        system, 
                        protocol, 
                        engine="SOMD",						# this can also be set to GROMACS if required.
                        work_dir="outputs/SOMD/"+pert_name
    )
```
Alternatively, free energy of solvation can also be computed using FEP. Whereas in the above code we have used bound and solvated legs to close our thermodynamic cycle, we can also use solvated and vacuum legs instead. This technique allows us to reduce computing time (but obviously does not give us binding affinity estimates). For such a workflow only the following code will have to be changed:

- instead of making a system and solvating it, the *merged* molecule object can be solvated directly
- instead of using BSS.FreeEnergy.Binding(), use BSS.FreeEnergy.Solvation()

#### Final remarks

In this chapter we have set up a FEP directory environment that is ready to be simulated using BioSimSpace. At this point, the directories setup in ```./outputs/SOMD/``` will have to be moved to a machine that contains GPUs, unless the above code has been run on such a system already. 

For this tutorial, we will assume that the above code has been run on a workstation that does not contain the hardware to run FEP simulations. Therefore, the next chapter in this tutorial will outline ways to run FEP on an external computing cluster. 