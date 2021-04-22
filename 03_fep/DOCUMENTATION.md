# Free Energy Perturbation using BioSimSpace

Computational chemists can support structure-activity relationship studies in medicinal chemistry by making computer models that can predict binding affinity of ligands to proteins. One of the most popular techniques for this is Free Energy Perturbation (FEP), which relies on simulation alchemical transformations of ligands in a congeneric series, simulating them both in a protein target and in just a waterbox. Relative free energies of binding (ΔΔG in kcal/mol) can then be computed by simply subtracting the ΔΔG (in protein) and the ΔΔG (in water). Some introductory reading:

[Best Practices for Alchemical Free Energy Calculations](https://www.livecomsjournal.org/article/18378-best-practices-for-alchemical-free-energy-calculations-article-v1-0) 

[Relative Binding Free Energy Calculations in Drug Discovery: Recent Advances and Practical Considerations](https://pubs.acs.org/doi/10.1021/acs.jcim.7b00564)

[Assessment of Binding Affinity via Alchemical Free-Energy Calculations](https://pubs.acs.org/doi/10.1021/acs.jcim.0c00165)

This documentation will outline the steps needed to:

- Select a series of transformations to simulate using LOMAP
- Use BioSimSpace to set up files needed for a standard FEP run in both SOMD and GROMACS
- Run FEP using BioSimSpace on a computing cluster
- Analyse FEP simulation results
- Compile all FEP results locally and perform data analyses

For this tutorial we will be using TYK2, a common benchmarking set in the FEP field, first used by Schrödinger in their [2015 FEP+ paper](https://pubs.acs.org/doi/abs/10.1021/ja512751q).

![image-20210416115353759](./inputs/tut_imgs/tyk2_protlig.png)

​															*Figure 1: Tyrosine kinase 2 (TYK2) structure with bound ligand (ejm_48).*



![suggested_workflow-white](./inputs/tut_imgs/fep_pipeline.png)

*Figure 2: Schematic of the FEP pipeline in this report. Whereas green boxes represent notebooks run on a local machine, purple boxes represent python scripts run sequentially on a computing cluster.*

# 1. Generating a perturbation network to run FEP on a congeneric series of ligands

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
```

Loading our ligands is slightly more complicated as we have 16 of them. Instead of writing 16 lines of code with individual file readings, we can use ```glob``` to search through the directory we specified earlier. Together with the ligand structures we will also derive the names of the ligands to use later on in the FEP workflow.

```python
# glob returns a list of the filepaths of all MOL2 files in the specified directory.
ligand_files = glob.glob("{}/*.mol2".format(path_to_ligands))

ligands = []
ligand_names = []

for filepath in ligand_files:
    # append the molecule object to a list.
    ligands.append(BSS.IO.readMolecules(filepath)[0])
    
    # append the molecule name to another list so that we can use the name of each molecule in our workflow.
    ligand_names.append(filepath.split("/")[-1].replace(".mol2",""))

tranformations, lomap_scores = BSS.Align.generateNetwork(ligands)
```

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

*Note that in state-of-the-art FEP software GUIs, users can view the suggested perturbation network with molecular structures. Although there is an option for BioSimSpace to output an image of this network, generating a well formatted network without overlapping molecular structures is difficult with python modules alone. We therefore recommend visualising molecular structures in external molecular viewing software, such as [PyMol](https://pymol.org/2/) or [VMD](https://www.ks.uiuc.edu/Research/vmd/).*

#### 1.3. Setting up files for our execution model

Because typically FEP simulations are run on a cluster of GPUs, this documentation will follow that same route (see figure 2). We need to generate several files that can be moved to the cluster that will instruct BioSimSpace to set up simulations in a certain way. We could use BioSimSpace to set up the simulations locally, but doing this cluster-side saves time because there are far fewer files to upload to the server.

First, write a file that simply contains a list of all our ligand names.

```python
# write ligands file.
with open("./execution_model/ligands.dat", "w") as ligands_file:
    writer = csv.writer(ligands_file)
    for lig in ligand_names:
        writer.writerow([lig])
```

Then, write an intermediate file where all the transformations (generated with LOMAP and our manual insertion) are allocated a standard seven lambda windows. 

```python
# write lambda windows per perturbation. Separating this cell and the next one allows the user to
# edit the windows file if >7 windows are required for some perturbations.
with open("./execution_model/lambdas_per_pert.dat", "w") as lambda_file:
    writer = csv.writer(lambda_file)
    writer.writerow(["lig1", "lig2", "n_lambda"])
    for pert in transformations_named:
        
        writer.writerow([pert[0], pert[1], 7])
```

At this point, the user can use a text editor to edit the file ```./execution_model/lambdas_per_pert.dat``` if certain perturbations are thought to require more (or less) sampling. This can be done by simply changing the integer (in our case 7) to another number (at least 3). Executing the following cell will read the intermediate file again, and write a second file ```./execution_model/network.dat``` which is read by BioSimSpace on the cluster. This file contains the lambda arrays formatted such that they are processed easily by scripts on the cluster.

```python
# write perts file. Base the lambda schedule on the file generated in the previous cell.
np.set_printoptions(formatter={'float': '{: 0.4f}'.format})

with open("./execution_model/network.dat", "w") as network_file, \
open("./execution_model/lambdas_per_pert.dat", "r") as lambda_file:
    writer = csv.writer(network_file, delimiter=" ")
    reader = csv.reader(lambda_file)
    
    # skip header.
    next(reader)
    
    for pert in reader:
        # given the number of allocated lambda windows, generate an array for parsing downstream.
        lam_array = str(np.around(np.linspace(0, 1, int(pert[2])), decimals=4))
        # make the array into a format readable by bash.
        lam_array = lam_array.replace("[ ", "").replace("]", "").replace("  ", ",")
        
        # write out both directions for this perturbation.
        writer.writerow([pert[0], pert[1], lam_array])
        writer.writerow([pert[1], pert[0], lam_array])
```

Finally, we want to communicate to BioSimSpace under what conditions we want to run our simulations. We can set options in the file ```./execution_model/protocol.dat```.

```python
# create protocol. 
protocol = [
    "ligand forcefield = GAFF2/AM1-BCC",
    "protein forcefield = ff14SB",
    "solvent = TIP3P",
    "box edges = 10*angstrom",
    "box type = orthorhombic",
    "protocol = default ",
    "sampling = 2*ns",
    "engine = SOMD"
]
# write protocol to file.
with open("./execution_model/protocol.dat", "w") as protocol_file:
    writer = csv.writer(protocol_file)

    for prot_line in protocol:
        
        writer.writerow([prot_line])
```

# 2. Running FEP on a computing cluster using BioSimSpace

this will mostly consist of lines that show how to run a certain bash/py script. The main lines are:

copy ./execution_model/ to a cluster using SCP

run processFEP-slurm.sh

run analysis scripts? or does that happen from above scripts?

copy everything back to local dir



# 3. Analysing FEP results 

############### HERE SOMEWHERE PUT CHAPTER 0 OF ANALYSIS NOTEBOOK -> CHANGE TO CHPT 1 in ANALYSIS NOTEBOOK

#### 3.1. Importing libraries and data files

Now that we have used our GPU hardware to run our FEP simulations, we can start analysing our results. For this we will use [FreeEnergyNetworkAnalysis](https://github.com/michellab/freenrgworkflows) functionalities in python.

```python
import networkanalysis.networkanalysis as n_graph
import networkanalysis.plotting as n_plot
import networkanalysis.experiments as n_ex
import networkanalysis.stats as n_stats
import networkanalysis
```

Additionally, we will use [matplotlib](https://matplotlib.org/), the main plotting library for python.

```python
import matplotlib.pyplot as plt
```

Finally, we will use [pandas](https://pandas.pydata.org/), a very helpful data analysis library for python that will help us work with our data more intuitively (instead of lists of numbers, using pandas we can work with a more table-like format).

```python
import pandas as pd
```

Now that we have the required modules, we can start our work by setting paths to our input files. On our computing cluster, BioSimSpace has appended FEP results to a .csv file that we can load with ```FreeEnergyNetworkAnalysis```. Additionally, we would like to compare our FEP predictions to experimental measures to see how reliable our predictions are. 

```python
# a CSV file output by SOMD that contains perturbation name,free energy,confidence.
results_filepath = './outputs/summary.csv'

# experimental values (e.g. ic50/ki) for all ligands in our set.
exp_filepath = './inputs/exp_data_tyk2.dat'
```

*Note that in advanced stages of lead optimisation FEP runs will consist mostly of 'new' ligands, for which no experimental data is abailable. In that situation we can omit supplying FreeEnergyNetworkAnalysis with experimental data, and all we will end up with is the free energy prediction per ligand.*

#### 3.2. Loading files into FreeEnergyNetworkAnalysis

We will use ```FreeEnergyNetworkAnalysis``` to analyse our predictions. Instead of just computing ΔΔG values for each transformation, we would like to estimate the ΔΔG value for each individual ligand. There are some rather involved algorithms needed for these steps which is what ```FreeEnergyNetworkAnalysis``` takes care of. 

```python
# Creating the perturbation network to do calculations on.
pG = n_graph.PerturbationGraph()

# Populate the network with our simulation results.
pG.populate_pert_graph(results_filepath)
```

*In case you have run replicates (which is recommended), you can add more results to the graph by using the function ```pG.add_data_to_graph('/path/to/additional/runs.csv')```*.

#### 3.3. Estimating Free Energies

Given our graph, we can now compute ΔΔG values predicted by FEP for each individual ligand. 

```python
####### REPLACE WITH WLS FUNCTION WHEN READY
pG.compute_weighted_avg_paths(target_compound)

# take the FEP predictions from the graph.
computed_relative_DDGs = pG.freeEnergyInKcal
```

Next, we want to have our experimental data in the same unit as our FEP predictions (kcal/mol), in the case of our TYK2 set the experimental data is in Ki. Again, ```FreeEnergyNetworkAnalysis``` can do this for us with only a few lines of code.

```python
experiments = n_ex.ExperimentalData()

####### REPLACE WITH WLS FUNCTION WHEN READY --> NO REF
experiments.compute_DDG_from_IC50s(exp_filepath, reference=target_compound)

experimental_DDGs = experiments.freeEnergiesInKcal
```



Now that we have both our binding affinities estimated by FEP and corresponding experimental values, we can work towards plotting our results. First, we will transform our data into a so-called ```pandas DataFrame```, which is an intuitive way of viewing data. 

```python
freenrg_dict = {}

# construct dict with experimental freenrg and error.
for item in experimental_DDGs:
    ligand = list(item.keys())[0]
    freenrg = list(item.values())[0]
    error = list(item.values())[1]
    freenrg_dict[ligand] = [freenrg, error]

# append computed freenrg and error.
for item in computed_relative_DDGs:
    ligand = list(item.keys())[0]
    freenrg = list(item.values())[0]
    error = list(item.values())[1]
    freenrg_dict[ligand].append(freenrg)
    freenrg_dict[ligand].append(error)

# from the newly structured dictionary, generate a pandas dataframe.
freenrg_df = pd.DataFrame(freenrg_dict, index=["freenrg_exp", "err_exp", "freenrg_fep", "err_fep"]).transpose()

# save our results to a file that can be opened in e.g. Excel.
freenrg_df.to_csv("outputs/fep_results_table.csv")
```

If we want to view our results we can call ```print(freenrg_df)```, which will output the following table:

*Table 1: table view of estimated free energies of binding versus experimental (all values in kcal/mol).*

![image-20210416160212921](/Users/jscheen/projects/BioSimSpaceTutorials/03_fep/inputs/tut_imgs/freenrg_analysis_table.png)

Although informative, viewing these results in a table is not very intuitive. Next, we will generate some plots common in FEP studies.

#### 3.4. Generating plots

```Matplotlib``` is a very popular plotting library in python. Although it takes a bit of time to learn the basics, the library offers great flexibility in plotting any kind of data when compared to e.g. MS Excel. 

First, we will generate a barplot of our predictions, complete with error bars and plot formatting.

```python
# initiate an empty figure with fixed dimensions.
fig, ax = plt.subplots(figsize=(10,5))

# determine positions for X axis labels.
x_locs = np.arange(len(freenrg_df))

# set bar width
width = 0.35  

# plot both our experimental and FEP free energies using an offset on the x position so bars don't overlap.
ax.bar(x_locs - width/2, height=freenrg_df["freenrg_exp"], width=width, yerr=freenrg_df["err_exp"],
                label='Experimental')
ax.bar(x_locs + width/2, height=freenrg_df["freenrg_fep"], width=width, yerr=freenrg_df["err_fep"],
                label='FEP')
 
# format the plot further.
plt.axhline(color="black")
plt.ylabel("$\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
plt.xticks(x_locs, freenrg_df.index, rotation=70, ha="right")
plt.legend()

plt.savefig("outputs/fep_vs_exp_barplot.png", dpi=300)
```

The plot will now be saved as a high-resolution image in ```./outputs/``` and will look like this:

![img](/Users/jscheen/projects/BioSimSpaceTutorials/03_fep/inputs/tut_imgs/barplot.png)

*Figure 1: barplot of free energies of binding predicted by FEP (orange) versus experimental (blue). All values depicted are in kcal/mol.*

############## REPLACE FIG WITH WLS FIG

Another common way of depicting FEP predictions is with scatterplots. FEP scatterplots are often annotated with 1- and 2 kcal/mol confidence regions to indicate regions in which a FEP prediction is deemed highly accurate (inner band) and reasonably accurate (outer band).

```python
plt.figure(figsize=(10,10))

plt.scatter(freenrg_df["freenrg_exp"], freenrg_df["freenrg_fep"], zorder=10)

# plot 1/2 kcal bounds:
plt.fill_between(
				x=[-15, 15], 
				y2=[-14.75,15.25],
				y1=[-15.25, 14.75],
				lw=0, 
				zorder=-10,
				alpha=0.3,
				color="grey")
# upper bound:
plt.fill_between(
				x=[-15, 15], 
				y2=[-14.5, 15.5],
				y1=[-14.75, 15.25],
				lw=0, 
				zorder=-10,
				color="grey", 
				alpha=0.2)
# lower bound:
plt.fill_between(
				x=[-15, 15], 
				y2=[-15.5,14.5],
				y1=[-14.75, 15.25],
				lw=0, 
				zorder=-10,
				color="grey", 
				alpha=0.2)

# plot error bars:
yerr = freenrg_df["err_fep"]
xerr = freenrg_df["err_exp"]

plt.errorbar(freenrg_df["freenrg_exp"], freenrg_df["freenrg_fep"], 
            yerr=yerr,
            xerr=xerr,   # comment this line to hide experimental error bars \
                         # as this can sometimes overcrowd the plot.
            ls="none",
            lw=0.5, 
            capsize=2,
            color="black",
            zorder=5
            )

# format the plot further.
plt.axhline(color="black", zorder=1)
plt.axvline(color="black", zorder=1)
plt.ylabel("Predicted $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
plt.xlabel("Experimental $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")

# get the bounds. This can be done with min/max or simply by hand.
all_freenrg_values = np.concatenate([freenrg_df["freenrg_exp"].values,freenrg_df["freenrg_fep"].values])
min_lim = min(all_freenrg_values)
max_lim = max(all_freenrg_values)

# for a scatterplot we want the axis ranges to be the same. 
plt.xlim(min_lim*1.3, max_lim*1.3)
plt.ylim(min_lim*1.3, max_lim*1.3)

plt.savefig("outputs/fep_vs_exp_scatterplot.png", dpi=300)
```

The plot will now be saved as a high-resolution image in ```./outputs/``` and will look like this:

![img](/Users/jscheen/projects/BioSimSpaceTutorials/03_fep/inputs/tut_imgs/scatterplot.png)

*Figure 2:  scatterplot of free energies of binding predicted by FEP (y axis) versus experimental (x axis). All values depicted are in kcal/mol. Depicted confidence bounds are of 1 kcal/mol (dark grey) and 2 kcal/mol (light grey).*

############## REPLACE FIG WITH WLS FIG

#### Final remarks

There are several other things we might do with ```FreeEnergyNetworkAnalysis```. 

For instance, in advanced stages of lead optimisation we might be working with a ligand (for which we have experimental affinity) and a series of new ligands (for which we have none). In such a case, it makes sense to use that ligand as a reference. This can be done simply by replacing the above code with 

```python
# Add a reference compound if needed. 
target_compound = 'ejm48'
pG.compute_weighted_avg_paths(target_compound)

computed_relative_DDGs = pG.freeEnergyInKcal

pG.write_free_energies(computed_relative_DDGs)
```

Sometimes, when working with ligand series that contain very dissimilar ligands, it is necessary to include one or more *intermediate* ligands. The purpose of these is to create more reliable transformations (at the cost of doing extra simulations). For analysis these can be discarded, as we are not necessarily interested in their binding affinity to the protein target. These can be excluded from the analysis graph by running:

```python
pG.format_free_energies(intermed_ID="name_of_intermediate")
```

*this should be run after computing free energies and before writing free energies, such that the intermediate energy is excluded.*

