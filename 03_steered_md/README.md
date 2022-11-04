# Steered molecular dynamics

Author: Adele Hardie

Email: adele.hardie@ed.ac.uk

#### Requirements:
* BioSimSpace
* AMBER or GROMACS compiled with PLUMED
* An equilibrated starting system
* A target structure

This tutorial covers how to set up and run steered MD (sMD) simulations with BioSimSpace. In this example, the sMD trajectory is used to obtained seeds for equilibirum MD simulations for building a Markov State Model (MSM). The overall protocol consists of 3 main steps:

### 1. Set up and run sMD
An example of how to set up and run sMD with BioSimSpace is provided in [01_setup_sMD](01_setup_sMD.md). It covers some background for sMD, and creating steering protocols via BioSimSpace. Examples of both a single CV and a multi CV protocol are provided. A [notebook version](01_setup_sMD_AMBER.ipynb) exists as well. Additionally, a [python script](01_run_sMD_simple.py) that can be run from command line is provided. Note that it is specific to using RMSD of all heavy atoms of specified residues as the CV, but it can be easily modified to work with other CVs [available in BioSimSpace](https://biosimspace.org/api/index_Metadynamics_CollectiveVariable.html). Example [LSF](01_run_sMD_LSF.sh) and [slurm](01_run_sMD_slurm.sh) submission scripts are available.

### 2. Analysing trajectory and seeded MD
Once an sMD run is complete, the trajectory can be analysed and snapshots for seeded MD extracted [here](02_trajectory_analysis.ipynb). In this case we only look at the `COLVAR` file produced by PLUMED, but, depending on the system, other criteria might be important to decide whether the sMD run was successful. After the snapshots are extracted, they are used as starting points for equilibrium MD simulations, the basics of which are covered in the [introduction](../01_introduction). Note that the [python script](02_run_seededMD.py) is specific to the example system of PTP1B due to additional parameters required, but the `load_system()` function can be easily modified to work woth other systems. Since this step is very computationally intensive, it is recommended to do this on an HPC system. An example array job [submission script](02_run_seeded_MD_LSF.sh) is available.

### 3. Build the MSM
After all the data for the MSM is gathered via seeded MD simulations, it is time to build it. Here we use [PyEMMA](http://emma-project.org/latest/), a python library with their own [tutorials](http://emma-project.org/latest/tutorial.html). There is a lot to MSMs that could span a few workshops on its own, but for illustration we provide [results](03_msm.md) for PTP1B. For those who are curious, a [full notebook](03_msm_full.ipynb) describes exactly how that particular MSM was built, but there are many ways of doing it and so we refer you to PyEMMA's website for more detail.
