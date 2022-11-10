# Steered molecular dynamics

### Authors
* Adele Hardie (adele.hardie@ed.ac.uk)

### Table of Contents
1. [Introduction](#1-Introduction)

## 1 Introduction

Allosteric inhibition can be a useful alternative to conventional protein targeting when the nature of the active site makes it difficult to design binders. One of such proteins is protein tyrosine phosphatase 1B (PTP1B), which will be used as an example system for this tutorial. The activity of PTP1B depends on the conformation of its WPD loop, which can be open (yellow) or closed (red):

<img src="figures/open-close.png" width=300>

PTP1B is difficult to drug due to the charged nature of its active site, and so is of interest for allosteric inhibition studies. However, with allostery, knowing that a molecule binds to the protein is  not enough. It also requires  assessment  of  whether  an  allosteric  binder  actually  has  an  effect  on protein function. This can be assessed by Markov State Models(MSMs). MSMs are transition matrixes that provide insight into the statistical ensemble of protein conformations, such as what is the probability that a protein will exist in a certain conormation. Comparing the probability of the target protein being catalytically active between models with and without an allsoteric binder indicates whether or not it has potential as an inhibitor. Furthermore, since the system is treated as memoryless, model building only requires local equilibrium. Therefore, it can make use of shorter MD simulations, allowing them to be run in parallel.

In order to have a more complete view of the protein ensemble, enhanced sampling methods are used, among them steered MD (sMD) (1). It introduces a bias potential that is added to the Hamiltonian, thus biasing the simulation towards a specified value of a chosen collective variable. Once the system has reached a certain conformation, those coordinates (2) can be used as starting points for equilibrium MD simulations (4) that can subsequently be used as data for constructing an MSM (4). An example summary of this is shown below:

<img src="figures/ensemble-md-protocol.png" width=450>

PLUMED is a library that, among other things, has enhanced sampling algorithms. It works with multiple MD engines, including GROMACS and AMBER. PLUMED uses a [moving restraint](https://www.plumed.org/doc-v2.5/user-doc/html/_m_o_v_i_n_g_r_e_s_t_r_a_i_n_t.html) that is calculated as follows:

$V(\vec{s},t) = \frac{1}{2} \kappa(t) ( \vec{s} - \vec{s}_0(t) )^2$     (Eq. 1)

where $\vec{s}_0$ and $\kappa$ are time dependent and specified in the PLUMED input. $\vec{s}_0$ is the target CV value and $\kappa$ is the force constant in kJ mol$^{-1}$. The values of both of them are set at specific steps, and linearly interpolated in between.

This tutorial focuses on running the prerequisite MD simulations using BioSimSpace.

## 2 Set up and run sMD
An example of how to set up and run sMD with BioSimSpace is provided in [jupyter notebook](02_setup_sMD.ipynb). It covers some background for sMD, and creating steering protocols via BioSimSpace. Examples of both a single CV and a multi CV protocol are provided. Additionally, a example scripts that can be run from command line is provided, both for [single CV](scripts/sMD_simple.py) and [multiple CV](scripts/sMD_multiCV.py). Note that it is specific to using RMSD of all heavy atoms of specified residues as the CV, but it can be easily modified to work with other CVs [available in BioSimSpace](https://biosimspace.org/api/index_Metadynamics_CollectiveVariable.html). Example [LSF](scripts/sMD_LSF.sh) and [slurm](scripts/sMD_slurm.sh) submission scripts are available.

## 3 Analysing trajectory and seeded MD
Once an sMD run is complete, the trajectory can be analysed and snapshots for seeded MD extracted [here](03_trajectory_analysis.ipynb). In this case we only look at the `COLVAR` file produced by PLUMED, but, depending on the system, other criteria might be important to decide whether the sMD run was successful. After the snapshots are extracted, they are used as starting points for equilibrium MD simulations, the basics of which are covered in the [introduction](../01_introduction). Note that the [python script](scripts/seededMD.py) is specific to the example system of PTP1B due to additional parameters required, but the `load_system()` function can be easily modified to work with other systems. Since this step is very computationally intensive, it is recommended to do this on an HPC system. An example array job [submission script](scripts/seededMD_LSF.sh) is available.

## 4 Build the MSM
After all the data for the MSM is gathered via seeded MD simulations, it is time to build it. Here we use [PyEMMA](http://emma-project.org/latest/), a python library with their own [tutorials](http://emma-project.org/latest/tutorial.html). For those who are curious, a [full notebook](04_msm.ipynb) describes exactly how the MSM for PTP1B was built, but there are many ways of doing it and so we refer you to PyEMMA's website for more detail.

#### Requirements:

* BioSimSpace
* AMBER or GROMACS compiled with PLUMED
* An equilibrated starting system
* A target structure
* pyemma (optional, to build MSMs)
* pytraj (optional, to featurise data for MSMs)


