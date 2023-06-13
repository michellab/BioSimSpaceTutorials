# Funnel metadynamics

Author: Dominykas Lukauskis

Email: dominykas.lukauskis.19@ucl.ac.uk

Requirements:
* BioSimSpace
* OpenMM or Gromacs/Amber patched with PLUMED

This workshop tutorial walks you through how to setup funnel metadynamics (fun-metaD) simulations to estimate the absolute binding free energy (ABFE), how to visualise the funnel restraints and how to analyse the simulations using [BioSimSpace](https://biosimspace.org/index.html).

The workshop is divided into three parts, each linked below:

## 1. Introduction and some theory of fun-metaD
As an introduction, look [here](01_Introduction.md) for an explanation of funnel metadynamics, and consult [this Jupyter NB](01_bss-fun-metaD-introduction.ipynb).

## 2. Simulation setup and restraint visualisation
For a step-by-step guide on how to visualise the funnel restraints and how to set up the fun-metaD simulations, follow [this Jupyter NB](02_bss-fun-metad-tutorial.ipynb).

## 3. Analysis
Now that the simulations have finished, we can have a look at how to analyse them [here](03_bss-fun-metad-analysis.ipynb), along with a discussion on how to think about convergence in fun-metaD simulations and how to get a robust ABFE estimate. 
