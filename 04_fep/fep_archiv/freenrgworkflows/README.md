[![PyPI version](https://badge.fury.io/py/freenrgworkflows.svg)](https://badge.fury.io/py/freenrgworkflows)
[![Build Status](https://travis-ci.org/michellab/freenrgworkflows.svg?branch=devel)](https://travis-ci.org/michellab/freenrgworkflows)
[![codecov](https://codecov.io/gh/michellab/freenrgworkflows/branch/devel/graph/badge.svg)](https://codecov.io/gh/michellab/freenrgworkflows)

# FreeEnergyNetworkAnalysis

## What is FreeEnergyNetworkAnalysis?


It is a python (3.5+) module for analysing alchemical free energy simulations based on a given pertubration map. 
The tool works well with [SOMD](https://github.com/michellab/Sire) output.

## Installation
with pip:   
   `pip install freenrgworkflows`

from github:  
 
```
git clone git@github.com:michellab/freenrgworkflows.git   
cd freenrgworkflows   
python setup.py install  
```

## Testing your installation

To test the installation from the github repository you can use the pytest frame work by typing:
```
pytest --verbosity 1 tests
```
in the root directory, i.e. freenrgworkflows

## Getting started
Networkanalysis provides both a simple API and a commandline tool in order to extract relative free energies
from a large number of perturbations. The main script can be found in the bin directory and run in the following way:
```
python run_networkanalysis.py --help
```
Gives information on the command line options available. 
An example usage can loook like this:
```
python run_networkanalysis.py ../tests/io/summary_r1.csv --target_compound FXR45 -o ~/Desktop/blub.dat -e ../tests/io/ic50_exp.dat --stats --generate_notebook
```
summary_r1.csv is an collated output of analyse_freenrg mbar from [SOMD](https://github.com/michellab/Sire), the target compound is randomly selected. The output file can be anything really, such as blub.dat written to the computers desktop. If you want to compare experiments to simulations you will need to provide a file containing IC50s on the commandline or kDs with the API. The flag --generate_notebook will generate a default jupyter notebook running the analysis from a notebook including plots comparing the experimental and computed free energies. The generated notebook can be opened using `jupyter notebook ~/Desktop/blub.ipynb`

## Authors
- Antonia Mey [antonia.mey@ed.ac.uk]

## Questions, Comments, Docs

Contact one of the authors.
