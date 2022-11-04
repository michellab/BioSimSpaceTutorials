#!/usr/bin/env python
# -*- coding: utf-8 -*-

# !/usr/bin/env python

# This file is part of freenrgworkflows.
#
# Copyright 2016,2017 Julien Michel Lab, University of Edinburgh (UK)
#
# freenrgworkflows is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

__author__ = "Antonia Mey"
__email__ = "antonia.mey@ed.ac.uk"

import numpy as np
import warnings
import pandas as pd
import networkx as nx
import scipy.stats
import copy
import sys


class ExperimentalData(object):
    """docstring for ExperimentalData"""

    def __init__(self, temperature=300.0):
        self._DG_in_kcal = None
        self._DG_in_kJ = None
        self._ic50s = None
        self._kD = None
        self._referenceCompound = None
        self._kTkcal = 0.0019872041 * temperature
        self._kTkJ = 0.0083144621 * temperature
        self._RTkJ = 8.314459848 * temperature
        self._RTkcal = 1.987203611 * temperature
        self._keys = None

    def compute_affinities(self, filename, data_type=None, comments=None, delimiter=','):
        r"""


        Parameters
        ----------
        filename : path
            Path to file with experimental data
            Data of type:
            Compound name, affinity, error

        data_type: string
            type of experimental data provided chose from: ['k_D','IC50', 'pIC50', 'k_i']

        Returns
        -------

        """
        if data_type not in ['k_D', 'IC50', 'pIC50', 'K_i']:
            raise ValueError(
                'The data_type you have specified is not recognised please choose from: k_D, IC50, pIC50, K_i')
        data = pd.read_csv(filename, comment=comments, delimiter=delimiter)
        error_list = []
        error_list_kj = []
        error = None
        self._keys = []
        self._ic50s = []
        self._DG_in_kcal = []
        self._DG_in_kJ = []

        if data_type == 'IC50':
            if len(data.columns) == 2:
                error = self._kTkcal * np.log(2)
                error_kj = self._kTkJ * np.log(2)
            dG = []
            dG_kj_mol = []
            for i in range(len(data)):
                r = float(data.iloc[i, 1] / data.iloc[0, 1])
                dG.append(self._kTkcal * np.log(r))
                dG_kj_mol.append(self._kTkJ * np.log(r))
                if error is None:
                    error_list.append(data.iloc[i, 2])
                    error_list_kj.append(data.iloc[i, 2])
                else:
                    error_list.append(error)
                    error_list_kj.append(error_kj)
            # Now we normalise the free energies
            dG = dG - np.mean(dG)
            dG_kj_mol = dG_kj_mol - np.mean(dG_kj_mol)
            for i in range(len(data)):
                self._DG_in_kcal.append({data.iloc[i, 0]: dG[i], 'error': error_list[i]})
                self._DG_in_kJ.append({data.iloc[i, 0]: dG_kj_mol[i], 'error': error_list_kj[i]})

        elif data_type == 'k_D':
            if len(data.columns) == 2:
                error = self._kTkcal * np.log(2)
                error_kj = self._kTkJ * np.log(2)

            dG = []
            dG_kj_mol = []
            for i in range(len(data)):
                r = float(data.iloc[i, 1] / data.iloc[0, 1])
                dG.append(self._kTkcal * np.log(r))
                dG_kj_mol.append(self._kTkJ * np.log(r))
                if error is None:
                    error_list.append(data.iloc[i, 2])
                    error_list_kj.append(data.iloc[i, 2])
                else:
                    error_list.append(error)
                    error_list_kj.append(error_kj)
            # Now we normalise the free energies
            dG = dG - np.mean(dG)
            dG_kj_mol = dG_kj_mol - np.mean(dG_kj_mol)
            for i in range(len(data)):
                self._DG_in_kcal.append({data.iloc[i, 0]: dG[i], 'error': error_list[i]})
                self._DG_in_kJ.append({data.iloc[i, 0]: dG_kj_mol[i], 'error': error_list_kj[i]})

        elif data_type == 'pIC50':
            raise NotImplementedError
        elif data_type == 'k_i':
            raise NotImplementedError

    def compute_DDG_from_IC50s(self, filename, reference=None, smiles_string=False):
        r"""
        filename : string
            file containing ic50 data, format - compound name, ic50 value, error
        """
        warnings.warn(
            'compute_DDG_from_IC50s is deprecated use compute_affinites instead.',
            DeprecationWarning, stacklevel=2)
        self._keys = []
        self._ic50s = []
        self._DG_in_kcal = []
        self._DG_in_kJ = []
        f = open(filename, 'r')
        for line in f.readlines():
            curr_ic50 = {}
            fields = line.split(',')
            curr_ic50[fields[0]] = float(fields[1].strip())
            if smiles_string and len(fields) > 2:
                curr_ic50['smiles'] = fields[2].strip()
            # note down the keys
            self._keys.append(fields[0])
            self._ic50s.append(curr_ic50)  # append to list of ic50 compounds.
        f.close()

        reference_index = 0
        self._referenceCompound = self._keys[reference_index]
        if reference is not None:
            self._referenceCompound = reference
            reference_index = self._keys.index(self._referenceCompound)
        else:
            self._referenceCompound = self._keys[0]

        for k in range(len(self._keys)):
            key = self._keys[k]
            ic50 = self._ic50s[k][key]
            r = float(ic50 / float(self._ic50s[reference_index][self._referenceCompound]))
            a_kcal = {}
            a_kcal[key] = self._kTkcal * np.log(r)
            a_kcal['error'] = self._kTkcal * np.log(2)
            self._DG_in_kcal.append(a_kcal)
            a_kJ = {}
            a_kJ[key] = self._kTkJ * np.log(r)
            a_kJ['error'] = self._kTkJ * np.log(2)
            self._DG_in_kJ.append(a_kJ)

    def compute_DDG_from_kD(self, filename, reference=None, delimiter=','):
        r"""Reads KDs from file and converts them to DDG values to a given reference compound
        Parameters:
        ----------
        filename : string
            The name of the file containing kD value in format compound_name,kD_value
        reference : string
            Name of the reference compound which will be used to compute DDGs for
        """
        warnings.warn(
            'compute_DDG_from_kD is deprecated use compute_affinites instead.',
            DeprecationWarning, stacklevel=2)
        self._kD = []
        self._DG_in_kcal = []
        self._DG_in_kJ = []
        self._keys = []
        with open(filename, 'r') as f:
            lines = f.readlines()
            for line in lines:
                fields = line.strip().split(delimiter)
                curr_kD = {}
                curr_kD[fields[0]] = float(fields[1])
                # note down the keys
                self._keys.append(fields[0])
                self._kD.append(curr_kD)  # append to list of ic50 compounds.

        if reference is not None:
            self._referenceCompound = reference
            reference_index = self._keys.index(self._referenceCompound)
        else:
            self._referenceCompound = self._keys[0]
            reference_index = 0

        for k in range(len(self._keys)):
            key = self._keys[k]
            kD = self._kD[k][key]
            r = float(kD / float(self._kD[reference_index][self._referenceCompound]))
            a_kcal = {}
            a_kcal[key] = self._RTkcal * np.log(r)
            a_kcal['error'] = self._RTkcal * np.log(2)

            # computation for KJ/mol
            self._DG_in_kcal.append(a_kcal)
            a_kJ = {}
            a_kJ[key] = self._RTkJ * np.log(r)
            a_kJ['error'] = self._RTkJ * np.log(2)
            self._DG_in_kJ.append(a_kJ)

    def read_free_energies(self, filename, kcal=True, comment='#'):
        r"""Read free energies from a file
        filename : string
            Filename containing free energies        

        """
        self._keys = []
        self._ic50s = []
        self._DG_in_kcal = []
        self._DG_in_kJ = []
        if not kcal:
            raise (NotImplementedError('This has not been implemented yet'))
        else:
            f = open(filename, 'r')
            for line in f.readlines():
                if line.startswith(comment):
                    continue
                curr_ic50 = {}
                fields = line.split(',')
                if fields[1] == 'NoPred':
                    continue
                F_kcal = {}
                F_kcal[fields[0]] = float(fields[1])
                F_kcal['error'] = float(fields[2].strip())
                # note down the keys
                self._keys.append(fields[0])
                self._DG_in_kcal.append(F_kcal)  # append to list of ic50 compounds.
            f.close()

    @property
    def ic50s(self):
        return self._ic50s

    @property
    def kD(self):
        return self._kD

    @property
    def freeEnergiesInKcal(self):
        return self._DG_in_kcal

    @property
    def freeEnergiesInKJmol(self):
        return self._DG_in_kJ

    @property
    def compoundList(self):
        return self._keys
