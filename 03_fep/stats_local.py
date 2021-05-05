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


__author__ = "Antonia Mey", "Julien Michel"
__email__ = "antonia.mey@ed.ac.uk"

import numpy as np
import networkx as nx
import scipy.stats
import copy
import warnings


class freeEnergyStats(object):
    """docstring for freeEnergyStats"""

    r"""
    Parameters
    ----------

    prediction : list of dictionaries
        list of dictionaries of computed free energies and their errors
        Default = None
    target : list of dictionaries
        list of dictionaries of experimental free energies and their errors
        Default = None
    compound_list : list of strings
        list should contain dictionary keys of compounds to be compared statistically
        Default = None
    verbose : boolean
        Be loud and noisy!
        Default = False
    """
    def __init__(self, prediction=None, target=None, compound_list=None, verbose=False):

        self.data_comp = None
        self.data_exp = None
        self._compound_list = None
        if prediction is not None:
            if target is None:
                raise ValueError("If you give a prediction you also need to give a target value")
            else:
                if compound_list is None:
                    cl_exp = set().union(*(d.keys() for d in target))
                    cl_comp = set().union(*(d.keys() for d in prediction))
                    compound_list = list(set(cl_exp).intersection(cl_comp))
                    if 'error' in compound_list:
                        index = compound_list.index('error')
                        compound_list.pop(index)
                    self._compound_list = compound_list
                else:
                        self._compound_list = compound_list
            for k in self._compound_list:
                comp = next(item for item in prediction if k in item)
                exp = next(item for item in target if k in item)
                val = comp[k]
                err = comp['error']
                self.data_comp.append([val, err])
                val = exp[k]
                self.data_exp.append(val)
            if verbose:
                print("Setup stats...")

        # Statistics from bootstrapping
        self._PI = None
        self._R = None
        self._R2 = None
        self._tau = None
        self._mue = None
        self._rmse = None
        self._rmse_error = None
        self._R_error = None
        self._R2_error = None
        self._tau_error = None
        self._mue_error = None
        self._confidence_interval = 0.68

        # Direct statistics from data
        self._R_from_data = None
        self._mue_from_data = None
        self._tau_from_data = None
        self._rmse_from_data = None

    def generate_statistics(self, comp_data, exp_data, compound_list=None, repeats=1000):
        r"""
        Parameters
        ----------

        comp_data : list of dictionaries
            list of dictionaries of computed free energies and their errors
        exp_data : list of dictionaries
            list of dictionaries of experimental free energies and their errors
        compound_list : list of strings
            list should contain dictionary keys of compounds to be compared statistically
        repeats : integer
            number of times new samples are drawn from the gaussian distribution of the computational data
        """
        self.data_comp = []
        self.data_exp = []

        # Setting up compound list
        if compound_list is None:
            cl_exp = set().union(*(d.keys() for d in exp_data))
            cl_comp = set().union(*(d.keys() for d in comp_data))
            compound_list = list(set(cl_exp).intersection(cl_comp))
            if 'error' in compound_list:
                index = compound_list.index('error')
                compound_list.pop(index)
            self._compound_list = compound_list
        else:
            self._compound_list = compound_list

        # Setting up computed and experimental dictionaries
        for k in self._compound_list:
            comp = next(item for item in comp_data if k in item)
            exp = next(item for item in exp_data if k in item)
            val = comp[k]
            err = comp['error']
            self.data_comp.append([val, err])
            val = exp[k]
            self.data_exp.append(val)

        new_data = np.array(self.data_comp)[:,0]
        self._R_from_data, p = scipy.stats.pearsonr(new_data, np.array(self.data_exp))
        self._tau_from_data = scipy.stats.kendalltau(new_data, np.array(self.data_exp))[0]
        self._rmse_from_data = self._calculate_rmse(new_data, np.array(self.data_exp))
        self._mue_from_data = self._calculate_mue(new_data, np.array(self.data_exp))

        self._R = []
        self._R2 = []
        self._tau = []
        self._mue = []
        self._rmse = []

        # Now generate the data
        for i in range(repeats):
            new_data = []
            for i in range(len(self.data_comp)):
                val = self.data_comp[i][0]
                err = self.data_comp[i][1]
                if err != 0.0:
                    val2 = np.random.normal(val, err)
                    new_data.append(val2)
                else:
                    new_data.append(val)
            R2, R = self._calculate_r2(new_data, self.data_exp)
            tau = self._calculate_tau(new_data, self.data_exp)
            mue = self._calculate_mue(new_data, self.data_exp)
            rmse = self._calculate_rmse(np.array(new_data), np.array(self.data_exp))
            self._R.append(R)
            self._R2.append(R2)
            self._tau.append(tau)
            self._mue.append(mue)
            self._rmse.append(rmse)

    def _calculate_predictive_index(self, prediction, target):
        '''r This function needs to be implemented properly'''
        raise NotImplementedError('Calculating predictive index not implemented yet.')
        '''sumwijcij = 0.0
        sumwij = 0.0

        keys = series1.keys()
        keys.sort()

        for i in range(0,len(keys)):
            keyi = keys[i]
            for j in range(i+1,len(keys)):
                keyj = keys[j]
                wij = abs(series1[keyj][0] - series1[keyi][0] )
                # print "series0 j %s series 0 i %s wij %s i %s j %s" % (series[0][j],series[0][i],wij,i,j)
                num =  (series1[keyj][0] - series1[keyi][0])
                den =  (series2[keyj][0] - series2[keyi][0] )
                #if den < 0.0001:
                #    den = 0.001
                #print num, den
                val = num / den
                # print val,serie[j],serie[i]
                if val > 0:
                    cij = 1.0
                elif val < 0:
                    cij = -1.0
                # print cij
                sumwijcij += wij*cij
                sumwij += wij
        PI = sumwijcij/sumwij
        return PI
        '''

    def _calculate_rmse(self, prediction, target):
            return np.sqrt(np.mean(np.square(target - prediction)))

    def _calculate_r2(self, prediction, target):
        r_value, p = scipy.stats.pearsonr(prediction, target)

        return r_value ** 2, r_value

    def _calculate_tau(self, prediction, target):
        tau = scipy.stats.kendalltau(prediction, target)
        return tau[0]

    def _calculate_mue(self, prediction, target):

        sumdev = 0.0
        for x in range(0, len(prediction)):
            sumdev += abs(prediction[x] - target[x])
        sumdev /= len(prediction)

        # print sumdev
        return sumdev

    def _confidence(self, data):
        sorted_data = np.sort(data)
        lower = int(np.floor((1 - self.confidence_interval) * len(sorted_data)))
        upper = int(np.ceil(self.confidence_interval * len(sorted_data)))
        return [sorted_data[lower], sorted_data[upper]]

    @property
    def confidence_interval(self):
        return self._confidence_interval

    @confidence_interval.setter
    def confidence_interval(self, confidence_interval):
        if confidence_interval < 0 or confidence_interval > 1:
            warnings.warn(UserWarning(
                'Confidence interval needs to be between 0 and 1, please try something like 0.68 for one sigma confidence'))
        self._confidence_interval = confidence_interval

    @property
    def rmse(self):
        """
        Returns:
        -------
        rmse of original input datasets
        """
        return self._rmse_from_data

    @property
    def rmse_mean(self):
        """
        Returns:
        -------
        Mean of sampled RMSE
        """
        return np.mean(self._rmse)

    @property
    def rmse_std(self):
        """
        Returns:
        -------
        standard deviation of sampled RMSE
        """
        return np.std(self._rmse)

    @property
    def rmse_confidence(self):
        self._rmse_error = self._confidence(self._rmse)
        self._rmse_error = np.concatenate([[np.median(self._rmse)], self._rmse_error])
        return self._rmse_error

    @property
    def R(self):
        """
        Returns:
        -------
        person R of original input datasets
        """
        return self._R_from_data

    @property
    def R_mean(self):
        """
        Returns:
        -------
        mean of sampled pearson R
        """
        return np.mean(self._R)

    @property
    def R_std(self):
        """
        Returns:
        -------
        standard deviation of sampled pearson R
        """
        return np.std(self._R)

    @property
    def R_confidence(self):
        """
        Returns:
        -------
        confidence : np.array
            [median, lower_bound, upper_bound]
        """
        self._R_error = self._confidence(self._R)
        self._R_error = np.concatenate([[np.median(self._R)], self._R_error])
        return self._R_error

    @property
    def R2_mean(self):
        return np.mean(self._R2)

    @property
    def R2_std(self):
        return np.std(self._R2)

    @property
    def R2_confidence(self):
        """
        Returns:
        -------
        confidence : np.array
            [median, lower_bound, upper_bound]
        """
        self._R2_error = self._confidence(self._R2)
        self._R2_error = np.concatenate([[np.median(self._R2)], self._R2_error])
        return self._R2_error

    @property
    def tau(self):
        """
        Returns:
        -------
        Kendall tau of original input datasets
        """
        return self._tau_from_data

    @property
    def tau_mean(self):
        """
        Returns:
        -------
        mean of sampled Kendall tau
        """
        return np.mean(self._tau)

    @property
    def tau_std(self):
        """
        Returns:
        -------
        standard deviation of sampled Kendall tau
        """
        return np.std(self._tau)

    @property
    def tau_confidence(self):
        """
        Returns:
        -------
        confidence : np.array
            [median, lower_bound, upper_bound]
        """
        self._tau_error = self._confidence(self._tau)
        self._tau_error = np.concatenate([[np.median(self._tau)], self._tau_error])
        return self._tau_error

    @property
    def mue(self):
        """
        Returns:
        -------
        Mean unsigned of original input datasets
        """
        return self._mue_from_data

    @property
    def mue_mean(self):
        """
        Returns:
        -------
        mean of sampled Mean unsigned error
        """
        return np.mean(self._mue)

    @property
    def mue_std(self):
        """
        Returns:
        -------
        standard deviation of sampled mean unsigned error
        """
        return np.std(self._mue)

    @property
    def mue_confidence(self):
        """
        Returns:
        -------
        confidence : np.array
            [median, lower_bound, upper_bound]
        """
        self._mue_error = self._confidence(self._mue)
        self._mue_error = np.concatenate([[np.median(self._mue)], self._mue_error])
        return self._mue_error
