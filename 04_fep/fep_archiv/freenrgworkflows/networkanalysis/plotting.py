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

import matplotlib.pylab as plt
import matplotlib
import seaborn as sns
import copy
import numpy as np

# sns.set_style("ticks")
# sns.set_context("notebook", font_scale=2)


class FreeEnergyPlotter(object):
    """FreeEnergyPlotter contains helper function that allows to plot two or more free energy estimates against each other"""

    def __init__(self, DDG_series1, DDG_series2, compound_list=None):
        r"""
        Parameters
        ----------
        DDG_series1 : list of dictionaries
            list of dictionaries containing free energies
        DDG_series2 : list of dictionaries
            list of second set of dictionaries containing free energies
        compound_list : list
            list containg compound names for which to plot free energies. 

        """

        self.dataseries1 = []
        self.dataseries2 = []
        self.labels = []
        self.compound_list = []

        cl_exp = set().union(*(d.keys() for d in DDG_series1))
        cl_comp = set().union(*(d.keys() for d in DDG_series2))
        ids = list(set(cl_exp).intersection(cl_comp))
        if 'error' in ids:
            index = ids.index('error')
            ids.pop(index)
        if compound_list is None:
            print (np.sort(ids))
            self.compound_list = np.sort(ids)
        else:
            if not set(compound_list).issubset(ids):
                print(
                    "The compound list you have supplied does not match the compounds in the compounds in the dataseries")
                sys.exit(1)
            self.compound_list = compound_list

        for e in self.compound_list:
            s1 = False
            s2 = False
            for i in DDG_series1:
                if e in i:
                    data = []
                    data.append(i[e])
                    data.append(i['error'])
                    self.dataseries1.append(data)
                    s1 = True
            for i in DDG_series2:
                if e in i:
                    data = []
                    data.append(i[e])
                    data.append(i['error'])
                    self.dataseries2.append(data)
                    s2 = True
            if s1 and s2:
                self.labels.append(e)
        self.dataseries2 = np.array(self.dataseries2)
        self.dataseries1 = np.array(self.dataseries1)

    def plot_bar_plot(self, legend=('experimental', 'computed'),
                      colors=[sns.xkcd_rgb["pale red"], sns.xkcd_rgb["denim blue"]]):
        r"""Do a barplot of two different free energy series
        Parameters
        ----------
        legend : list
            contains the labels for the two dataseries, e.g. experimental and computed
        colors : list
            list of second set of dictionaries containing free energies
        """

        N = len(self.dataseries1) * 2

        ind = np.arange(0, N, 2)  # the x locations for the groups
        width = 0.35 * 2  # the width of the bars
        fig, ax = plt.subplots(figsize=(8, 6))

        rects1 = ax.bar(ind, self.dataseries1[:, 0], width, yerr=self.dataseries1[:, 1], color=colors[0], ecolor='k')
        rects2 = ax.bar(ind + width, self.dataseries2[:, 0], width, yerr=self.dataseries2[:, 1], color=colors[1],
                        ecolor='k')

        # add some text for labels, title and axes ticks
        ax.set_ylabel(r'$\Delta \Delta G$ in [kcal/mol]', fontsize=15)
        ax.set_xticks(ind + width)
        ax.set_xticklabels(self.labels, fontsize=15, rotation='vertical')

        ax.legend((rects1[0], rects2[0]), legend, fontsize=15, loc='best')
        sns.despine()
        return ax, fig

    def plot_scatter_plot(self, xlabel=r'experimental $\Delta \Delta G$ in [kcal/mol]',
                          ylabel=r'computed $\Delta \Delta G$ in [kcal/mol]', color=sns.xkcd_rgb["pale red"]):
        r"""Do a barplot of two different free energy series
        Parameters
        ----------
        colors : string
            color string for matplot lib for the data color
        """

        # TODO:
        # Add args for plotting in general rather than setting exact defaults.
        plt.errorbar(self.dataseries1[:, 0], self.dataseries2[:, 0], yerr=self.dataseries2[:, 1],
                     color=sns.xkcd_rgb["pale red"], linewidth=0, elinewidth=1, marker='o')
        a = self.dataseries1[:, 0]
        b = self.dataseries2[:, 0]
        all_min = min(min(a), min(b))
        all_max = max(max(a), max(b))
        x = np.linspace(all_min, all_max, 10)
        plt.plot(x, x, '--')
        plt.xlim(np.min(self.dataseries1[:, 0]) - 1, np.max(self.dataseries1[:, 0]) + 1)
        plt.ylim(np.min(self.dataseries2[:, 0]) - 1, np.max(self.dataseries2[:, 0]) + 1)
        # add some text for labels, title and axes ticks
        plt.ylabel(ylabel, fontsize=15)
        plt.xlabel(xlabel, fontsize=15)

        sns.despine()

    def _plot_master_plot(self, x, y, title='',
                          xerr=None, yerr=None,
                          method_name='', target_name='', quantity=rf'$\Delta \Delta$ G',
                          xlabel=f'Experimental', ylabel=f'Calculated', units=r'$\mathrm{kcal\,mol^{-1}}$',
                          guidelines=True, origins=True, color=None,
                          statistics=['RMSE', 'MUE'], filename=None, centralizing=True,
                          shift=0., figsize=6):
        """ Handles the aesthetics of the plots in one place.
        Parameters
        ----------
        x : list
            Values to plot on the x axis
        y : list
            Values to plot on the y axis
        title : string, default = ''
            Title for the plot
        xerr : list , default = None
            Error bars for x values
        yerr : list , default = None
            Error bars for y values
        method_name : string, optional
            name of method associated with results, e.g. 'perses'
        target_name : string, optional
            name of system for results, e.g. 'Thrombin'
        quantity : str, default = '$\Delta \Delta$ G'
            metric that is being plotted
        xlabel : str, default = 'Experimental'
            label for xaxis
        ylabel : str, default = 'Calculated'
            label for yaxis
        units : str, default = r'$\mathrm{kcal\,mol^{-1}}$'
            string value of units to label axis
        guidelines : bool, default = True
            toggles plotting of grey 0.5 and 1 kcal/mol error zone
        origins : bool, default = True
            toggles plotting of x and y axis
        color : str, default = None
            if None, will be coloured according to distance from unity
        statistics : list(str), default = ['RMSE',  'MUE']
            list of statistics to calculate and report on the plot
        filename : str, default = None
            filename for plot
        centralizing : bool, default = True
            ofset the free energies
        shift : float, default = 0.
            shift both the x and y axis by a constant
        figsize : float, default = 3.25
            size of figure for matplotlib
        Returns
        -------
        """
        nsamples = len(x)
        # aesthetics
        plt.rcParams['xtick.labelsize'] = 12
        plt.rcParams['ytick.labelsize'] = 12
        plt.rcParams['font.size'] = 12

        f = plt.figure(figsize=(figsize, figsize))
        plt.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.8)

        plt.xlabel(f'{xlabel} {quantity} / ' + units)
        plt.ylabel(f'{ylabel} {quantity} / ' + units)

        ax_min = min(min(x), min(y)) - 0.5
        ax_max = max(max(x), max(y)) + 0.5
        scale = [ax_min, ax_max]

        plt.xlim(scale)
        plt.ylim(scale)

        # plots x-axis and y-axis
        if origins:
            plt.plot([0, 0], scale, 'gray')
            plt.plot(scale, [0, 0], 'gray')

        # plots x=y line
        plt.plot(scale, scale, 'k:')
        if guidelines:
            small_dist = 0.5
            # plots grey region around x=y line
            plt.fill_between(scale, [ax_min - small_dist, ax_max - small_dist],
                             [ax_min + small_dist, ax_max + small_dist],
                             color='grey', alpha=0.2)
            plt.fill_between(scale, [ax_min - small_dist * 2, ax_max - small_dist * 2],
                             [ax_min + small_dist * 2, ax_max + small_dist * 2],
                             color='grey', alpha=0.2)
        # actual plotting
        cm = plt.get_cmap('coolwarm')

        if color is None:
            color = np.abs(x - y)
            # 2.372 kcal / mol = 4 RT
            color = cm(color / 2.372)
        plt.errorbar(x, y, xerr=xerr, yerr=yerr, color='gray', linewidth=0., elinewidth=2., zorder=1)
        plt.scatter(x, y, color=color, s=10, marker='o', zorder=2)

        # stats and title

        long_title = f'{title} \n {target_name} '

        plt.title(long_title, fontsize=12, loc='right', horizontalalignment='right', family='monospace')

        if filename is None:
            plt.show()
        else:
            plt.savefig(filename, bbox_inches='tight')
        return f

    def plot_DGs(self, method_name='', target_name='', title='', filename=None, **kwargs):

        self._plot_master_plot(self.dataseries1[:, 0], self.dataseries2[:, 0],
                          xerr=self.dataseries1[:, 1], yerr=self.dataseries2[:, 1], filename=filename,
                          title=title, method_name=method_name, target_name=target_name, **kwargs)

    def _plot_bar_plot_no_dic_graph2(self, r1_weight, r2_weight, r3_weight, keys):
        r""" Obsolete do not use!!!

        """
        labels = keys

        N = len(r1_weight) * 2

        ind = np.arange(0, N, 2)  # the x locations for the groups
        width = 0.25 * 2  # the width of the bars
        fig, ax = plt.subplots(figsize=(7, 7))
        print (r1_weight)
        data1 = r2_weight[:, 0]
        data2 = r3_weight[:, 0]
        rects1 = ax.bar(ind, r1_weight, width, color=sns.xkcd_rgb["pale red"])
        rects2 = ax.bar(ind + width, data1, width, yerr=r2_weight[:, 1], color=sns.xkcd_rgb["denim blue"],
                        error_kw={'ecolor': 'black',  # error-bars colour
                                  'linewidth': 2})
        rects3 = ax.bar(ind + 2 * width, data2, width, yerr=r3_weight[:, 1], color='#759D70',
                        error_kw={'ecolor': 'black',  # error-bars colour
                                  'linewidth': 2})

        # add some text for labels, title and axes ticks
        ax.set_ylabel(r'$\Delta \Delta G$ in [kcal/mol]', fontsize=20)
        ax.set_xticks(ind + width + width * 0.5)
        ax.set_xticklabels(labels, fontsize=20)

        # ax.legend((rects1[0], rects2[0]), ('experimental', 'computational'), fontsize=20)
        sns.despine()

    def _plot_hist(self, edges, hist, label, color, alpha):
        r""" Obsolete do not use!

        """
        halfwidth = (edges[1] - edges[0]) / 2
        centers = edges[:-1] + halfwidth
        fig = plt.figure(figsize=(8, 4))
        plt.plot(centers, hist, color=color)
        plt.fill_between(centers, hist, facecolor=color, alpha=alpha)
        plt.xlabel(label)
        plt.ylabel('P(%s)' % label)
        sns.despine()
