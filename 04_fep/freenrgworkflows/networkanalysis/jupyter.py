#!/usr/bin/env python
# -*- coding: utf-8 -*-

# !/usr/bin/env python

# This file is part of freenrgworkflows.
#
# Copyright 2016,2017 2018 Julien Michel Lab, University of Edinburgh (UK)
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

try:
    import nbformat as nbf
except ImportError:
    nbf = None


if nbf:
    class JupyterNotebookCreator(object):
        def __init__(self, nbname, networkfile=None, experimentalfile=None, custom_heading=None):
            r"""
            Parameters:
            -----------
            nbname : string
                Name of the notebook to be generated
            networkfiles : string
                string containing network pertubration edges,
            experimentalfile : string
                filename for a file either containing either IC50s or K_d
            """
            self._notebook_name = nbname
            if networkfile is None:
                self._networkfile = ['tests/io/ic50_exp.dat']
            else:
                self._networkfile = networkfile
            if experimentalfile is None:
                self._experimentalfile = 'tests/io/ic50_exp.dat'
            else:
                self._experimentalfile = experimentalfile
            self._custom_heading = custom_heading
            pass

        def _generate_heading(self, custom_heading=None):
            if not custom_heading:
                heading = """\
    # Default Perturbation network analysis notebook
    This notebook was automatically generated using freenrgworkflows
    Author: Antonia Mey
    Email: antonia.mey@ed.ac.uk"""
            else:
                heading = custom_heading
            return nbf.v4.new_markdown_cell(heading)

        def _generate_imports(self, custom_imports=None):
            if not custom_imports:
                imports = """\
    %pylab inline
    import networkanalysis.networkanalysis as n_graph
    import networkanalysis.plotting as n_plot
    import networkanalysis.experiments as n_ex
    import networkanalysis.stats as n_stats
    import networkanalysis
    networkanalysis.__version__"""
            else:
                imports = custom_imports
            return nbf.v4.new_code_cell(imports)

        def _generate_custom_code_cell(self, code="""\ #This is code"""):
            return nbf.v4.new_code_cell(code)

        def _generate_custom_markdown_cell(self, markdown="""\ I am markdown """):
            return nbf.v4.new_markdown_cell(markdown)

        def _generate_default_notebook(self):
            nb = nbf.v4.new_notebook()
            cell_list = []
            cell_list.append(self._generate_heading())
            cell_list.append(self._generate_imports())

            pG = """\
    # Creating and populating the perturbation network
    pG = n_graph.PerturbationGraph()
    # Change the path below to the csv file containing the individual perturbations
    pG.populate_pert_graph('%s')
    # Uncomment below if you have run multiple runs for some perturbations and add file path
    #pG.add_data_to_graph('/path/to/additional/runs.csv')
    target_compound = pG.compoundList[0] #change this to your target compound
    pG.compute_weighted_avg_paths(target_compound)
    pG.format_free_energies(merge_BM=True,intermed_ID='INT')
    computed_relative_DDGs = pG.freeEnergyInKcal
    print ("Free energies computed from the perturbation network are: ")
    print ("---------------------------------------- ")
    pG.write_free_energies(computed_relative_DDGs)""" % self._networkfile

            cell_list.append(self._generate_custom_code_cell(pG))
            exp_markdown = """
    ### Experimental data
    It is useful to compare computed free energies to experimental data.
    The cells below will read in your experimental data. Just replace the path to you IC50 data in the
    `IC_50_file` variable """
            cell_list.append(self._generate_custom_markdown_cell(exp_markdown))

            exp_code = """\
    experiments = n_ex.ExperimentalData()
    IC_50_file = '%s'
    experiments.compute_DDG_from_IC50s(IC_50_file, reference=target_compound)
    experimental_DDGs = experiments.freeEnergiesInKcal
    print ("Free energies computed from IC50 data: ")
    print ("---------------------------------------- ")
    pG.write_free_energies(experimental_DDGs)""" % self._experimentalfile
            cell_list.append(self._generate_custom_code_cell(exp_code))

            plots_markdown = """
    ### Typical plots
    Below a bar plot and scatter plot template for comparing experimental and computed free energy values"""
            cell_list.append(self._generate_custom_markdown_cell(plots_markdown))

            plot_bar_code = """\
    plotter = n_plot.FreeEnergyPlotter(experimental_DDGs, computed_relative_DDGs)
    ax,fig = plotter.plot_bar_plot(legend=('experimental', 'computed'))"""
            cell_list.append(self._generate_custom_code_cell(plot_bar_code))

            plot_scatter_code = """\
    plotter.plot_scatter_plot() """
            cell_list.append(self._generate_custom_code_cell(plot_scatter_code))

            stats_markdown = """
    ### Error analysis on typical statistical measures: R_mean, MUE and Kendall tau_mean
    Below are examples of how to re-sample from the data in order to obtain error bars on correlation coefficients,
    mean unsigned errors and Kendall tau. Returned are confidence intervals of 65% and the median of the distribution.
     However, standard deviations and mean can also
    be returned, though less likely to give good information as these distributions are often heavily skewed and not
    normally distributed. """
            cell_list.append(self._generate_custom_markdown_cell(stats_markdown))

            stats_code = """\
    stats = n_stats.freeEnergyStats()
    stats.generate_statistics(computed_relative_DDGs,experimental_DDGs,repeats=10000)
    r_confidence = stats.R_confidence
    tau_confidence = stats.tau_confidence
    mue_confidence = stats.mue_confidence
    print ("R confidence is: %.2f < %.2f < %.2f" %(r_confidence[1], r_confidence[0], r_confidence[2]))
    print ("Mue confidence is: %.2f < %.2f < %.2f" %(mue_confidence[1], mue_confidence[0], mue_confidence[2]))
    print ("tau confidence is: %.2f < %.2f < %.2f" %(tau_confidence[1], tau_confidence[0], tau_confidence[2]))"""
            cell_list.append(self._generate_custom_code_cell(stats_code))
            nb['cells'] = cell_list
            return nb

        def _generate_custom_notebook(self):
            pass

        def write_notebook(self):
            nb = None
            if self._networkfile is None or self._experimentalfile is None:
                nb = self._generate_default_notebook()
            else:
                nb = self._generate_default_notebook()
            with open(self._notebook_name, 'w') as f:
                nbf.write(nb, f)
