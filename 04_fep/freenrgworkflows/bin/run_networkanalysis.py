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

####################################################################################################
#
#   IMPORTS
#
####################################################################################################

from networkanalysis.networkanalysis import *
from networkanalysis.experiments import *
from networkanalysis.stats import *
from networkanalysis.jupyter import *
import networkanalysis
from argparse import ArgumentParser, FileType
import numpy as np
import os
import warnings

####################################################################################################
#
#   MAIN PART
#
####################################################################################################

if '__main__' == __name__:

    ############################################################################
    #
    #   capture the command line arguments
    #
    ############################################################################

    parser = ArgumentParser()
    parser.add_argument(
        'files',
        help='Networkx compatible csv file/files of the computed free energies with Sire',
        nargs='*',
        metavar='FILE'
    )
    parser.add_argument(
        '-o',
        '--network_output',
        help='File to write final free energies differences to',
        metavar='FILE',
        default=None
    )
    parser.add_argument(
        '-e',
        '--experiments',
        help='File containing experimental IC50 data',
        default=None,
        metavar='FILE'
    )
    parser.add_argument(
        "--target_compound",
        help="Name of the reference compound with respect to which the free energy should be computed",
        metavar='STRING'
    )
    parser.add_argument(
        "--intermed_ID",
        help="String identifier for intermediates, e.g. INT_01 which should not be retained in the final output",
        metavar='STRING'
    )
    parser.add_argument(
        "--stats",
        help="Print correclation statistics between computated and experimental data, "
             "this will only work if and experimental data file was given",
        action='store_true'
    )
    parser.add_argument(
        "--comments",
        help="Identifier used to mark comments in the input network files",
        metavar='STRING',
        default='#'
    )
    parser.add_argument(
        "--delimiter",
        help="delimiter used in the input network files",
        metavar='STRING',
        default=','
    )
    parser.add_argument(
        "--merge_BM",
        help="Merge binding modes, assuming they are identified as compoun_BM1 and compound_BM2",
        metavar='BOOLEAN',
        default='True'
    )
    parser.add_argument(
        "--weighted",
        help="Compute weighted path averages when true and unweighted path averages when false",
        metavar='BOOLEAN',
        default='True'
    )
    parser.add_argument(
        "--generate_notebook",
        help="Autogenerates a jupyter notebook showing the working of the anaysis and useful plots. "
             "The filename is the arguemtn of -o with a .ipynb extension.",
        action='store_true'
    )

    args = parser.parse_args()

    ############################################################################
    #
    #   check mandatory command line arguments
    #
    ############################################################################
    if 1 > len(args.files):
        raise OSError('You must give at least one networkanalysis networkx compatible file')

    ############################################################################
    #
    #   write header
    #
    ############################################################################
    print (
                "\n\n################# NETWORKANALYSIS v. %s WITH NETWORKX ################################" % networkanalysis.__version__)
    print ("\n\n########################## Parameters ######################################")
    print ("filelist: \t\t\t\t%s" % args.files)
    print ("file comment: \t\t\t\t%s" % args.comments)
    print ("file delimiter: \t\t\t%s" % args.delimiter)
    print ("target compound: \t\t\t%s" % args.target_compound)
    print ("intermed_ID: \t\t\t\t%s" % args.intermed_ID)
    print ("Network computed free energies file: \t%s" % args.network_output)
    print ("IC50s datafile: \t\t\t%s" % args.experiments)
    print ("Weidghted averages:\t\t\t%s" % args.weighted)
    print ("Merge binding modes:\t\t\t%s" % args.merge_BM)
    print ("#############################################################################\n\n")

    # Do the network analysis
    pG = PerturbationGraph()
    pG.populate_pert_graph(args.files[0], delimiter=args.delimiter, comments=args.comments)
    if len(args.files) > 1:
        for f in args.files[1:]:
            pG.add_data_to_graph(f, delimiter=args.delimiter, comments=args.comments)
    target_compound = args.target_compound
    if target_compound == None:
        target_compound = list(pG.graph.nodes)[0]
        warnings.warn(
            UserWarning("No target compound given, using the first compound in the node list: %s" % target_compound))
    if args.weighted == False:
        pG.compute_avg_paths(target_compound)
    else:
        pG.compute_weighted_avg_paths(target_compound)
    pG.format_free_energies(merge_BM=args.merge_BM, intermed_ID=args.intermed_ID, weighted=args.weighted)
    comp_DDG = pG.freeEnergyInKcal

    if args.network_output != None:
        pG.write_free_energies(comp_DDG, filename=args.network_output)
    else:
        pG.write_free_energies(comp_DDG)

    # Read experimental data
    if args.experiments != None and args.stats:
        ex = ExperimentalData()
        ex.compute_DDG_from_IC50s(args.experiments, reference=target_compound)
        exp_DDG = ex.freeEnergiesInKcal
        stats = freeEnergyStats()
        stats.generate_statistics(comp_DDG, exp_DDG, repeats=1000)

        print("\n########################## Statistics ######################################")
        print(" R and std = %f +/- %f" % (stats.R_mean, stats.R_std))
        print(" R2 and std = %f +/- %f" % (stats.R2_mean, stats.R2_std))
        print(" tau and std = %f +/- %f" % (stats.tau_mean, stats.tau_std))
        print(" MUE and std = %f +/- %f" % (stats.mue_mean, stats.mue_std))
        print("#############################################################################\n\n")

    if args.generate_notebook:
        try:
            JupyterNotebookCreator
        except NameError:
            warnings.warn(UserWarning("The Jupyter notebook module is not available, "
                "so generating a notebook is not possible"))
            args.generate_notebook = False
        else:
            if args.network_output != None:
                nbname = os.path.splitext(args.network_output)[0] + '.ipynb'
                print(nbname)
            else:
                nbname = "Default_Analysis.ipynb"
            print("\n###########################Generating jupyter notebook#######################")
            book = JupyterNotebookCreator(nbname, networkfile=args.files[0], experimentalfile=args.experiments)
            book.write_notebook()
            print("#                       Notebook written to %s" % nbname)
            print("##############################################################################\n\n")

    ############################################################################
    #
    #   say good bye
    #
    ############################################################################
    print("\n#################################################################################################\n#")
    print("#                 That's it, now it's time to put the kettle on ")
    print("#                Thank you for using the network analysis package!")
    print("#\n################################################################################################\n\n")
