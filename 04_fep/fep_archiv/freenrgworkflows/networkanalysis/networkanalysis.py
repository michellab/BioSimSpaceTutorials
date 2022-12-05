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
import networkx as nx
import sys
import warnings
import pandas as pd


class NetworkAnalyser(object):

    def __init__(self, use_weights=True, target_compound=None, iterations=10000, verbose=False, 
                balance_hysteresis=True):
        self._weights = {}
        self._ddG_edges = {}
        self._compoundList = None
        self._free_energies = []
        self.use_weights = use_weights
        self._nlinks = 0
        self.target_compound = target_compound
        self.iterations = iterations
        self._verbose = verbose
        self._graph = None

        # if True, use relative hysteresis penalty. If False, use absolute.
        self.balance_hysteresis = balance_hysteresis

    def _identify_header(self, path, n=5, th=0.9, comments=None):
        """

        Parameters
        ----------
        comments : String
            passing comments string along to deal with comment headers
        """
        df1 = pd.read_csv(path, header='infer', nrows=n, comment=comments)
        df2 = pd.read_csv(path, header=None, nrows=n, comment=comments)
        sim = (df1.dtypes.values == df2.dtypes.values).mean()
        columns = len(df1.columns)
        return 'infer' if sim < th else None, columns

    def read_perturbations_pandas(self, filename, delimiter=',', comments=None, source='lig_1', target='lig_2',
                                  edge_attr=['freenrg', 'error'], save_graph=True):
        """ Reads a networkx compatible csv file using pandas dataframes

        Parameters:
        -----------
        filename : path
            path to the csv file containing free energies
            Usually of the type:
            lig1,lig2,dG,ddG,engine
            a,b,-10,0.3,SOMD

        delimiter: string
            delimiter used for csv file, default is ','

        comments : string
            comment lines can be identified with a String, e.g. '#'

        source : string
            title of first column if inferred from header overridden

        target : string
            title of the second column if inferred from header overridden

        edge_atrr : list of strings
            titles of the 3rd and 4th column
        """

        header, n_cols = self._identify_header(filename, comments=comments)

        data = None
        if header is None:
            col_head = [source, target, edge_attr[0], edge_attr[1]]
            col_diff = n_cols - len(col_head)
            if col_diff == 0:
                data = pd.read_csv(filename, delimiter=delimiter, comment=comments, names=col_head)
            elif col_diff == 1:
                col_head.append('engine')
                data = pd.read_csv(filename, delimiter=delimiter, comment=comments, names=col_head)
            elif col_diff > 1:
                col_head.append('engine')
                for x in range(col_diff - 1):
                    col_head.append('not_needed_' + str(x))
                data = pd.read_csv(filename, delimiter=delimiter, comment=comments, names=col_head)
            elif col_diff < 0:
                raise ValueError("You don't have enough columns in your free energy perturbation results file")

        else:
            data = pd.read_csv(filename, delimiter=delimiter, comment=comments, header=header)

        # Getting rid of any NANs, this may make a network disconnected
        data = data.dropna()

        # print(data)

        # Now convert the pandas data frame to a networkx graph
        graph = nx.from_pandas_edgelist(data, source=source, target=target, edge_attr=edge_attr,
                                        create_using=nx.DiGraph())

        # We want to know what the largest component is, so we know if we may not be able to estimate certain free
        # energies
        largest = max(nx.strongly_connected_components(graph), key=len)

        # populate compound list:
        self._compoundList = list(graph.nodes())
        self._compoundList.sort()
        if len(largest) < len(self._compoundList):
            warnings.warn('Provided network is disconnected. Doing analysis on subgraph.') # FIXME

        # Doing analysis for all the nodes
        for node in self._compoundList:

            if node not in self._ddG_edges:
                self._ddG_edges[node] = {}
                self._weights[node] = {}
            out_edges = list(graph.out_edges(node))
            for e in out_edges:
                out_edge = e[1]
                # print(out_edge)
                if out_edge not in self._ddG_edges:
                    self._ddG_edges[out_edge] = {}
                    self._weights[out_edge] = {}
                edge_info = graph.get_edge_data(node, out_edge)

                self._ddG_edges[node][out_edge] = float(edge_info[edge_attr[0]])
                err = float(edge_info[edge_attr[1]])
                self._weights[node][out_edge] = 1 / (float(err) * float(err))
                self._nlinks += 1
            in_edges = list(graph.in_edges(node))
            for e in in_edges:
                in_edge = e[0]
                # print(in_edge)
                if in_edge not in self._ddG_edges:
                    self._ddG_edges[in_edge] = {}
                    self._weights[in_edge] = {}
                edge_info = graph.get_edge_data(in_edge, node)

                self._ddG_edges[in_edge][node] = float(edge_info[edge_attr[0]])
                err = float(edge_info[edge_attr[1]])
                self._weights[in_edge][node] = 1 / (float(err) * float(err))
                self._nlinks += 1
        if save_graph:
            self._graph = graph

    def add_data_to_graph_pandas(self, filename, delimiter=',', comments=None, source='lig_1', target='lig_2',
                                  edge_attr=['freenrg', 'error'], save_graph=True):
        r"""
        Adds data to an existing graph from a csv file using pandas dataframe

        Parameters:
        -----------
        filename : path
            path to the csv file containing free energies
            Usually of the type:
            lig1,lig2,dG,ddG,engine
            a,b,-10,0.3,SOMD

        delimiter: string
            delimiter used for csv file, default is ','

        comments : string
            comment lines can be identified with a String, e.g. '#'

        source : string
            title of first column if inferred from header overridden

        target : string
            title of the second column if inferred from header overridden

        edge_atrr : list of strings
            titles of the 3rd and 4th column
        """
        header, n_cols = self._identify_header(filename, comments=comments)

        data = None
        if header is None:
            col_head = [source, target, edge_attr[0], edge_attr[1]]
            col_diff = n_cols - len(col_head)
            if col_diff == 0:
                data = pd.read_csv(filename, delimiter=delimiter, comment=comments, names=col_head)
            elif col_diff == 1:
                col_head.append('engine')
                data = pd.read_csv(filename, delimiter=delimiter, comment=comments, names=col_head)
            elif col_diff > 1:
                col_head.append('engine')
                for x in range(col_diff - 1):
                    col_head.append('not_needed_' + str(x))
                data = pd.read_csv(filename, delimiter=delimiter, comment=comments, names=col_head)
            elif col_diff < 0:
                raise ValueError("You don't have enough columns in your free energy perturbation results file")

        else:
            data = pd.read_csv(filename, delimiter=delimiter, comment=comments, header=header)

        # Getting rid of any NANs, this may make a network disconnected
        data = data.dropna()

        # Now convert the pandas data frame to a networkx graph
        newGraph = nx.from_pandas_edgelist(data, source=source, target=target, edge_attr=edge_attr,
                                        create_using=nx.DiGraph())

        averaged_edge_counter, added_edge_counter = 0, 0

        if self._graph != None:
            for u, v, w in newGraph.edges(data=True):
                if self._graph.has_edge(u, v):

                    # compute average freenrg and propagate error.
                    z = self._graph.get_edge_data(u, v)
                    mean_edge = np.mean([z['freenrg'], w['freenrg']])
                    prop_error = 0.5 * np.sqrt(z['error'] ** 2 + w['error'] ** 2)

                    # replace the edge with new one.
                    self._ddG_edges[u][v] = mean_edge
                    self._weights[u][v] = prop_error

                    averaged_edge_counter += 1

                else:
                    self._ddG_edges[u][v] = w['freenrg']
                    self._weights[u][v] = w['error']
                    added_edge_counter += 1
        else:
            raise ValueError("No graph present to add data to. Use read_perturbations_pandas() instead.")

        print(f"Added additional data to {averaged_edge_counter} edges; added {added_edge_counter} new edges.")

    def read_perturbations(self, filename, delimiter=',', comments='#', nodetype=str,
                           data=(('weight', float), ('error', float))):
        r""" Read perturbations from networkx graph file

        Parameters
        ----------
        filename : path
            filename for computed free energies
        delimiter : string
            delimiter of csv file, Default: ','
        comments : string
            comment character, Default: '#'
        nodetype : string

        data : tuple

        Returns
        -------
        None
        """

        warnings.warn(
            'read_perturbations is deprecated, please use read_perturbations_pandas',
            DeprecationWarning, stacklevel=2)
        graph = nx.read_edgelist(filename, delimiter=delimiter, comments=comments, create_using=nx.DiGraph(),
                                 nodetype=nodetype, data=data)

        # populate compound list:
        self._compoundList = list(graph.nodes())
        self._compoundList.sort()
        if self._verbose:
            print('The graph is:')
            print(graph.nodes())
            print('done')

        f = open(filename)
        lines = f.readlines()
        for line in lines:
            if line.find(",") != -1 and not line.startswith("#"):
                mol1, mol2, nrg, err = line.split(",")
                if mol1 not in self._ddG_edges:
                    self._ddG_edges[mol1] = {}
                    self._weights[mol1] = {}
                if mol2 not in self._ddG_edges:
                    self._ddG_edges[mol2] = {}
                    self._weights[mol2] = {}
                self._ddG_edges[mol1][mol2] = float(nrg)
                self._weights[mol1][mol2] = 1 / (float(err) * float(err))
                self._nlinks += 1
        f.close()
        if self._verbose:
            print("Processed ", len(self._ddG_edges), " molecules")

    def _get_avg_nrg(self, mol1, mol2):
        """ Two molecular IDs return an average energy
        Parameters:
        -----------
        mol1 : string
            molecular ID key
        mol2 : string
            molecular ID key

        Returns:
        -------
        energy : float
            the average free energy from two IDs
        """
        eng1 = None
        try:
            eng1 = self._ddG_edges[mol1][mol2]
        except KeyError:
            eng1 = None
        try:
            eng2 = -self._ddG_edges[mol2][mol1]
        except KeyError:
            eng2 = None
        if eng1 is not None and eng2 is not None:
            return (eng1 + eng2) / 2.0
        elif eng1 is not None:
            return eng1
        return eng2

    def _compute_weight_matrix(self):
        """ Diagonal matrix containing weights that correspond to the adjacency matrix
        Returns:
        -------
        W : 2D numpy array
            Diagonal matrix containing weights
        """
        w = [1]
        for mol1 in self._ddG_edges:
            for mol2 in self._ddG_edges[mol1]:
                # Only handle links one way round: if both links present then only take one of them
                if mol1 < mol2 or mol2 not in self._ddG_edges or mol1 not in self._ddG_edges[mol2]:
                    if self.use_weights:
                        w.append(self._get_avg_weight(mol1, mol2))
                    else:
                        w.append(1.0)
        W = np.diag(w)
        return W

    def _compute_vector(self):
        """ Vector containing the pairwise DDG values
        Returns:
        -------
        b : numpy array
            array containing pairwise DDGs
        """
        b = [0]
        for mol1 in self._ddG_edges:
            for mol2 in self._ddG_edges[mol1]:
                # Only handle links one way round: if both links present then only take one of them
                if mol1 < mol2 or mol2 not in self._ddG_edges or mol1 not in self._ddG_edges[mol2]:
                    b.append(self._get_avg_nrg(mol1, mol2))
        return b

    def _compute_adjacency_matrix(self):

        firstrow = [0] * len(self._compoundList)
        firstrow[0] = 1
        A = np.array([firstrow], dtype="float64")

        for name1 in self._ddG_edges:
            for name2 in self._ddG_edges[name1]:
                # Only handle links one way round: if both links present then only take one of them
                if name1 < name2 or name2 not in self._ddG_edges or name1 not in self._ddG_edges[name2]:
                    row = [0] * len(self._compoundList)
                    row[self._compoundList.index(name1)] = -1
                    row[self._compoundList.index(name2)] = 1
                    A = np.append(A, [row], axis=0)
        return A

    def _compute_free_energies(self):
        """
            Solves the least square graph problem and populates _free_energies.

            Errors are computed using boostrapping and the
            Run 'self.iterations' iterations of bootstrap error
            analysis and return the standard deviation of each value as well.
        """

        self._free_energies = []

        # Now we want so solve A'*W*A dG = A'*W b - see the Wikipedia entry on weighted
        # least squares for the gory details
        A = self._compute_adjacency_matrix()
        W = self._compute_weight_matrix()
        b = self._compute_vector()
        if self._verbose:
            print("Adjacency is:  ")
            print(A)
            print("---------")
            print("Weight is: ")
            print(W)
            print("----------")
            print("b is: ")
            print(b)
            print('----------')
        AT_W_A = np.dot((np.dot(np.transpose(A), W)), A)
        AT_W_b = np.dot((np.dot(np.transpose(A), W)), b)
        dG = np.linalg.solve(AT_W_A, AT_W_b)

        dG = dG - np.mean(dG)

        h = self._error_estimate()

        err = np.zeros(shape=[len(dG), self.iterations], dtype="float64")
        for r in range(self.iterations):
            # Compute eb, which is b with random noise sized by the hysteresis
            # Given two values, the hysteresis is the difference between them,
            # and the standard deviation is sqrt(2)/2 times this
            eb = b + np.random.normal(0, h)

            # Solve with eb
            AT_W_eb = np.dot((np.dot(np.transpose(A), W)), eb)
            ex = np.linalg.solve(AT_W_A, AT_W_eb)

            ex = ex - np.mean(ex)

            # Append these results as a new column in err
            err[:, r] = ex

        # Compute standard deviations for each row of err:
        # use ddof=1 for sample estimate rather than population estimate
        std_dG = [np.std(err[i, :], ddof=1) for i in range(len(dG))]

        for c_idx in range(len(self._compoundList)):
            entry = {self._compoundList[c_idx]: dG[c_idx], 'error': std_dG[c_idx]}
            self._free_energies.append(entry)

    def _error_estimate(self, minh=0.4):
        """
            Create the hysteresis vector holding the pairwise hysteresis values
            The minimum hysteresis value is minh (also used if a link is unidirectional)
        """
        h = [0]
        for mol1 in self._ddG_edges:
            for mol2 in self._ddG_edges[mol1]:
                # Only handle links one way round: if both links present then only take one of them
                if (mol1 < mol2 or mol2 not in self._ddG_edges or mol1 not in self._ddG_edges[mol2]):
                    h.append((self._get_hysteresis(mol1, mol2, minh)))
        return h

    def _get_hysteresis(self, mol1, mol2, minh=0.4):
        """ Given two keys, return the hysteresis of the links, or minh
            if that is higher. If either link is missing returns minh.
        """
        try:
            eng1 = self._ddG_edges[mol1][mol2]
        except KeyError:
            return minh
        try:
            eng2 = self._ddG_edges[mol2][mol1]
        except KeyError:
            return minh

        # Compute relative hysteresis. This way high ddG edges result in the same hysteresis
        # penalty as low ddG edges, instead of high ddG edges being over-penalised.
        hys = eng1 + eng2

        rel_hys = hys/max(abs(eng1), abs(eng2))
        
        if self.balance_hysteresis:
            return max(minh, abs(rel_hys))
        else:
            return max(minh, abs(hys))

    def _get_avg_weight(self, mol1, mol2):
        """ Given two keys, return the average weight of
            the links, otherwise just whichever value exists
        """
        try:
            wt1 = self._weights[mol1][mol2]
        except KeyError:
            wt1 = None
        try:
            wt2 = self._weights[mol1][mol2]
        except KeyError:
            wt2 = None
        if wt1 and wt2:
            return (wt1 + wt2) / 2.0
        elif wt1:
            return wt1
        return wt2

    def write_free_energies(self, freeEnergies, filename=None, fmt=None):
        r"""Either write free energies to a file or std out
        Parameters
        ----------
        freeEnergies : list of dictionaries
            contains dictionaries with free energies and their errors
        filename : string
            file to which free energies should be written
            default = None
        fmt : string
            format string for the free energies, e.g. '%s, %f, %f\n'
            Default = None
        """
        if filename is not None:
            f = open(filename, 'w')
        else:
            print ('#FREE ENERGIES ARE:')
        for d in freeEnergies:
            for k, v in iter(d.items()):
                if k == 'error':
                    error = v
                else:
                    r_energy_k = k
                    r_energy_v = v
            if filename is not None:
                if fmt is None:
                    f.write('%s, %f, %f\n' % (r_energy_k, r_energy_v, error))
                else:
                    f.write(fmt % (r_energy_k, r_energy_v, error))
            else:
                if fmt is None:
                    print('{:10s} {:5.3f} +/- {:5.3f}'.format(r_energy_k, r_energy_v, error))
                else:
                    print (fmt % (r_energy_k, r_energy_v, error))
        if filename is not None:
            f.close()

    @property
    def weights(self):
        return self._weights

    @property
    def freeEnergyInKcal(self, balance_hysteresis=True):
        """ Return the free energies as a list of dictionaries
        """

        self._compute_free_energies()
        return self._free_energies

    @property
    def compoundList(self):
        return self._compoundList


class PerturbationGraph(object):
    """Populates a directed free energy perturbation graph"""

    def __init__(self):
        self._graph = None
        self._pathAverages = []
        self._weightedPathAverages = []
        self._weighted_paths = None
        self._compoundList = []
        self._free_energies = []
        warnings.warn(
            'PerturbationGraph is deprecated use the NetworkAnnalyser class instead.',
            DeprecationWarning, stacklevel=2)

    def populate_pert_graph(self, filename, delimiter=',', comments='#', nodetype=str,
                            data=(('weight', float), ('error', float))):
        r"""
        Reads data from a correctly formatted csv file into a networkx digraph
        Parameters
        ----------
        filename : String
            filename of the forward and backward perturbation generated from simulation output
            File structure should be:
            node1,node2,DG,eDG,other_attributes
        delimiter : String
            delimiter for network file 
            Default = ','
        comments : String
            Symbol used for comments in network file
            Default = '#'
        nodetype : String
            All nodes are usually identified by the compound name
        data : list
            Default, weight and error on Free energies of node
        """
        if self._graph is None:
            graph = nx.read_edgelist(filename, delimiter=delimiter, comments=comments, create_using=nx.DiGraph(),
                                     nodetype=nodetype, data=data)
            self._graph = self._symmetrize_graph(graph)
            self._compoundList = np.sort(self._graph.nodes())
        else:
            warnings.warn(UserWarning(
                "Warning...........Use the method add_data_to_graph, to add further data to an existing graph"))
            return 1

    def populate_graph(self, filename, delimiter=',', comments='#'):
        r"""alternative way of populating graph

        """
        g = nx.DiGraph()
        # add bit to check you can read file
        f = open(filename, 'r')
        data = f.readlines()
        f.close()

        # now populate the graph

        for d in data:
            l = d.strip().split(delimiter)
            if g.has_edge(l[0], l[1]):
                print('Must do something with edge: %s,%s' % (l[0], l[1]))
                a = g.get_edge_data(l[0], l[1])
                a['weight_list'].append(float(l[2]))
                g.edges[l[0], l[1]]['weight_list'] = a['weight_list']
            else:
                g.add_edge(l[0], l[1], weight=float(l[2]), error=float(l[3]), weight_list=[float(l[2])])
        for e in g.edges:
            e_data = g.get_edge_data(e[0], e[1])
            if 'weight_list' in e_data:
                w_list = e_data['weight_list']
                if len(w_list) > 1:
                    print(e[0], e[1])
                    print((w_list))
                    g.edges[e[0], e[1]]['weight'] = np.mean(w_list)
                    g.edges[e[0], e[1]]['error'] = np.std(w_list)
        self._graph = self._symmetrize_graph(g)
        self._compoundList = np.sort(self._graph.nodes())

    def add_data_to_graph(self, filename, delimiter=',', comments='#', nodetype=str,
                          data=(('weight', float), ('error', float))):
        r"""
        Adds data to an existing graph from a csv file in the right networkx format
        Parameters
        ----------
        filename : String
            filename of the forward and backward perturbation generated from simulation output
            File structure should be:
            node1,node2,DG,eDG,other_attributes
        delimiter : String
            delimiter for network file 
            Default = ','
        comments : String
            Symbol used for comments in network file
            Default = '#'
        nodetype : String
            All nodes are usually identified by the compound name
        data : list
            Default, weight and error on Free energies of node
        """
        newGraph = nx.read_edgelist(filename, delimiter=delimiter, comments=comments, create_using=nx.DiGraph(),
                                    nodetype=nodetype, data=data)
        newGraph = self._symmetrize_graph(newGraph)
        if self._graph != None:
            for u, v, w in newGraph.edges(data=True):
                if self._graph.has_edge(u, v):
                    z = self._graph.get_edge_data(u, v)
                    mean_edge = np.mean([z['weight'], w['weight']])
                    error = 0.5 * np.sqrt(z['error'] ** 2 + w['error'] ** 2)
                    self._graph.remove_edge(u, v)
                    self._graph.add_edge(u, v, weight=mean_edge, error=error)
                else:
                    self._graph.add_edge(u, v, w)
        else:
            self._graph = newGraph

    def remove_compound_from_graph(self, compound):
        r""" removes a node from the current graph
        Parameters
        ----------
        compound : string
            name of the compound to be removed from the graph

        """
        self._graph.remove_node(compound)
        self._compoundList = self._graph.nodes()

    def _symmetrize_graph(self, graph):
        r"""symmetrizes the graph and computes backward and forward averages where  given.
        Parameters
        ----------
        graph : networkx graph
            directed networkx graph

        Returns
        -------
        graph : networkx graph
            returns directed graph where, if not both a forward and backward edge are present a symmetrized reverse edge
             is included
        """
        symmetrizedGraph = nx.DiGraph()
        for u, v, w_forward in graph.edges(data=True):
            if graph.has_edge(v, u):
                w_backward = graph.get_edge_data(v, u)
                avg_weight_forw = np.mean([w_forward['weight'], -w_backward['weight']])
                avg_weight_back = -avg_weight_forw
                error = np.std([w_forward['weight'], -w_backward['weight']]) / np.sqrt(2.0)
                if error == 0.0:
                    error = np.mean([w_forward['error'], w_backward['error']])
                symmetrizedGraph.add_edge(u, v, weight=avg_weight_forw, error=error)
                symmetrizedGraph.add_edge(v, u, weight=avg_weight_back, error=error)
            else:
                symmetrizedGraph.add_edge(u, v, weight=w_forward['weight'], error=w_forward['error'])
        for u, v, w in symmetrizedGraph.edges(data=True):
            if not symmetrizedGraph.has_edge(v, u):
                assymetric_w = -w['weight']
                assymetric_e = -w['error']
                symmetrizedGraph.add_edge(v, u, weight=assymetric_w, error=assymetric_e)
        return symmetrizedGraph

    def format_free_energies(self, merge_BM=False, kT=0.594, intermed_ID=None, compound_order=None, weighted=True,
                             path_dictionary=None):
        r"""
         Parameters
        ----------
        fmt : string
            format string for the free energies, e.g. '%s, %f, %f\n'
            Default = None
        merge_BM : boolean
            true or false for binding modes using identified xxx_BMyyy, where xxx is the molecule name and yyy is the number of the binding mode
            Default = False
        kT : float
            simulation temperature times Boltzmann constant in [kcal/mol]
            Default = 0.594
        intermed_ID : string
            string identifier of intermediate simulated compounds, e.g 'INT'
            Default = None
        compound_order : list
            list of compounds
        weighted : boolean
            use weighted or none error weighted paths
        """
        if self._free_energies:
            self._free_energies = []
        mols = {}
        if weighted:
            if not self._weightedPathAverages and path_dictionary == None:
                print('compute weighted path averages for network first in order to format free energies')
                sys.exit(1)
            elif path_dictionary:
                freeEnergies = path_dictionary
            else:
                freeEnergies = self._weightedPathAverages
        else:
            if not self._pathAverages:
                print('compute path averages for network first in order to format free energies')
                sys.exit(1)
            else:
                freeEnergies = self._pathAveages

        for data in freeEnergies:
            keys = list(data.keys())
            if keys[0] != 'error':
                mol = keys[0]
            else:
                mol = keys[1]
            nrg = data[mol]
            err = data['error']
            if merge_BM:
                elems = mol.split("_BM")
                moln = elems[0]
            else:
                moln = mol
            try:
                mols[moln]
            except KeyError:
                mols[moln] = []
            mols[moln].append([nrg, err])
        ids = list(mols.keys())
        ids.sort()
        if compound_order is not None:
            if set(compound_order).issubset(ids):
                ids = compound_order
            else:
                print ("The list of compounds you provided does not match the ones stored in the perturbation network")
                print ("Compounds are:")
                print (ids)
                sys.exit(1)
        for mol in ids:
            if intermed_ID is not None:
                if mol.startswith(intermed_ID):
                    continue
            nrgtot = 0.0
            errtot = 0.0
            for nrg, err in mols[mol]:
                nrgtot += np.exp(-nrg / kT)
                errtot += err ** 2
            nrgtot = -kT * np.log(nrgtot)
            errtot = np.sqrt(errtot)
            a = {}
            a[mol] = nrgtot
            a['error'] = errtot
            self._free_energies.append(a)

    def write_free_energies(self, freeEnergies, filename=None, fmt=None):
        r"""Either write free energies to a file or std out
        Parameters
        ----------
        freeEnergies : list of dictionaries
            contains dictionaries with free energies and their errors
        filename : string
            file to which free energies should be written
            default = None
        fmt : string
            format string for the free energies, e.g. '%s, %f, %f\n'
            Default = None
        """
        if filename is not None:
            f = open(filename, 'w')
        else:
            print ('#FREE ENERGIES ARE:')
        for d in freeEnergies:
            for k, v in iter(d.items()):
                if k == 'error':
                    error = v
                else:
                    r_energy_k = k
                    r_energy_v = v
            if filename is not None:
                if fmt is None:
                    f.write('%s, %f, %f\n' % (r_energy_k, r_energy_v, error))
                else:
                    f.write(fmt % (r_energy_k, r_energy_v, error))
            else:
                if fmt is None:
                    print('{:10s} {:5.3f} +/- {:5.3f}'.format(r_energy_k, r_energy_v, error))
                else:
                    print (fmt % (r_energy_k, r_energy_v, error))
        if filename is not None:
            f.close()

    def shift_free_energies(self, shift_value=0.0):
        for d in self._free_energies:
            for k, v in iter(d.items()):
                if k != 'error':
                    d[k] = d[k] - shift_value

    def compute_average_paths(self, target_node):
        r"""
        Parameters
        ----------
        target_node : string
            node to which all possible paths are computed
        """
        # Get all relative free energies with respect to node x
        self._weighted_paths = False
        self._pathAverages = []
        for n in self._compoundList:
            paths = list(nx.all_simple_paths(self._graph, target_node, n, cutoff=12))
            err_list = []
            sum_list = []
            for p in paths:
                sum = 0
                error = 0.0
                for node in range(len(p) - 1):
                    sum = sum + self._graph.get_edge_data(p[node], p[node + 1])['weight']
                    error = error + self._graph.get_edge_data(p[node], p[node + 1])['error'] ** 2
                sum_list.append(sum)
                err_list.append(error)
                error = np.sqrt(error)
                # print r'DDG for path %s is %f +/- %f kcal/mol' %(p, sum, error)
            avg_sum = np.mean(np.array(sum_list))
            avg_err = np.mean(np.array(err_list))
            avg_std = np.std(np.array(sum_list))
            # print ("Average sum for path to %s is %f " %(n,avg_sum))
            a = {str(n): avg_sum}
            a['error'] = avg_std
            # a['error']=sqrt(avg_err)
            self._pathAverages.append(a)

    def compute_weighted_avg_paths(self, target_node):
        r""" computes all possible paths to a target node and returns a weighted average based on the errors along the edges of the path
        Parameters
        ----------
        target_node : string
            string name of the target node as defined in the networkx graph
        """
        # Get all relative free energies with respect to node x
        self._weighted_paths = True
        self._weightedPathAverages = []
        a = {target_node: 0.0}
        a['error'] = 0.0
        self._weightedPathAverages.append(a)
        for n in self._compoundList:
            if n == target_node:
                continue
            paths = list(nx.shortest_simple_paths(self._graph, target_node, n))
            err_list = []
            sum_list = []
            for p in paths:
                summing = 0
                error = 0.0
                for node in range(len(p) - 1):
                    summing = summing + self._graph.get_edge_data(p[node], p[node + 1])['weight']
                    error = error + self._graph.get_edge_data(p[node], p[node + 1])['error'] ** 2
                sum_list.append(summing)
                error = np.sqrt(error)
                err_list.append(error)
            err_list = np.array(err_list)
            sum_weights = np.sum(1.0 / err_list)
            path_weights = (1.0 / err_list) / sum_weights
            avg_sum = 0.0
            avg_err = 0.0
            for i in range(len(sum_list)):
                s = sum_list[i]
                avg_sum = avg_sum + (path_weights[i] * s)
                avg_err = avg_err + path_weights[i] * err_list[i] ** 2
            avg_err = np.sqrt(avg_err)
            a = {str(n): avg_sum}
            a['error'] = avg_err
            self._weightedPathAverages.append(a)

    def get_cycles(self, max_length=4, closure_threshold=1.0, print_all=False):
        r"""
        TODO: elaborate and find good way of saving this information 
        """
        # cycle closure
        cyc = nx.simple_cycles(self._graph)
        for c in cyc:
            sum = 0
            error = 0
            if len(c) > 2:
                sum = self._graph.get_edge_data(c[-1], c[0])['weight']
                error = (self._graph.get_edge_data(c[-1], c[0])['error']) ** 2
                for node in range(len(c) - 1):
                    sum = sum + self._graph.get_edge_data(c[node], c[node + 1])['weight']
                    error = error + (self._graph.get_edge_data(c[node], c[node + 1])['error']) ** 2
                error = np.sqrt(error)
                if len(c) <= max_length and not print_all:
                    if sum > closure_threshold:
                        print ('DDG for cycle %s is %.2f +/- %.2f kcal/mol' % (c, sum, error))
                if print_all:
                    print ('DDG for cycle %s is %.2f +/- %.2f kcal/mol' % (c, sum, error))

    def rename_compounds(self):
        warnings.warn(NotImplementedError)('This function is not implemented yet')
        sys.exit(1)

    @property
    def graph(self):
        return self._graph

    @property
    def pathAverages(self):
        r"""
        Return
        ------
        pathAverages : dictionary
            dictionary containing averaged free energies for each path, with paths weighted in the same way 
        """
        return self._pathAverages

    @property
    def weightedPathAverages(self):
        return self._weightedPathAverages

    @property
    def freeEnergyInKcal(self):
        if self._free_energies:
            return self._free_energies
        else:
            if self._weighted_paths:
                return self._weightedPathAverages
            else:
                return self._pathAverages

    @property
    def compoundList(self):
        return self._compoundList
