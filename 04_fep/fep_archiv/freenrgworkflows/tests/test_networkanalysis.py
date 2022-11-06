import pytest
import warnings
import networkx as nx
from networkanalysis.networkanalysis import *


@pytest.fixture
def nA():
    return NetworkAnalyser()


def test_compoundList(nA):
    nA.read_perturbations_pandas('tests/io/graph.csv')
    assert ('FXR17' in nA.compoundList)


def test_dG_simple(nA):
    nA.read_perturbations_pandas('tests/io/simple.csv',comments='#')
    energies = nA.freeEnergyInKcal
    x, y = convert_energy_list(energies)
    validate_results(x, [-0.5, 0.5])
    validate_results(y, [0.4, 0.4], 0.05)


def test_perfectcycle(nA):
    nA.read_perturbations_pandas('tests/io/perfectcycle.csv',comments='#')
    energies = nA.freeEnergyInKcal
    x, y = convert_energy_list(energies)
    validate_results(x, [-0.5, 0.5, 0])
    validate_results(y, [0.4, 0.4, 0.4], 0.05)


def test_inconsistentcycle(nA):
    nA.read_perturbations_pandas('tests/io/inconsistentcycle.csv',comments='#')
    energies = nA.freeEnergyInKcal
    x, y = convert_energy_list(energies)
    validate_results(x, [0, 0, 0])
    validate_results(y, [0.4, 0.4, 0.4], 0.05)


def test_inconsistentcycle_weights(nA):
    nA.read_perturbations_pandas("tests/io/inconsistentcycle_weights.csv",comments='#')
    energies = nA.freeEnergyInKcal
    x, y = convert_energy_list(energies)
    validate_results(x, [-0.5, 0.5, 0])
    validate_results(y, [0.44, 0.44, 0.38], 0.05)


def test_large_hysteresis(nA):
    nA.read_perturbations_pandas("tests/io/large_hysteresis.csv",comments='#')
    energies = nA.freeEnergyInKcal
    x, y = convert_energy_list(energies)
    validate_results(x, [-0.5, 0.5, 0])
    validate_results(y, [0.66, 0.66, 0.66], 0.05)


def test_vlarge_hysteresis(nA):
    nA.read_perturbations_pandas("tests/io/vlarge_hysteresis.csv",comments='#')
    energies = nA.freeEnergyInKcal
    x, y = convert_energy_list(energies)
    validate_results(x, [-0.5, 0.5, 0])
    validate_results(y, [1.60, 1.60, 1.60], 0.05)


def test_large_cycle(nA):
    nA.read_perturbations_pandas("tests/io/large_cycle.csv",comments='#')
    energies = nA.freeEnergyInKcal
    x, y = convert_energy_list(energies)
    x = [t - x[0] for t in x]  # Set the first mol to zero
    validate_results(x, [0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0])
    validate_results(y, [0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65], 0.05)


def test_large_cycle_poor_link(nA):
    nA.read_perturbations_pandas("tests/io/large_cycle_poor_link.csv",comments='#')
    energies = nA.freeEnergyInKcal
    x, y = convert_energy_list(energies)
    x = [t - x[0] for t in x]  # Set the first mol to zero
    validate_results(x, [0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0], 0.01)
    # These are the error values when the "poor" link is deleted: we should get the same when it
    # is present but has a very low weight
    validate_results(y, [1.1772, 0.9524, 0.7728, 0.6652, 0.6612, 0.7720, 0.9576, 1.1777], 0.05)


def test_DG_noise_A(nA):
    nA.read_perturbations_pandas("tests/io/noise0.csv",comments='#')
    energies = nA.freeEnergyInKcal
    x, y = convert_energy_list(energies)
    x = [t - x[0] for t in x]  # Set the first mol to zero
    validate_results(x, [0, 1, 2, 3, 4, 4])
    # Errors: given the network, m4 should have the lowest error and m6 the highest
    validate_results(y, [0.54, 0.53, 0.43, 0.35, 0.53, 0.73], 0.05)


def test_DG_noise_B(nA):
    # Same network, random error with std dev 0.5 added to all deltag measurements
    nA.read_perturbations_pandas("tests/io/noise0.5.csv", comments='#')
    energies = nA.freeEnergyInKcal
    x, y = convert_energy_list(energies)
    print (x)
    x = [t - x[0] for t in x]  # Set the first mol to zero
    validate_results(x, [0.0, 0.935, 2.425, 3.365, 4.64, 4.57])
    # Errors: given the network, m4 should have the lowest error and m6 the highest
    validate_results(y, [0.845, 0.902, 0.591, 0.564, 0.779, 0.864], 0.05)


def test_value_error(nA):
    with pytest.raises(ValueError):
        nA.read_perturbations_pandas("tests/io/too_few_columns.csv")


def validate_results(x, xcorrect, delta=0.01):
    print("Checking", x, "against", xcorrect, "delta", delta)
    assert len(x) == len(xcorrect)
    for i in range(len(x)):
        assert x[i] == pytest.approx(xcorrect[i], abs=delta)


def convert_energy_list(energies):
    x = []
    y = []
    for e in energies:
        keys = list(e.keys())
        idx = keys.index('error')
        if idx == 1:
            x.append(e[keys[0]])
        else:
            x.append(e[keys[1]])
        y.append(e['error'])
    return x, y

def test_add_replicates(nA):
    nA.read_perturbations_pandas('tests/io/replicate_1.csv',comments='#')
    energies1 = nA.freeEnergyInKcal
    x1, y1 = convert_energy_list(energies1)

    nA.add_data_to_graph_pandas('tests/io/replicate_2.csv',comments='#')
    energies2 = nA.freeEnergyInKcal
    x2, y2 = convert_energy_list(energies1)

    assert (x1 != x2) and (y1 != y2)


@pytest.fixture
def pG():
    return PerturbationGraph()


def test_populate_pert_graph_not_None(pG):
    pG.populate_pert_graph('tests/io/graph.csv')
    assert (pG.graph != None)


def test_compoundList(pG):
    pG.populate_pert_graph('tests/io/graph.csv')
    assert ('FXR17' in pG.compoundList)


def test_double_call_pert_graph_not_None(pG):
    pG.populate_pert_graph('tests/io/graph.csv')
    warn_string = 'Warning...........Use the method add_data_to_graph, to add further data to an existing graph'
    with pytest.warns(UserWarning) as warnmessg:
        pG.populate_pert_graph('tests/io/graph.csv')
    assert len(warnmessg) == 1
    assert warnmessg[0].message.args[0] == warn_string


def test_perfect_graph(pG):
    pG.populate_pert_graph('tests/io/test_perfect_graph.csv')
    G = pG.graph
    for edge in G.edges:
        data = G.get_edge_data(edge[0], edge[1])
        assert data['error'] is not 0.0
# def test_symmetrize_graph():
#    newGraph = nx.read_edgelist('tests/io/graph.csv', delimiter=',', comments='#', nodetype=str, data=(('weight', float),('error',float)))

# def test_cycle_closure(pG, capsys):


# def test_compute_weighted_avg_paths():

# def test_compute_average_paths():

# def test_remove_compound_from_graph():

# def test_format_free_energies():
