import pytest
import warnings
import networkx as nx
import os
from networkanalysis.experiments import *


@pytest.fixture
def ExpData300():
    return ExperimentalData()


@pytest.mark.parametrize('temperature', [(300.0), (400.0)])
def test_temperature(ExpData300, temperature):
    expData = ExperimentalData(temperature)
    assert (expData._kTkcal == 0.0019872041 * temperature)
    assert (expData._kTkJ == 0.0083144621 * temperature)
    if temperature == 300.0:
        assert (ExpData300._kTkcal == expData._kTkcal)


def test_read_comp_data_kj(ExpData300):
    filename = os.path.join(os.getcwd(), 'tests', 'io', 'comp.dat')
    err_string = 'This has not been implemented yet'
    with pytest.raises(NotImplementedError) as errmessg:
        ExpData300.read_free_energies(filename, kcal=False)


def test_read_comp_data_kj(ExpData300):
    filename = os.path.join(os.getcwd(), 'tests', 'io', 'comp.dat')
    ExpData300.read_free_energies(filename, kcal=True)
    assert (ExpData300.freeEnergiesInKcal[-1]['FXR100'] == 7.01744151024)
    assert (ExpData300.freeEnergiesInKcal[-1]['error'] == 0.206262061369)


def test_keys(ExpData300):
    filename = os.path.join(os.getcwd(), 'tests', 'io', 'comp.dat')
    ExpData300.read_free_energies(filename, kcal=True)
    keys = []
    f = open(filename, 'r')
    for line in f.readlines():
        if line.startswith('#'):
            continue
        fields = line.split(',')
        if fields[1] == 'NoPred':
            continue
        keys.append(fields[0])
    f.close()
    assert (keys == ExpData300.compoundList)


def test_from_IC50s_keys(ExpData300):
    filename = os.path.join(os.getcwd(), 'tests', 'io', 'ic50_exp.dat')
    ExpData300.compute_DDG_from_IC50s(filename, 'FXR17')
    cList = ['FXR45', 'FXR93', 'FXR91', 'FXR96', 'FXR46', 'FXR95', 'FXR98', 'FXR49', 'FXR17', 'FXR48', 'FXR47', 'FXR99',
             'FXR102', 'FXR101', 'FXR100']
    assert (ExpData300.compoundList == cList)


def test_from_kD_keys(ExpData300):
    filename = os.path.join(os.getcwd(), 'tests', 'io', 'kD_exp.dat')
    ExpData300.compute_DDG_from_IC50s(filename, 'FXR17')
    cList = ['FXR45', 'FXR93', 'FXR91', 'FXR96', 'FXR46', 'FXR95', 'FXR98', 'FXR49', 'FXR17', 'FXR48', 'FXR47', 'FXR99',
             'FXR102', 'FXR101', 'FXR100']
    assert (ExpData300.compoundList == cList)
