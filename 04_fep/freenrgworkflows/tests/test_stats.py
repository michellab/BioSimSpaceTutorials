import pytest
import numpy as np
import warnings
from networkanalysis.stats import *
from networkanalysis.jupyter import *


@pytest.fixture
def stats():
    return freeEnergyStats()


def test_calculate_r2(stats):
    a = np.array([2, 4, 5, 6, 4, 6, 8.9, 7])
    b = np.array([5, 4, 6, 7, 8, 7, 7.8, 6])
    r2 = 0.2931274858030396
    r = 0.5414124913622141
    test_r2, test_r = stats._calculate_r2(a, b)
    assert (pytest.approx(test_r2) == r2)
    assert (pytest.approx(test_r) == r)


def test_calculate_tau(stats):
    a = np.array([2, 4, 5, 6, 4, 6, 8.9, 7])
    b = np.array([5, 4, 6, 7, 8, 7, 7.8, 6])
    tau = 0.3461538461538462
    test_tau = stats._calculate_tau(a, b)
    assert (pytest.approx(test_tau) == tau)


def test_calculate_mue(stats):
    a = np.array([2, 4, 5, 6, 4, 6, 8.9, 7])
    b = np.array([5, 4, 6, 7, 8, 7, 7.8, 6])
    mue = 1.5125000000000002
    test_mue = stats._calculate_mue(a, b)
    assert (pytest.approx(test_mue) == mue)


def test_calculate_rmse(stats):
    a = np.array([2, 4, 5, 6, 4, 6, 8.9, 7])
    b = np.array([5, 4, 6, 7, 8, 7, 7.8, 6])
    rmse = 1.9432575742808775
    test_mue = stats._calculate_rmse(a, b)
    assert (pytest.approx(test_mue) == rmse)


@pytest.mark.parametrize('boundaries', [(-72), (4)])
def test_confidence_warnings(stats, boundaries):
    warn_string = 'Confidence interval needs to be between 0 and 1, please try something like 0.68 for one sigma confidence'
    with pytest.warns(UserWarning) as warnmessg:
        stats.confidence_interval = boundaries
        # warnings.warn(warn_string, UserWarning)
    assert len(warnmessg) == 1
    assert warnmessg[0].message.args[0] == warn_string


def test_confidence(stats):
    data = np.array([5, 4, 6, 7, 8, 7, 7.8, 6])
    assert (stats._confidence(data) == [6.0, 7.8])


@pytest.mark.parametrize('repeat', [(100), (20)])
def test_repeats(stats, repeat):
    exp_dat = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    comp = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    stats.generate_statistics(comp, exp_dat, repeats=repeat)
    assert (len(stats._R) == repeat)
    assert (len(stats._mue) == repeat)
    assert (len(stats._tau) == repeat)
    assert (len(stats._R2) == repeat)
    assert (len(stats._rmse) == repeat)


def test_array_conversion(stats):
    exp_dat = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    comp = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    stats.generate_statistics(comp, exp_dat, repeats=1)
    data_comp = np.array(stats.data_comp)
    data_exp = np.array(stats.data_exp)
    assert (np.all(data_comp[:, 0]) == np.all(data_exp))


def test_statistics(stats):
    exp_dat = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    comp = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    stats.generate_statistics(comp, exp_dat, repeats=1)
    comp_array = np.array(stats.data_comp)
    r = stats._calculate_r2(comp_array[:, 0], stats.data_exp)
    tau = stats._calculate_tau(comp_array[:, 0], stats.data_exp)
    mue = stats._calculate_mue(comp_array[:, 0], stats.data_exp)
    assert (pytest.approx(r) == (1.0, 1.0))
    assert (pytest.approx(tau) == 1.0)
    assert (pytest.approx(mue) == 0.0)


def test_properties(stats):
    exp_dat = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    comp = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    stats.generate_statistics(comp, exp_dat, repeats=100)
    assert (np.mean(stats._R) == stats.R_mean)
    assert (np.mean(stats._R2) == stats.R2_mean)
    assert (np.mean(stats._tau) == stats.tau_mean)
    assert (np.mean(stats._mue) == stats.mue_mean)
    R,p = scipy.stats.pearsonr(np.array(stats.data_comp)[:,0], np.array(stats.data_exp))
    assert (R == stats.R)


def test_property_errors(stats):
    stats._R = np.loadtxt('tests/io/R.dat')
    stats.confidence_interval = 0.68
    err1 = stats.R_confidence
    stats.confidence_interval = 0.95
    err2 = stats.R_confidence
    assert (err1 != err2)


@pytest.mark.parametrize('boundaries', [(0.95)])
def test_property_errors(stats, boundaries):
    stats._R = np.loadtxt('tests/io/R.dat')
    stats._tau = np.loadtxt('tests/io/tau.dat')
    stats._mue = np.loadtxt('tests/io/mue.dat')
    stats._R2 = stats._R ** 2
    stats.confidence_interval = boundaries
    np.testing.assert_allclose(stats.R_confidence, np.array([0.89136482, 0.8682079, 0.9338724]), rtol=1e-8, atol=0)
    np.testing.assert_allclose(stats.tau_confidence, np.array([1.0, 1.0, 1.0]), rtol=1e-8, atol=0)
    np.testing.assert_allclose(stats.mue_confidence, np.array([0.43437384, 0.35623369, 0.53095783]), rtol=1e-6, atol=0)
    np.testing.assert_allclose(stats.R2_confidence, np.array([0.79453165, 0.75378495, 0.87211766]), rtol=1e-8, atol=0)
