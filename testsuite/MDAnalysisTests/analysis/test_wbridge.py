from __future__ import print_function, absolute_import
from six import StringIO
from collections import defaultdict

import numpy as np
from numpy.testing import (
    assert_equal, assert_array_equal, assert_almost_equal,
    assert_array_almost_equal, assert_allclose,)
import pytest

import MDAnalysis
import MDAnalysis.analysis.hbonds
from MDAnalysis.analysis.hbonds.wbridge_analysis import WaterBridgeAnalysis
from MDAnalysisTests.datafiles import PDB_helix, GRO, XTC, waterPSF, waterDCD

# For type guessing:
from MDAnalysis.topology.core import guess_atom_type
from MDAnalysis.core.topologyattrs import Atomtypes

def test_import_from_hbonds():
    try:
        from MDAnalysis.analysis.hbonds import WaterBridgeAnalysis
    except ImportError:
        raise AssertionError("Issue #2064 not fixed: "
                             "importing WaterBridgeAnalysis from "
                             "MDAnalysis.analysis.hbonds failed.'")

class TestWaterBridgeAnalysis(object):
    def test_donor_accepter(self):
        '''Test zeroth order donor to acceptor hydrogen bonding'''
        grofile = '''Test gro file
3
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    4ALA      O    3   0.300   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)', order=0)
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal(list(network.keys())[0][:4], (1, 2, ('ALA', 1, 'H'), ('ALA', 4, 'O')))

    def test_donor_accepter_pbc(self):
            '''Test zeroth order donor to acceptor hydrogen bonding'''
            grofile = '''Test gro file
3
    1ALA      N    1   0.800   0.000   0.000
    1ALA      H    2   0.900   0.000   0.000
    4ALA      O    3   0.100   0.000   0.000
  1.0   1.0   1.0'''
            u = MDAnalysis.Universe(StringIO(grofile), format='gro')
            wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)', order=0, pbc=True)
            wb.run(verbose=False)
            network = wb._network[0]
            assert_equal(list(network.keys())[0][:4], (1, 2, ('ALA', 1, 'H'), ('ALA', 4, 'O')))

    def test_accepter_donor(self):
        '''Test zeroth order acceptor to donor hydrogen bonding'''
        grofile = '''Test gro file
3
    1ALA      O    1   0.000   0.000   0.000
    4ALA      H    2   0.200   0.000   0.000
    4ALA      N    3   0.300   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)', order=0)
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal(list(network.keys())[0][:4], (0, 1, ('ALA', 1, 'O'), ('ALA', 4, 'H')))

    def test_acceptor_water_accepter(self):
        '''Test case where the hydrogen bond acceptor from selection 1 form
        water bridge with hydrogen bond acceptor from selection 2'''
        grofile = '''Test gro file
5
    1ALA      O    1   0.000   0.000   0.000
    2SOL     OW    2   0.300   0.000   0.000
    2SOL    HW1    3   0.200   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)')
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal(list(network.keys())[0][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O')))
        assert_equal(second[list(second.keys())[0]], None)

    def test_donor_water_accepter(self):
        '''Test case where the hydrogen bond donor from selection 1 form
        water bridge with hydrogen bond acceptor from selection 2'''
        grofile = '''Test gro file
5
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    2SOL     OW    3   0.300   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)')
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal(list(network.keys())[0][:4], (1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW')))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O')))
        assert_equal(second[list(second.keys())[0]], None)

    def test_acceptor_water_donor(self):
        '''Test case where the hydrogen bond acceptor from selection 1 form
        water bridge with hydrogen bond donor from selection 2'''
        grofile = '''Test gro file
5
    1ALA      O    1   0.000   0.000   0.000
    2SOL     OW    2   0.300   0.000   0.000
    2SOL    HW1    3   0.200   0.000   0.000
    4ALA      H    4   0.500   0.000   0.000
    4ALA      N    5   0.600   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)')
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal(list(network.keys())[0][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (1, 3, ('SOL', 2, 'OW'), ('ALA', 4, 'H')))
        assert_equal(second[list(second.keys())[0]], None)

    def test_donor_water_donor(self):
        '''Test case where the hydrogen bond donor from selection 1 form
        water bridge with hydrogen bond donor from selection 2'''
        grofile = '''Test gro file
5
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    2SOL     OW    3   0.300   0.000   0.000
    4ALA      H    4   0.500   0.000   0.000
    4ALA      N    5   0.600   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)')
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal(list(network.keys())[0][:4], (1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW')))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (2, 3, ('SOL', 2, 'OW'), ('ALA', 4, 'H')))
        assert_equal(second[list(second.keys())[0]], None)

    def test_empty(self):
        '''Test case where no water bridge exists'''
        grofile = '''Test gro file
5
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    2SOL     OW    3   3.000   0.000   0.000
    4ALA      H    4   0.500   0.000   0.000
    4ALA      N    5   0.600   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein', 'protein')
        wb.run(verbose=False)
        assert_equal(wb._network[0], defaultdict(dict))

    def test_same_selection(self):
        '''
        This test tests that if the selection 1 and selection 2 are both protein.
        However, the protein only forms one hydrogen bond with the water.
        This entry won't be included.
        :return:
        '''
        grofile = '''Test gro file
3
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    2SOL     OW    3   0.300   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein', 'protein')
        wb.run(verbose=False)
        assert_equal(wb._network[0], defaultdict(dict))

    def test_acceptor_2water_accepter(self):
        '''Test case where the hydrogen bond acceptor from selection 1 form second order
        water bridge with hydrogen bond acceptor from selection 2'''
        grofile = '''Test gro file
7
    1ALA      O    1   0.000   0.000   0.000
    2SOL     OW    2   0.300   0.000   0.000
    2SOL    HW1    3   0.200   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    3SOL     OW    5   0.600   0.000   0.000
    3SOL    HW1    6   0.700   0.000   0.000
    4ALA      O    7   0.900   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        # test first order
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)')
        wb.run(verbose=False)
        assert_equal(wb._network[0], defaultdict(dict))
        # test second order
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)', order=2)
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal(list(network.keys())[0][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW')))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 6, ('SOL', 3, 'HW1'), ('ALA', 4, 'O')))
        assert_equal(third[list(third.keys())[0]], None)
        # test third order
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)', order=3)
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal(list(network.keys())[0][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW')))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 6, ('SOL', 3, 'HW1'), ('ALA', 4, 'O')))
        assert_equal(third[list(third.keys())[0]], None)

    def test_acceptor_3water_accepter(self):
        '''Test case where the hydrogen bond acceptor from selection 1 form third order
        water bridge with hydrogen bond acceptor from selection 2'''
        grofile = '''Test gro file
9
    1ALA      O    1   0.000   0.000   0.000
    2SOL     OW    2   0.300   0.000   0.000
    2SOL    HW1    3   0.200   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    3SOL     OW    5   0.600   0.000   0.000
    3SOL    HW1    6   0.700   0.000   0.000
    4SOL     OW    7   0.900   0.000   0.000
    4SOL    HW1    8   1.000   0.000   0.000
    5ALA      O    9   1.200   0.000   0.000
  10.0   10.0   10.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 5)', order=2)
        wb.run(verbose=False)
        assert_equal(wb._network[0], defaultdict(dict))

        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 5)', order=3)
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal(list(network.keys())[0][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW')))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 6, ('SOL', 3, 'HW1'), ('SOL', 4, 'OW')))
        fourth = third[list(third.keys())[0]]
        assert_equal(list(fourth.keys())[0][:4], (7, 8, ('SOL', 4, 'HW1'), ('ALA', 5, 'O')))
        assert_equal(fourth[list(fourth.keys())[0]], None)

        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 5)', order=4)
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal(list(network.keys())[0][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW')))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 6, ('SOL', 3, 'HW1'), ('SOL', 4, 'OW')))
        fourth = third[list(third.keys())[0]]
        assert_equal(list(fourth.keys())[0][:4], (7, 8, ('SOL', 4, 'HW1'), ('ALA', 5, 'O')))
        assert_equal(fourth[list(fourth.keys())[0]], None)

    def test_acceptor_4water_accepter(self):
        '''Test case where the hydrogen bond acceptor from selection 1 form fourth order
        water bridge with hydrogen bond acceptor from selection 2'''
        grofile = '''Test gro file
11
    1ALA      O    1   0.000   0.000   0.000
    2SOL     OW    2   0.300   0.000   0.000
    2SOL    HW1    3   0.200   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    3SOL     OW    5   0.600   0.000   0.000
    3SOL    HW1    6   0.700   0.000   0.000
    4SOL     OW    7   0.900   0.000   0.000
    4SOL    HW1    8   1.000   0.000   0.000
    5SOL     OW    9   1.200   0.000   0.000
    5SOL    HW1   10   1.300   0.000   0.000
    6ALA      O   11   1.400   0.000   0.000
  10.0   10.0   10.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 6)', order=3)
        wb.run(verbose=False)
        assert_equal(wb._network[0], defaultdict(dict))

        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 6)', order=4)
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal(list(network.keys())[0][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW')))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 6, ('SOL', 3, 'HW1'), ('SOL', 4, 'OW')))
        fourth = third[list(third.keys())[0]]
        assert_equal(list(fourth.keys())[0][:4], (7, 8, ('SOL', 4, 'HW1'), ('SOL', 5, 'OW')))
        fifth = fourth[list(fourth.keys())[0]]
        assert_equal(list(fifth.keys())[0][:4], (9, 10, ('SOL', 5, 'HW1'), ('ALA', 6, 'O')))
        assert_equal(fifth[list(fifth.keys())[0]], None)

        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 6)', order=5)
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal(list(network.keys())[0][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW')))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 6, ('SOL', 3, 'HW1'), ('SOL', 4, 'OW')))
        fourth = third[list(third.keys())[0]]
        assert_equal(list(fourth.keys())[0][:4], (7, 8, ('SOL', 4, 'HW1'), ('SOL', 5, 'OW')))
        fifth = fourth[list(fourth.keys())[0]]
        assert_equal(list(fifth.keys())[0][:4], (9, 10, ('SOL', 5, 'HW1'), ('ALA', 6, 'O')))
        assert_equal(fifth[list(fifth.keys())[0]], None)

    def test_acceptor_22water_accepter(self):
        '''Test case where the hydrogen bond acceptor from selection 1 form a second order
        water bridge with hydrogen bond acceptor from selection 2
        and the last water is linked to two residues in selection 2'''
        grofile = '''Test gro file
9
    1ALA      O    1   0.000   0.000   0.000
    2SOL     OW    2   0.300   0.000   0.000
    2SOL    HW1    3   0.200   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    3SOL     OW    5   0.600   0.000   0.000
    3SOL    HW1    6   0.700   0.000   0.000
    3SOL    HW2    7   0.600   0.100   0.000
    4ALA      O    8   0.900   0.000   0.000
    5ALA      O    9   0.600   0.300   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4 or resid 5)', order=2)
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal(list(network.keys())[0][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW')))
        third = second[list(second.keys())[0]]
        assert_equal([(5, 7, ('SOL', 3, 'HW1'), ('ALA', 4, 'O')), (6, 8, ('SOL', 3, 'HW2'), ('ALA', 5, 'O'))],
                     sorted([key[:4] for key in list(third.keys())]))

    def test_timeseries(self):
        '''Test if the time series data is correctly generated'''
        grofile = '''Test gro file
9
    1ALA      O    1   0.000   0.000   0.000
    2SOL     OW    2   0.300   0.000   0.000
    2SOL    HW1    3   0.200   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    3SOL     OW    5   0.600   0.000   0.000
    3SOL    HW1    6   0.700   0.000   0.000
    3SOL    HW2    7   0.600   0.100   0.000
    4ALA      O    8   0.900   0.000   0.000
    5ALA      O    9   0.600   0.300   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4 or resid 5)', order=2)
        wb.run(verbose=False)
        wb.timeseries
        timeseries = sorted(wb._timeseries[0])
        assert_equal(timeseries[0][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        assert_equal(timeseries[1][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        assert_equal(timeseries[2][:4], (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW')))
        assert_equal(timeseries[3][:4], (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW')))
        assert_equal(timeseries[4][:4], (5, 7, ('SOL', 3, 'HW1'), ('ALA', 4, 'O')))
        assert_equal(timeseries[5][:4], (6, 8, ('SOL', 3, 'HW2'), ('ALA', 5, 'O')))

    def test_acceptor_12water_accepter(self):
        '''Test of independent first order and second can be recognised correctely'''
        grofile = '''Test gro file
12
    1ALA      O    1   0.000   0.000   0.000
    2SOL     OW    2   0.300   0.000   0.000
    2SOL    HW1    3   0.200   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
    5ALA      O    6   0.000   1.000   0.000
    6SOL     OW    7   0.300   1.000   0.000
    6SOL    HW1    8   0.200   1.000   0.000
    6SOL    HW2    9   0.400   1.000   0.000
    7SOL     OW   10   0.600   1.000   0.000
    7SOL    HW1   11   0.700   1.000   0.000
    8ALA      O   12   0.900   1.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1 or resid 5)', 'protein and (resid 4 or resid 8)', order=1)
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal(list(network.keys())[0][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O')))
        assert_equal(second[list(second.keys())[0]], None)
        network = wb._network[0]
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1 or resid 5)', 'protein and (resid 4 or resid 8)', order=2)
        wb.run(verbose=False)
        network = wb._network[0]
        assert_equal([(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')), (5, 7, ('ALA', 5, 'O'), ('SOL', 6, 'HW1'))],
                     sorted([key[:4] for key in list(network.keys())]))

    def test_count_by_type_single_link(self):
        '''
        This test tests the simplest water bridge to see if count_by_type() works.
        :return:
        '''
        grofile = '''Test gro file
5
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    2SOL     OW    3   0.300   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
 1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)')
        wb.run(verbose=False)
        assert_equal(wb.count_by_type(), [(1, 4, 'ALA', 1, 'H', 'ALA', 4, 'O',  1.)])

    def test_count_by_type_multiple_link(self):
        '''
        This test tests if count_by_type() can give the correct result for more than 1 links.
        :return:
        '''
        grofile = '''Test gro file
12
    1ALA      O    1   0.000   0.000   0.000
    2SOL     OW    2   0.300   0.000   0.000
    2SOL    HW1    3   0.200   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
    5ALA      O    6   0.000   1.000   0.000
    6SOL     OW    7   0.300   1.000   0.000
    6SOL    HW1    8   0.200   1.000   0.000
    6SOL    HW2    9   0.400   1.000   0.000
    7SOL     OW   10   0.600   1.000   0.000
    7SOL    HW1   11   0.700   1.000   0.000
    8ALA      O   12   0.900   1.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1 or resid 5)', 'protein and (resid 4 or resid 8)', order=2)
        wb.run(verbose=False)
        assert_equal(sorted(wb.count_by_type()), [[0, 4, 'ALA', 1, 'O', 'ALA', 4, 'O', 1.0], [5, 11, 'ALA', 5, 'O', 'ALA', 8, 'O', 1.0]])


    def test_count_by_type_multiple_frame(self):
        '''
        This test tests if count_by_type() works in multiply situations.
        :return:
        '''
        grofile = '''Test gro file
5
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    2SOL     OW    3   0.300   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
 1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)', order=4)
        # Build an dummy WaterBridgeAnalysis object for testing
        wb._network = []
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (1, 3, ('SOL', 2, 'OW'), ('ALA', 4, 'H'), 2.0, 180.0): None}})
        wb._network.append({(1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 180.0): {
            (2, 3, ('SOL', 2, 'OW'), ('ALA', 4, 'H'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('ALA', 4, 'O'), 2.0, 180.0): None}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('SOL', 4, 'OW'), 2.0, 180.0): {
                    (7, 8, ('SOL', 4, 'HW1'), ('ALA', 5, 'O'), 2.0, 180.0): None}}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('SOL', 4, 'OW'), 2.0, 180.0): {
                    (7, 8, ('SOL', 4, 'HW1'), ('SOL', 5, 'OW'), 2.0, 180.0): {
                        (9, 10, ('SOL', 5, 'HW1'), ('ALA', 6, 'O'), 1.0, 180.0): None}}}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 7, ('SOL', 3, 'HW1'), ('ALA', 4, 'O'), 2.0, 180.0): None,
                (6, 8, ('SOL', 3, 'HW2'), ('ALA', 5, 'O'), 2.0, 180.0): None}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None},
            (5, 7, ('ALA', 5, 'O'), ('SOL', 6, 'HW1'), 2.0, 180.0): {
                (8, 9, ('SOL', 6, 'HW2'), ('SOL', 7, 'OW'), 2.0, 180.0): {
                    (10, 11, ('SOL', 7, 'HW1'), ('ALA', 8, 'O'), 2.0, 180.0): None}}})
        result = [(0, 3, 'ALA', 1, 'O', 'ALA', 4, 'H', 0.1),
                  (0, 4, 'ALA', 1, 'O', 'ALA', 4, 'O', 0.3),
                  (0, 6, 'ALA', 1, 'O', 'ALA', 4, 'O', 0.1),
                  (0, 7, 'ALA', 1, 'O', 'ALA', 4, 'O', 0.1),
                  (0, 8, 'ALA', 1, 'O', 'ALA', 5, 'O', 0.2),
                  (0, 10, 'ALA', 1, 'O', 'ALA', 6, 'O', 0.1),
                  (1, 3, 'ALA', 1, 'H', 'ALA', 4, 'H', 0.1),
                  (1, 4, 'ALA', 1, 'H', 'ALA', 4, 'O', 0.1),
                  (5, 11, 'ALA', 5, 'O', 'ALA', 8, 'O', 0.1)]
        assert_equal(sorted(wb.count_by_type()), result)

    def test_count_by_type_filter(self):
        '''
        This test tests if modifying analysis_func
        allows some results to be filtered out in count_by_type().
        :return:
        '''
        grofile = '''Test gro file
5
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    2SOL     OW    3   0.300   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
 1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)', order=4)
        # Build an dummy WaterBridgeAnalysis object for testing
        wb._network = []
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (1, 3, ('SOL', 2, 'OW'), ('ALA', 4, 'H'), 2.0, 180.0): None}})
        wb._network.append({(1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 180.0): {
            (2, 3, ('SOL', 2, 'OW'), ('ALA', 4, 'H'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('ALA', 4, 'O'), 2.0, 180.0): None}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('SOL', 4, 'OW'), 2.0, 180.0): {
                    (7, 8, ('SOL', 4, 'HW1'), ('ALA', 5, 'O'), 2.0, 180.0): None}}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('SOL', 4, 'OW'), 2.0, 180.0): {
                    (7, 8, ('SOL', 4, 'HW1'), ('SOL', 5, 'OW'), 2.0, 180.0): {
                        (9, 10, ('SOL', 5, 'HW1'), ('ALA', 6, 'O'), 1.0, 180.0): None}}}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 7, ('SOL', 3, 'HW1'), ('ALA', 4, 'O'), 2.0, 180.0): None,
                (6, 8, ('SOL', 3, 'HW2'), ('ALA', 5, 'O'), 2.0, 180.0): None}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None},
            (5, 7, ('ALA', 5, 'O'), ('SOL', 6, 'HW1'), 2.0, 180.0): {
                (8, 9, ('SOL', 6, 'HW2'), ('SOL', 7, 'OW'), 2.0, 180.0): {
                    (10, 11, ('SOL', 7, 'HW1'), ('ALA', 8, 'O'), 2.0, 180.0): None}}})

        def analysis(current, output):
            s1_index, to_index, (s1_resname, s1_resid, s1_name), (to_resname, to_resid, to_name), dist, angle = \
                current[0]
            from_index, s2_index, (from_resname, from_resid, from_name), (s2_resname, s2_resid, s2_name), dist, angle = \
                current[-1]
            key = (s1_index, s2_index, s1_resname, s1_resid, s1_name, s2_resname, s2_resid, s2_name)
            if s2_name == 'H':
                output[key] += 1
        result = [((0, 3, 'ALA', 1, 'O', 'ALA', 4, 'H'), 0.1),
                  ((1, 3, 'ALA', 1, 'H', 'ALA', 4, 'H'), 0.1)]
        assert_equal(sorted(wb.count_by_type(analysis_func=analysis)), result)

    def test_count_by_type_merge(self):
        '''
        This test tests if modifying analysis_func
        allows some same residue to be merged in count_by_type().
        :return:
        '''
        grofile = '''Test gro file
5
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    2SOL     OW    3   0.300   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
 1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)', order=4)
        # Build an dummy WaterBridgeAnalysis object for testing
        wb._network = []
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (1, 3, ('SOL', 2, 'OW'), ('ALA', 4, 'H'), 2.0, 180.0): None}})
        wb._network.append({(1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 180.0): {
            (2, 3, ('SOL', 2, 'OW'), ('ALA', 4, 'H'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('ALA', 4, 'O'), 2.0, 180.0): None}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('SOL', 4, 'OW'), 2.0, 180.0): {
                    (7, 8, ('SOL', 4, 'HW1'), ('ALA', 5, 'O'), 2.0, 180.0): None}}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('SOL', 4, 'OW'), 2.0, 180.0): {
                    (7, 8, ('SOL', 4, 'HW1'), ('SOL', 5, 'OW'), 2.0, 180.0): {
                        (9, 10, ('SOL', 5, 'HW1'), ('ALA', 6, 'O'), 1.0, 180.0): None}}}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 7, ('SOL', 3, 'HW1'), ('ALA', 4, 'O'), 2.0, 180.0): None,
                (6, 8, ('SOL', 3, 'HW2'), ('ALA', 5, 'O'), 2.0, 180.0): None}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None},
            (5, 7, ('ALA', 5, 'O'), ('SOL', 6, 'HW1'), 2.0, 180.0): {
                (8, 9, ('SOL', 6, 'HW2'), ('SOL', 7, 'OW'), 2.0, 180.0): {
                    (10, 11, ('SOL', 7, 'HW1'), ('ALA', 8, 'O'), 2.0, 180.0): None}}})
        def analysis(current, output):
            s1_index, to_index, (s1_resname, s1_resid, s1_name), (to_resname, to_resid, to_name), dist, angle = \
                current[0]
            from_index, s2_index, (from_resname, from_resid, from_name), (s2_resname, s2_resid, s2_name), dist, angle = \
                current[-1]
            key = (s1_resname, s1_resid, s2_resname, s2_resid)
            output[key] += 1
        result = [(('ALA', 1, 'ALA', 4), 0.8),
                  (('ALA', 1, 'ALA', 5), 0.2),
                  (('ALA', 1, 'ALA', 6), 0.1),
                  (('ALA', 5, 'ALA', 8), 0.1)]
        assert_equal(sorted(wb.count_by_type(analysis_func=analysis)), result)

    def test_count_by_type_order(self):
        '''
        This test tests if modifying analysis_func
        allows the order of water bridge to be separated in count_by_type().
        :return:
        '''
        grofile = '''Test gro file
5
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    2SOL     OW    3   0.300   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
 1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)', order=4)
        # Build an dummy WaterBridgeAnalysis object for testing
        wb._network = []
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (1, 3, ('SOL', 2, 'OW'), ('ALA', 4, 'H'), 2.0, 180.0): None}})
        wb._network.append({(1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 180.0): {
            (2, 3, ('SOL', 2, 'OW'), ('ALA', 4, 'H'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('ALA', 4, 'O'), 2.0, 180.0): None}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('SOL', 4, 'OW'), 2.0, 180.0): {
                    (7, 8, ('SOL', 4, 'HW1'), ('ALA', 5, 'O'), 2.0, 180.0): None}}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('SOL', 4, 'OW'), 2.0, 180.0): {
                    (7, 8, ('SOL', 4, 'HW1'), ('SOL', 5, 'OW'), 2.0, 180.0): {
                        (9, 10, ('SOL', 5, 'HW1'), ('ALA', 6, 'O'), 1.0, 180.0): None}}}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 7, ('SOL', 3, 'HW1'), ('ALA', 4, 'O'), 2.0, 180.0): None,
                (6, 8, ('SOL', 3, 'HW2'), ('ALA', 5, 'O'), 2.0, 180.0): None}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None},
            (5, 7, ('ALA', 5, 'O'), ('SOL', 6, 'HW1'), 2.0, 180.0): {
                (8, 9, ('SOL', 6, 'HW2'), ('SOL', 7, 'OW'), 2.0, 180.0): {
                    (10, 11, ('SOL', 7, 'HW1'), ('ALA', 8, 'O'), 2.0, 180.0): None}}})
        def analysis(current, output):
            s1_index, to_index, (s1_resname, s1_resid, s1_name), (to_resname, to_resid, to_name), dist, angle = \
                current[0]
            from_index, s2_index, (from_resname, from_resid, from_name), (s2_resname, s2_resid, s2_name), dist, angle = \
                current[-1]
            key = (s1_resname, s1_resid, s2_resname, s2_resid, len(current)-1)
            output[key] += 1
        result = [(('ALA', 1, 'ALA', 4, 1), 0.6),
                  (('ALA', 1, 'ALA', 4, 2), 0.2),
                  (('ALA', 1, 'ALA', 5, 2), 0.1),
                  (('ALA', 1, 'ALA', 5, 3), 0.1),
                  (('ALA', 1, 'ALA', 6, 4), 0.1),
                  (('ALA', 5, 'ALA', 8, 2), 0.1)]
        assert_equal(sorted(wb.count_by_type(analysis_func=analysis)), result)

    def test_count_by_time(self):
        '''
        This test tests if count_by_times() works.
        :return:
        '''
        grofile = '''Test gro file
5
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    2SOL     OW    3   0.300   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
 1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)', order=4)

        # Build an dummy WaterBridgeAnalysis object for testing
        wb._network = []
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (1, 3, ('SOL', 2, 'OW'), ('ALA', 4, 'H'), 2.0, 180.0): None}})
        wb._network.append({(1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 180.0): {
            (2, 3, ('SOL', 2, 'OW'), ('ALA', 4, 'H'), 2.0, 180.0): None}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('ALA', 4, 'O'), 2.0, 180.0): None}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('SOL', 4, 'OW'), 2.0, 180.0): {
                    (7, 8, ('SOL', 4, 'HW1'), ('ALA', 5, 'O'), 2.0, 180.0): None}}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 6, ('SOL', 3, 'HW1'), ('SOL', 4, 'OW'), 2.0, 180.0): {
                    (7, 8, ('SOL', 4, 'HW1'), ('SOL', 5, 'OW'), 2.0, 180.0): {
                        (9, 10, ('SOL', 5, 'HW1'), ('ALA', 6, 'O'), 1.0, 180.0): None}}}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW'), 2.0, 180.0): {
                (5, 7, ('SOL', 3, 'HW1'), ('ALA', 4, 'O'), 2.0, 180.0): None,
                (6, 8, ('SOL', 3, 'HW2'), ('ALA', 5, 'O'), 2.0, 180.0): None}}})
        wb._network.append({(0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1'), 2.0, 180.0): {
            (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 180.0): None},
            (5, 7, ('ALA', 5, 'O'), ('SOL', 6, 'HW1'), 2.0, 180.0): {
                (8, 9, ('SOL', 6, 'HW2'), ('SOL', 7, 'OW'), 2.0, 180.0): {
                    (10, 11, ('SOL', 7, 'HW1'), ('ALA', 8, 'O'), 2.0, 180.0): None}}})
        wb.timesteps = range(len(wb._network))
        assert_equal(wb.count_by_time(), [(0,1), (1,1), (2,1), (3,1), (4,1), (5,1), (6,1), (7,1), (8,2), (9,2)])


    def test_count_by_time_weight(self):
        '''
        This test tests if modyfing the analysis_func allows the weight to be changed
        in count_by_type().
        :return:
        '''
        grofile = '''Test gro file
12
    1ALA      O    1   0.000   0.000   0.000
    2SOL     OW    2   0.300   0.000   0.000
    2SOL    HW1    3   0.200   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
    5ALA      O    6   0.000   1.000   0.000
    6SOL     OW    7   0.300   1.000   0.000
    6SOL    HW1    8   0.200   1.000   0.000
    6SOL    HW2    9   0.400   1.000   0.000
    7SOL     OW   10   0.600   1.000   0.000
    7SOL    HW1   11   0.700   1.000   0.000
    8ALA      O   12   0.900   1.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1 or resid 5)', 'protein and (resid 4 or resid 8)', order=2)
        wb.run(verbose=False)
        def analysis(current, output):
            s1_index, to_index, (s1_resname, s1_resid, s1_name), (to_resname, to_resid, to_name), dist, angle = \
                current[0]
            from_index, s2_index, (from_resname, from_resid, from_name), (s2_resname, s2_resid, s2_name), dist, angle = \
                current[-1]
            key = (s1_resname, s1_resid, s2_resname, s2_resid)
            output[key] += len(current)-1
        assert_equal(wb.count_by_time(analysis_func=analysis), [(0,3), ])

    def test_count_by_time_empty(self):
        '''
        See if count_by_type() can handle zero well.
        :return:
        '''
        grofile = '''Test gro file
12
    1ALA      O    1   0.000   0.000   0.000
    2SOL     OW    2   0.300   0.000   0.000
    2SOL    HW1    3   0.200   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
    5ALA      O    6   0.000   1.000   0.000
    6SOL     OW    7   0.300   1.000   0.000
    6SOL    HW1    8   0.200   1.000   0.000
    6SOL    HW2    9   0.400   1.000   0.000
    7SOL     OW   10   0.600   1.000   0.000
    7SOL    HW1   11   0.700   1.000   0.000
    8ALA      O   12   0.900   1.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1 or resid 5)', 'protein and (resid 4 or resid 8)', order=2)
        wb.run(verbose=False)
        def analysis(current, output):
            pass
        assert_equal(wb.count_by_time(analysis_func=analysis), [(0,0), ])

def guess_types(names):
    """GRO doesn't supply types, this returns an Attr"""
    return Atomtypes(np.array([guess_atom_type(name) for name in names], dtype=object))


class TestHydrogenBondAnalysis(object):
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(PDB_helix)

    @staticmethod
    @pytest.fixture(scope='class')
    def values(universe):
        return {
            'num_bb_hbonds': universe.atoms.n_residues - universe.select_atoms('resname PRO').n_residues - 4,
            'donor_resid': np.array([5, 6, 8, 9, 10, 11, 12, 13]),
            'acceptor_resnm': np.array(['ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'PRO', 'ALA', 'ALA'], dtype='U4'),
        }

    kwargs = {
        'selection1': 'protein',
        'selection2': 'protein',
        'detect_hydrogens': "distance",
        'distance': 3.0,
        'angle': 150.0,
    }

    @pytest.fixture(scope='class')
    def h(self, universe):
        kw = self.kwargs.copy()
        # kw.update(kwargs)
        h = MDAnalysis.analysis.hbonds.WaterBridgeAnalysis(universe, order=0, **kw)
        # remove in 1.0
        if kw['detect_hydrogens'] == 'heuristic':
            with pytest.warns(DeprecationWarning):
                h.run(verbose=False)
        else:
            h.run(verbose=False)
        return h

    def test_helix_backbone(self, values, h):
        assert len(h.timeseries[0]) == values['num_bb_hbonds'], "wrong number of backbone hydrogen bonds"
        assert h.timesteps, [0.0]

    def test_generate_table(self, values, h):

        h.generate_table()
        assert len(h.table) == values['num_bb_hbonds'], "wrong number of backbone hydrogen bonds in table"
        assert isinstance(h.table, np.core.records.recarray)

        assert_array_equal(sorted(h.table.donor_resid), values['donor_resid'])
        assert_array_equal(h.table.acceptor_resnm, values['acceptor_resnm'])

    def test_atoms_too_far(self):
        pdb = '''
ATOM      1  N   LEU     1      32.310  13.778  14.372  1.00  0.00      SYST N 0
ATOM      2  OW  SOL     2       3.024   4.456   4.147  1.00  0.00      SYST H 0'''

        u = MDAnalysis.Universe(StringIO(pdb), format="pdb")
        h = WaterBridgeAnalysis(u, 'resname SOL', 'protein', order=0)
        h.run(verbose=False)
        assert h.timeseries == [[]]

    def test_acceptor_OC1_OC2(self):
        gro = '''test
3
    1ALA    OC1    1   0.000   0.000   0.000
    2ALA      N    2   0.300   0.000   0.000
    2ALA     H1    3   0.200   0.000   0.000
7.29748 7.66094 9.82962'''

        u = MDAnalysis.Universe(StringIO(gro), format="gro")
        h = WaterBridgeAnalysis(u, 'protein', 'protein', order=0)
        h.run(verbose=False)
        assert h.timeseries[0][0][2] == 'ALA2:H1'

    def test_true_traj(self):
        u = MDAnalysis.Universe(GRO, XTC)
        u.add_TopologyAttr(guess_types(u.atoms.names))
        h = WaterBridgeAnalysis(u, 'protein', 'resname ASP', distance=3.0, angle=120.0, order=0)
        h.run()
        assert len(h.timeseries) == 10

    def test_count_by_time(self, values, h):

        c = h.count_by_time()
        assert c, [(0.0, values['num_bb_hbonds'])]

    def test_count_by_type(self, values, h):

        c = h.count_by_type()
        assert_equal([line[-1] for line in c][:8], values['num_bb_hbonds'] * [1.0])

    def test_timesteps_by_type(self, values, h):

        t = h.timesteps_by_type()
        assert_equal([i[-1] for i in t][:8], values['num_bb_hbonds'] * [0.0])


class TestHydrogenBondAnalysisPBC(TestHydrogenBondAnalysis):
    # This system is identical to above class
    # But has a box introduced, and atoms moved into neighbouring images
    # The results however should remain identical if PBC is used
    # If pbc:True in kwargs is changed, these tests should fail
    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        u = MDAnalysis.Universe(PDB_helix)
        # transfer to memory to changes to coordinates are reset
        u.transfer_to_memory()
        # place in huge oversized box
        # real system is about 30A wide at most
        boxsize = 150.
        box = np.array([boxsize, boxsize, boxsize, 90., 90., 90.])
        u.dimensions = box

        # then scatter the atoms throughout various images of this box
        u.atoms[::4].translate([boxsize * 2, 0, 0])
        u.atoms[1::4].translate([0, boxsize * 4, 0])
        u.atoms[2::4].translate([-boxsize * 5, 0, -boxsize * 2])

        return u

    kwargs = {
        'selection1': 'protein',
        'selection2': 'protein',
        'detect_hydrogens': "distance",
        'distance': 3.0,
        'angle': 150.0,
        'pbc': True,
    }

class TestHydrogenBondAnalysisTIP3P(object):
    @staticmethod
    @pytest.fixture()
    def universe():
        return MDAnalysis.Universe(waterPSF, waterDCD)

    kwargs = {
        'selection1': 'all',
        'selection2': 'all',
        'detect_hydrogens': "distance",
        'distance': 3.0,
        'angle': 120.0,
    }

    @pytest.fixture()
    def h(self, universe):
        h = WaterBridgeAnalysis(universe, order=0, **self.kwargs)
        h.run(verbose=False)
        h.generate_table()
        return h

    @pytest.fixture()
    def normalized_timeseries(self, h):
        # timeseries in normalized form: (t, d_indx1, a_indx1, d_index0, a_index0, donor, acceptor, dist, angle)
        #                   array index:  0     1        2        3         4        5      6        7      8
        timeseries = [[t] + item
                      for t, hframe in zip(h.timesteps, h.timeseries)
                      for item in hframe]
        return timeseries

    # keys are the names in the h.table
    reference = {
        'distance': {'mean': 2.0208776, 'std': 0.31740859},
        'angle': {'mean': 155.13521, 'std': 12.98955},
    }

    @pytest.fixture()
    def reference_table(self, normalized_timeseries):
        # reference values for the table only
        return {
            'donor_resnm': ["TIP3"] * len(normalized_timeseries),
            'acceptor_resnm': ["TIP3"] * len(normalized_timeseries),
        }

    # index into timeseries (ADJUST ONCE donor_idx and acceptor_ndx are removed)
    # with keys being field names in h.table
    columns = {
        'time': 0,
        'donor_index': 1,
        'acceptor_index': 2,
        'distance': 5,
        'angle': 6,
    }

    # hackish way to allow looping over self.reference and generating tests
    _functions = {
        'mean': np.mean,
        'std': np.std,
    }

    def test_timeseries(self, h, normalized_timeseries):
        h = h
        assert len(h.timeseries) == 10
        assert len(normalized_timeseries) == 29

        for observable in self.reference:
            idx = self.columns[observable]
            for quantity, reference in self.reference[observable].items():
                func = self._functions[quantity]
                assert_allclose(
                    func([item[idx] for item in normalized_timeseries]), reference,
                    rtol=1e-5, atol=0,
                    err_msg="{quantity}({observable}) does not match reference".format(**vars())
                )

    def test_table_atoms(self, h, normalized_timeseries, reference_table):
        h = h
        table = h.table

        assert len(h.table) == len(normalized_timeseries)

        # test that timeseries and table agree on index data and
        # hydrogen bond information at atom level
        for name, idx in self.columns.items():
            assert_array_almost_equal(table.field(name), [data[idx] for data in normalized_timeseries],
                                      err_msg="table[{name}] and timeseries[{idx} do not agree".format(**vars()))

        # test at residue level (issue #801
        # https://github.com/MDAnalysis/mdanalysis/issues/801)
        for name, ref in reference_table.items():
            assert_array_equal(h.table.field(name), ref, err_msg="resname for {0} do not match (Issue #801)")