from __future__ import print_function, absolute_import
from six import StringIO
from collections import defaultdict

from numpy.testing import assert_equal, assert_almost_equal

import MDAnalysis
from MDAnalysis.analysis.hbonds.wbridge_analysis import WaterBridgeAnalysis

def test_import_from_hbonds():
    try:
        from MDAnalysis.analysis.hbonds import WaterBridgeAnalysis
    except ImportError:
        raise AssertionError("Issue #2064 not fixed: "
                             "importing WaterBridgeAnalysis from "
                             "MDAnalysis.analysis.hbonds failed.'")

class TestWaterBridgeAnalysis(object):
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
        timeseries = wb.timeseries[0]
        assert_equal(timeseries[0][0][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        assert_equal(timeseries[0][1][:4], (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW')))
        assert_equal(timeseries[1][0][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        assert_equal(timeseries[1][1][:4], (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW')))
        assert_equal([(5, 7, ('SOL', 3, 'HW1'), ('ALA', 4, 'O')), (6, 8, ('SOL', 3, 'HW2'), ('ALA', 5, 'O'))],
                     sorted([line[2][:4] for line in timeseries]))

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
        result = [(0, 3, 'ALA', 1, 'O', 'ALA', 4, 'H', 0.1),
                  (1, 3, 'ALA', 1, 'H', 'ALA', 4, 'H', 0.1)]
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
        result = [('ALA', 1, 'ALA', 4, 0.8),
                  ('ALA', 1, 'ALA', 5, 0.2),
                  ('ALA', 1, 'ALA', 6, 0.1),
                  ('ALA', 5, 'ALA', 8, 0.1)]
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
        result = [('ALA', 1, 'ALA', 4, 1, 0.6),
                  ('ALA', 1, 'ALA', 4, 2, 0.2),
                  ('ALA', 1, 'ALA', 5, 2, 0.1),
                  ('ALA', 1, 'ALA', 5, 3, 0.1),
                  ('ALA', 1, 'ALA', 6, 4, 0.1),
                  ('ALA', 5, 'ALA', 8, 2, 0.1)]
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
        assert_equal(wb.count_by_time(), [1, 1, 1, 1, 1, 1, 1, 1, 2, 2])


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
        assert_equal(wb.count_by_time(analysis_func=analysis), [3, ])

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
        assert_equal(wb.count_by_time(analysis_func=analysis), [0, ])