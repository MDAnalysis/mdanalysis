from __future__ import print_function, absolute_import
from six import StringIO

from numpy.testing import assert_equal

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
        timeseries = wb._timeseries
        assert_equal(timeseries[0][0][:4], (2, 0, ('SOL', 2, 'HW1'), ('ALA', 1, 'O')))
        assert_equal(timeseries[0][1][:4], (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O')))

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
        timeseries = wb._timeseries
        assert_equal(timeseries[0][0][:4], (1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW')))
        assert_equal(timeseries[0][1][:4], (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O')))

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
        timeseries = wb._timeseries
        assert_equal(timeseries[0][0][:4], (2, 0, ('SOL', 2, 'HW1'), ('ALA', 1, 'O')))
        assert_equal(timeseries[0][1][:4], (3, 1, ('ALA', 4, 'H'), ('SOL', 2, 'OW')))

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
        timeseries = wb._timeseries
        assert_equal(timeseries[0][0][:4], (1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW')))
        assert_equal(timeseries[0][1][:4], (3, 2, ('ALA', 4, 'H'), ('SOL', 2, 'OW')))

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
        assert_equal(wb._timeseries, [[]])

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
        assert_equal(wb._timeseries, [[]])

    def test_water_network(self):
        '''
        This test tests if the internal water object is generated correctly.
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
        water_network = wb._water_network[0]
        assert_equal(list(water_network.keys()), [('SOL', 2)])
        values = list(water_network.values())
        assert_equal(list(values[0][0])[0][:4], (1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW')))
        assert_equal(list(values[0][1])[0][:4], (3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O')))

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
        assert_equal(wb.count_by_type().tolist(), [(1, 4, 'ALA', 1, 'H', 'ALA', 4, 'O',  1.)])

    def test_count_by_type_multiple_link(self):
        '''
        This test tests if count_by_type() can assemble linkage from water network.
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
        # Build an dummy WaterBridgeAnalysis object for testing
        wb._timeseries = True
        wb.timesteps = [0]
        wb._water_network = [{('SOL', 2): [{(2, 0, ('SOL', 2, 'HW1'), ('ALA', 1, 'O'), 2.0, 179.99998),
                                            (1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 179.99998)},
                                           {(3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 179.99998),
                                            (3, 2, ('ALA', 4, 'H'), ('SOL', 2, 'OW'), 2.0, 179.99998)}]}]
        result = [(0, 3, 'ALA', 1, 'O', 'ALA', 4, 'H', 1.0),
                  (0, 4, 'ALA', 1, 'O', 'ALA', 4, 'O', 1.0),
                  (1, 3, 'ALA', 1, 'H', 'ALA', 4, 'H', 1.0),
                  (1, 4, 'ALA', 1, 'H', 'ALA', 4, 'O', 1.0)]
        assert_equal(sorted(wb.count_by_type().tolist()), result)

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
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)')
        # Build an dummy WaterBridgeAnalysis object for testing
        wb._timeseries = True
        wb.timesteps = [0, 1, 2, 3, 4, 5]
        wb._water_network = [# a 2 * 2 water netwrok consists of all four links
                             {('SOL', 2): [{(2, 0, ('SOL', 2, 'HW1'), ('ALA', 1, 'O'), 2.0, 179.99998),
                                            (1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 179.99998)},
                                           {(3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 179.99998),
                                            (3, 2, ('ALA', 4, 'H'), ('SOL', 2, 'OW'), 2.0, 179.99998)}]},
                             # a repeat
                             {('SOL', 2): [{(2, 0, ('SOL', 2, 'HW1'), ('ALA', 1, 'O'), 2.0, 179.99998),
                                            (1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 179.99998)},
                                           {(3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 179.99998),
                                            (3, 2, ('ALA', 4, 'H'), ('SOL', 2, 'OW'), 2.0, 179.99998)}]},
                             # single link 1
                             {('SOL', 2): [{(2, 0, ('SOL', 2, 'HW1'), ('ALA', 1, 'O'), 2.0, 179.99998)},
                                           {(3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 179.99998)}]},
                             # single link 2
                             {('SOL', 2): [{(2, 0, ('SOL', 2, 'HW1'), ('ALA', 1, 'O'), 2.0, 179.99998)},
                                           {(3, 1, ('ALA', 4, 'H'), ('SOL', 2, 'OW'), 2.0, 179.99998)}]},
                             # single link 3
                             {('SOL', 2): [{(1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 179.99998)},
                                           {(3, 4, ('SOL', 2, 'HW2'), ('ALA', 4, 'O'), 2.0, 179.99998)}]},
                             # single link 4
                             {('SOL', 2): [{(1, 2, ('ALA', 1, 'H'), ('SOL', 2, 'OW'), 2.0, 179.99998)},
                                           {(3, 2, ('ALA', 4, 'H'), ('SOL', 2, 'OW'), 2.0, 179.99998)}]}]
        result = [(0, 3, 'ALA', 1, 'O', 'ALA', 4, 'H', 0.5),
                  (0, 4, 'ALA', 1, 'O', 'ALA', 4, 'O', 0.5),
                  (1, 3, 'ALA', 1, 'H', 'ALA', 4, 'H', 0.5),
                  (1, 4, 'ALA', 1, 'H', 'ALA', 4, 'O', 0.5)]
        assert_equal(sorted(wb.count_by_type().tolist()), result)
