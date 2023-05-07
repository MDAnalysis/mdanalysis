from io import StringIO
from collections import defaultdict

from numpy.testing import (
    assert_equal, assert_array_equal,)
import pytest

import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.wbridge_analysis import (
    WaterBridgeAnalysis, )


class TestWaterBridgeAnalysis(object):
    @staticmethod
    @pytest.fixture(scope='class')
    def universe_empty():
        '''A universe with no hydrogen bonds'''
        grofile = '''Test gro file
5
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    2SOL     OW    3   3.000   0.000   0.000
    4ALA      H    4   0.500   0.000   0.000
    4ALA      N    5   0.600   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def universe_DA():
        '''A universe with one hydrogen bond acceptor bonding to a hydrogen bond
        donor'''
        grofile = '''Test gro file
3
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    4ALA      O    3   0.300   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def universe_DA_PBC():
        '''A universe with one hydrogen bond acceptor bonding to a hydrogen bond
        donor but in a PBC condition'''
        grofile = '''Test gro file
3
    1ALA      N    1   0.800   0.000   0.000
    1ALA      H    2   0.900   0.000   0.000
    4ALA      O    3   0.100   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def universe_AD():
        '''A universe with one hydrogen bond donor bonding to a hydrogen bond
        acceptor'''
        grofile = '''Test gro file
3
    1ALA      O    1   0.000   0.000   0.000
    4ALA      H    2   0.200   0.000   0.000
    4ALA      N    3   0.300   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def universe_loop():
        '''A universe with one hydrogen bond acceptor bonding to a water which
        bonds back to the first hydrogen bond  acceptor and thus form a loop'''
        grofile = '''Test gro file
5
    1ALA      O    1   0.000   0.001   0.000
    2SOL     OW    2   0.300   0.001   0.000
    2SOL    HW1    3   0.200   0.002   0.000
    2SOL    HW2    4   0.200   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def universe_DWA():
        '''A universe with one hydrogen bond donor bonding to a hydrogen bond
        acceptor through a water'''
        grofile = '''Test gro file
5
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    2SOL     OW    3   0.300   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def universe_DWD():
        '''A universe with one hydrogen bond donor bonding to a hydrogen bond
        donor through a water'''
        grofile = '''Test gro file
5
    1ALA      N    1   0.000   0.000   0.000
    1ALA      H    2   0.100   0.000   0.000
    2SOL     OW    3   0.300   0.000   0.000
    4ALA      H    4   0.500   0.000   0.000
    4ALA      N    5   0.600   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def universe_AWA():
        '''A universe with two hydrogen bond acceptor are joined by a water'''
        grofile = '''Test gro file
5
    1ALA      O    1   0.000   0.000   0.000
    2SOL     OW    2   0.300   0.000   0.000
    2SOL    HW1    3   0.200   0.000   0.000
    2SOL    HW2    4   0.400   0.000   0.000
    4ALA      O    5   0.600   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def universe_AWD():
        '''A universe with one hydrogen bond acceptor bonding to a hydrogen
        bond donor through a water'''
        grofile = '''Test gro file
5
    1ALA      O    1   0.000   0.000   0.000
    2SOL     OW    2   0.300   0.000   0.000
    2SOL    HW1    3   0.200   0.000   0.000
    4ALA      H    4   0.500   0.000   0.000
    4ALA      N    5   0.600   0.000   0.000
  1.0   1.0   1.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def universe_AWWA():
        '''A universe with one hydrogen bond acceptor bonding to a hydrogen bond
        acceptor through two waters'''
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
        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def universe_AWWWA():
        '''A universe with one hydrogen bond acceptor bonding to a hydrogen bond
        acceptor through three waters'''
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
        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def universe_AWWWWA():
        '''A universe with one hydrogen bond acceptor bonding to a hydrogen bond
        acceptor through three waters'''
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
        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def universe_branch():
        '''A universe with one hydrogen bond acceptor bonding to two hydrogen
        bond acceptor in selection 2'''
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
        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def universe_AWA_AWWA():
        '''A universe with one hydrogen bond acceptors are bonded through one or
        two water'''
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
        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def wb_multiframe():
        '''A water bridge object with multipley frames'''
        grofile = '''Test gro file
 13
    1ALA      O    1   0.000   0.000   0.000
    1ALA      H    2   0.000   0.000   0.000
    2SOL     OW    3   0.300   0.000   0.000
    2SOL    HW1    4   0.200   0.000   0.000
    2SOL    HW2    5   0.400   0.000   0.000
    3SOL     OW    6   0.600   0.000   0.000
    3SOL    HW1    7   0.700   0.000   0.000
    4SOL     OW    8   0.900   0.000   0.000
    4SOL    HW1    9   1.000   0.000   0.000
    5SOL     OW   10   1.200   0.000   0.000
    5SOL    HW1   11   1.300   0.000   0.000
    6ALA      H   12   1.400   0.000   0.000
    6ALA      O   13   1.400   0.000   0.000
   10.0   10.0   10.0'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'protein and (resid 1)', 'protein and (resid 4)',
                                 order=4)
        # Build an dummy WaterBridgeAnalysis object for testing
        wb.results.network = []
        wb.results.network.append({(1, 0, 12, None, 2.0, 180.0): None})
        wb.results.network.append({(0, None, 12, 13, 2.0, 180.0): None})
        wb.results.network.append({(1, 0, 3, None, 2.0, 180.0):
                            {(4, 2, 12, None, 2.0, 180.0): None}})
        wb.results.network.append({(0, None, 3, 2, 2.0, 180.0):
                            {(4, 2, 5, None, 2.0, 180.0):
                             {(5, None, 11, 12, 2.0, 180.0): None}}})
        wb.timesteps = range(len(wb.results.network))
        return wb

    def test_nodata(self, universe_DA):
        '''Test if the funtions can run when there is no data.
        This is achieved by not runing the run() first.'''
        wb = WaterBridgeAnalysis(universe_DA, 'protein and (resid 1)',
                                 'protein and (resid 4)', order=0)
        wb.generate_table()
        assert_equal(wb.timesteps_by_type(), None)
        assert_equal(wb.count_by_time(), None)
        assert_equal(wb.count_by_type(), None)

    def test_selection_type_error(self, universe_DA):
        '''Test the case when the wrong selection1_type is given'''
        try:
            wb = WaterBridgeAnalysis(universe_DA, 'protein and (resid 1)',
                'protein and (resid 4)', order=0, selection1_type='aaa')
        except ValueError:
            pass
        else:
            raise pytest.fail("selection_type aaa should rasie error")

    def test_distance_type_error(self, universe_DA):
        '''Test the case when the wrong selection1_type is given'''
        with pytest.raises(ValueError, match="Only 'hydrogen' and 'heavy' are allowed for option `distance_type'"):
            WaterBridgeAnalysis(universe_DA, 'protein and (resid 1)',
                                'protein and (resid 4)', order=0,
                                selection1_type='aaa', distance_type='aaa')

    def test_selection2_type_error(self, universe_DA):
        '''Test the case when the wrong selection1_type is given'''
        with pytest.raises(ValueError, match="`selection2_type` is not a keyword argument."):
            WaterBridgeAnalysis(universe_DA, 'protein and (resid 1)',
                                'protein and (resid 4)', order=0,
                                selection1_type='aaa', selection2_type='aaa')


    def test_empty_selection(self, universe_DA):
        '''Test the case when selection yields empty result'''
        wb = WaterBridgeAnalysis(universe_DA, 'protein and (resid 9)',
                                 'protein and (resid 10)', order=0)
        wb.run()
        assert wb.results.network == [{}]

    def test_loop(self, universe_loop):
        '''Test if loop can be handled correctly'''
        wb = WaterBridgeAnalysis(universe_loop, 'protein and (resid 1)',
                                 'protein and (resid 1 or resid 4)')
        wb.run()
        assert_equal(len(wb.results.network[0].keys()), 2)

    @pytest.mark.parametrize('distance_type', ["hydrogen", "heavy"])
    def test_donor_accepter(self, universe_DA, distance_type):
        '''Test zeroth order donor to acceptor hydrogen bonding'''
        wb = WaterBridgeAnalysis(universe_DA, 'protein and (resid 1)',
                                 'protein and (resid 4)',
                                 order=0,
                                 update_selection=True,
                                 debug=True,
                                 distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (1, 0, 2, None))

    @pytest.mark.parametrize('distance_type', ["hydrogen", "heavy"])
    def test_donor_accepter_pbc(self, universe_DA_PBC, distance_type):
        '''Test zeroth order donor to acceptor hydrogen bonding in PBC conditions'''
        wb = WaterBridgeAnalysis(universe_DA_PBC,
                                 'protein and (resid 1)',
                                 'protein and (resid 4)',
                                 order=0,
                                 pbc=True,
                                 distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (1, 0, 2, None))

    @pytest.mark.parametrize('distance_type', ["hydrogen", "heavy"])
    def test_accepter_donor(self, universe_AD, distance_type):
        '''Test zeroth order acceptor to donor hydrogen bonding'''
        wb = WaterBridgeAnalysis(universe_AD, 'protein and (resid 1)',
                                 'protein and (resid 4)', order=0,
                                 distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 1, 2))

    @pytest.mark.parametrize('distance_type', ["hydrogen", "heavy"])
    def test_acceptor_water_accepter(self, universe_AWA, distance_type):
        '''Test case where the hydrogen bond acceptor from selection 1 form
        water bridge with hydrogen bond acceptor from selection 2'''
        wb = WaterBridgeAnalysis(universe_AWA, 'protein and (resid 1)',
                                 'protein and (resid 4)', distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        assert_equal(second[list(second.keys())[0]], None)

    @pytest.mark.parametrize('distance_type', ["hydrogen", "heavy"])
    def test_donor_water_accepter(self, universe_DWA, distance_type):
        '''Test case where the hydrogen bond donor from selection 1 form
        water bridge with hydrogen bond acceptor from selection 2'''
        wb = WaterBridgeAnalysis(universe_DWA, 'protein and (resid 1)',
                                 'protein and (resid 4)', distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (1, 0, 2, None))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 2, 4, None))
        assert_equal(second[list(second.keys())[0]], None)

    @pytest.mark.parametrize('distance_type', ["hydrogen", "heavy"])
    def test_acceptor_water_donor(self, universe_AWD, distance_type):
        '''Test case where the hydrogen bond acceptor from selection 1 form
        water bridge with hydrogen bond donor from selection 2'''
        wb = WaterBridgeAnalysis(universe_AWD, 'protein and (resid 1)',
                                 'protein and (resid 4)', distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (1, None, 3, 4))
        assert_equal(second[list(second.keys())[0]], None)

    @pytest.mark.parametrize('distance_type', ["hydrogen", "heavy"])
    def test_donor_water_donor(self, universe_DWD, distance_type):
        '''Test case where the hydrogen bond donor from selection 1 form
        water bridge with hydrogen bond donor from selection 2'''
        wb = WaterBridgeAnalysis(universe_DWD, 'protein and (resid 1)',
                                 'protein and (resid 4)', distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (1, 0, 2, None))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (2, None, 3, 4))
        assert_equal(second[list(second.keys())[0]], None)

    def test_empty(self, universe_empty):
        '''Test case where no water bridge exists'''
        wb = WaterBridgeAnalysis(universe_empty, 'protein', 'protein')
        wb.run(verbose=False)
        assert_equal(wb.results.network[0], defaultdict(dict))

    def test_same_selection(self, universe_DWA):
        '''
        This test tests that if the selection 1 and selection 2 are both protein.
        However, the protein only forms one hydrogen bond with the water.
        This entry won't be included.
        '''
        wb = WaterBridgeAnalysis(universe_DWA, 'protein and resid 1',
                                 'protein and resid 1')
        wb.run(verbose=False)
        assert_equal(wb.results.network[0], defaultdict(dict))

    @pytest.mark.parametrize('distance_type', ["hydrogen", "heavy"])
    def test_acceptor_2water_accepter(self, universe_AWWA, distance_type):
        '''Test case where the hydrogen bond acceptor from selection 1 form second order
        water bridge with hydrogen bond acceptor from selection 2'''
        # test first order
        wb = WaterBridgeAnalysis(universe_AWWA, 'protein and (resid 1)',
                                 'protein and (resid 4)',
                                 distance_type=distance_type)
        wb.run(verbose=False)
        assert_equal(wb.results.network[0], defaultdict(dict))
        # test second order
        wb = WaterBridgeAnalysis(universe_AWWA,
                                 'protein and (resid 1)',
                                 'protein and (resid 4)',
                                 order=2,
                                 distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 4, 6, None))
        assert_equal(third[list(third.keys())[0]], None)
        # test third order
        wb = WaterBridgeAnalysis(universe_AWWA, 'protein and (resid 1)',
                                 'protein and (resid 4)', order=3,
                                 distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 4, 6, None))
        assert_equal(third[list(third.keys())[0]], None)

    @pytest.mark.parametrize('distance_type', ["hydrogen", "heavy"])
    def test_acceptor_3water_accepter(self, universe_AWWWA, distance_type):
        '''Test case where the hydrogen bond acceptor from selection 1 form third order
        water bridge with hydrogen bond acceptor from selection 2'''
        wb = WaterBridgeAnalysis(universe_AWWWA, 'protein and (resid 1)',
                                 'protein and (resid 5)', order=2,
                                 distance_type=distance_type)
        wb.run(verbose=False)
        assert_equal(wb.results.network[0], defaultdict(dict))

        wb = WaterBridgeAnalysis(universe_AWWWA, 'protein and (resid 1)',
                                 'protein and (resid 5)', order=3,
                                 distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 4, 6, None))
        fourth = third[list(third.keys())[0]]
        assert_equal(list(fourth.keys())[0][:4], (7, 6, 8, None))
        assert_equal(fourth[list(fourth.keys())[0]], None)

        wb = WaterBridgeAnalysis(universe_AWWWA, 'protein and (resid 1)',
                                 'protein and (resid 5)', order=4,
                                 distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 4, 6, None))
        fourth = third[list(third.keys())[0]]
        assert_equal(list(fourth.keys())[0][:4], (7, 6, 8, None))
        assert_equal(fourth[list(fourth.keys())[0]], None)

    @pytest.mark.parametrize('distance_type', ["hydrogen", "heavy"])
    def test_acceptor_4water_accepter(self, universe_AWWWWA, distance_type):
        '''Test case where the hydrogen bond acceptor from selection 1 form fourth order
        water bridge with hydrogen bond acceptor from selection 2'''
        wb = WaterBridgeAnalysis(universe_AWWWWA, 'protein and (resid 1)',
                                 'protein and (resid 6)', order=3,
                                 distance_type=distance_type)
        wb.run(verbose=False)
        assert_equal(wb.results.network[0], defaultdict(dict))

        wb = WaterBridgeAnalysis(universe_AWWWWA, 'protein and (resid 1)',
                                 'protein and (resid 6)', order=4,
                                 distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 4, 6, None))
        fourth = third[list(third.keys())[0]]
        assert_equal(list(fourth.keys())[0][:4], (7, 6, 8, None))
        fifth = fourth[list(fourth.keys())[0]]
        assert_equal(list(fifth.keys())[0][:4], (9, 8, 10, None))
        assert_equal(fifth[list(fifth.keys())[0]], None)

        wb = WaterBridgeAnalysis(universe_AWWWWA, 'protein and (resid 1)',
                                 'protein and (resid 6)', order=5,
                                 distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 4, 6, None))
        fourth = third[list(third.keys())[0]]
        assert_equal(list(fourth.keys())[0][:4], (7, 6, 8, None))
        fifth = fourth[list(fourth.keys())[0]]
        assert_equal(list(fifth.keys())[0][:4], (9, 8, 10, None))
        assert_equal(fifth[list(fifth.keys())[0]], None)

    @pytest.mark.parametrize('distance_type', ["hydrogen", "heavy"])
    def test_acceptor_22water_accepter(self, universe_branch, distance_type):
        '''Test case where the hydrogen bond acceptor from selection 1 form a second order
        water bridge with hydrogen bond acceptor from selection 2
        and the last water is linked to two residues in selection 2'''
        wb = WaterBridgeAnalysis(universe_branch, 'protein and (resid 1)',
                                 'protein and (resid 4 or resid 5)', order=2,
                                 distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        third = second[list(second.keys())[0]]
        assert_equal([(5, 4, 7, None), (6, 4, 8, None)],
                     sorted([key[:4] for key in list(third.keys())]))

    def test_timeseries_wba(self, universe_branch):
        '''Test if the time series data is correctly generated in water bridge analysis format'''
        wb = WaterBridgeAnalysis(universe_branch, 'protein and (resid 1)',
                                 'protein and (resid 4 or resid 5)', order=2)
        wb.output_format = 'sele1_sele2'
        wb.run(verbose=False)
        timeseries = sorted(wb.results.timeseries[0])

        assert_equal(timeseries[0][:4], (0, 2, ('ALA', 1, 'O'), ('SOL', 2, 'HW1')))
        assert_equal(timeseries[1][:4], (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW')))
        assert_equal(timeseries[2][:4], (5, 7, ('SOL', 3, 'HW1'), ('ALA', 4, 'O')))
        assert_equal(timeseries[3][:4], (6, 8, ('SOL', 3, 'HW2'), ('ALA', 5, 'O')))

    def test_timeseries_hba(self, universe_branch):
        '''Test if the time series data is correctly generated in hydrogen bond analysis format'''
        wb = WaterBridgeAnalysis(universe_branch, 'protein and (resid 1)',
                                 'protein and (resid 4 or resid 5)', order=2)
        wb.output_format = 'donor_acceptor'
        wb.run(verbose=False)
        timeseries = sorted(wb.results.timeseries[0])

        assert_equal(timeseries[0][:4], (2, 0, ('SOL', 2, 'HW1'), ('ALA', 1, 'O')))
        assert_equal(timeseries[1][:4], (3, 4, ('SOL', 2, 'HW2'), ('SOL', 3, 'OW')))
        assert_equal(timeseries[2][:4], (5, 7, ('SOL', 3, 'HW1'), ('ALA', 4, 'O')))
        assert_equal(timeseries[3][:4], (6, 8, ('SOL', 3, 'HW2'), ('ALA', 5, 'O')))

    @pytest.mark.parametrize('distance_type', ["hydrogen", "heavy"])
    def test_acceptor_12water_accepter(self, universe_AWA_AWWA, distance_type):
        '''Test of independent first order and second can be recognised correctely'''
        wb = WaterBridgeAnalysis(universe_AWA_AWWA, 'protein and (resid 1 or resid 5)',
                                 'protein and (resid 4 or resid 8)', order=1,
                                 distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        assert_equal(second[list(second.keys())[0]], None)
        network = wb.results.network[0]
        wb = WaterBridgeAnalysis(universe_AWA_AWWA, 'protein and (resid 1 or resid 5)',
                                 'protein and (resid 4 or resid 8)', order=2,
                                 distance_type=distance_type)
        wb.run(verbose=False)
        network = wb.results.network[0]
        assert_equal([(0, None, 2, 1), (5, None, 7, 6)],
                     sorted([key[:4] for key in list(network.keys())]))

    def test_count_by_type_single_link(self, universe_DWA):
        '''
        This test tests the simplest water bridge to see if count_by_type() works.
        '''
        wb = WaterBridgeAnalysis(universe_DWA, 'protein and (resid 1)',
                                 'protein and (resid 4)')
        wb.run(verbose=False)
        assert_equal(wb.count_by_type(), [(1, 4, 'ALA', 1, 'H', 'ALA', 4, 'O',  1.)])

    def test_count_by_type_multiple_link(self, universe_AWA_AWWA):
        '''
        This test tests if count_by_type() can give the correct result for more than 1 links.
        '''
        wb = WaterBridgeAnalysis(universe_AWA_AWWA, 'protein and (resid 1 or resid 5)',
                                 'protein and (resid 4 or resid 8)', order=2)
        wb.run(verbose=False)
        assert_equal(sorted(wb.count_by_type()),
        [[0, 4, 'ALA', 1, 'O', 'ALA', 4, 'O', 1.0],
         [5, 11, 'ALA', 5, 'O', 'ALA', 8, 'O', 1.0]])


    def test_count_by_type_multiple_frame(self, wb_multiframe):
        '''
        This test tests if count_by_type() works in multiply situations.
        :return:
        '''
        result = [[0, 11, 'ALA', 1, 'O', 'ALA', 6, 'H', 0.25],
                  [0, 12, 'ALA', 1, 'O', 'ALA', 6, 'O', 0.25],
                  [1, 12, 'ALA', 1, 'H', 'ALA', 6, 'O', 0.5]]

        assert_equal(sorted(wb_multiframe.count_by_type()), result)

    def test_count_by_type_filter(self, wb_multiframe):
        '''
        This test tests if modifying analysis_func
        allows some results to be filtered out in count_by_type().
        :return:
        '''
        def analysis(current, output, u):
            sele1_index, sele1_heavy_index, atom2, heavy_atom2, dist, angle = current[0]
            atom1, heavy_atom1, sele2_index, sele2_heavy_index, dist, angle = current[-1]
            sele1 = u.atoms[sele1_index]
            sele2 = u.atoms[sele2_index]
            (s1_resname, s1_resid, s1_name) = (sele1.resname, sele1.resid, sele1.name)
            (s2_resname, s2_resid, s2_name) = (sele2.resname, sele2.resid, sele2.name)

            key = (sele1_index, sele2_index, s1_resname, s1_resid, s1_name, s2_resname, s2_resid, s2_name)
            if s2_name == 'H':
                output[key] += 1
        result = [((0, 11, 'ALA', 1, 'O', 'ALA', 6, 'H'), 0.25)]
        assert_equal(sorted(wb_multiframe.count_by_type(analysis_func=analysis)), result)

    def test_count_by_type_merge(self, wb_multiframe):
        '''
        This test tests if modifying analysis_func
        allows some same residue to be merged in count_by_type().
        '''
        def analysis(current, output, u):
            sele1_index, sele1_heavy_index, atom2, heavy_atom2, dist, angle = current[0]
            atom1, heavy_atom1, sele2_index, sele2_heavy_index, dist, angle = current[-1]
            sele1 = u.atoms[sele1_index]
            sele2 = u.atoms[sele2_index]
            (s1_resname, s1_resid, s1_name) = (sele1.resname, sele1.resid, sele1.name)
            (s2_resname, s2_resid, s2_name) = (sele2.resname, sele2.resid, sele2.name)

            key = (s1_resname, s1_resid, s2_resname, s2_resid)
            output[key] = 1
        result = [(('ALA', 1, 'ALA', 6), 1.0)]
        assert_equal(sorted(wb_multiframe.count_by_type(analysis_func=analysis)), result)

    def test_count_by_type_order(self, wb_multiframe):
        '''
        This test tests if modifying analysis_func
        allows the order of water bridge to be separated in count_by_type().
        :return:
        '''
        def analysis(current, output, u):
            sele1_index, sele1_heavy_index, atom2, heavy_atom2, dist, angle = current[0]
            atom1, heavy_atom1, sele2_index, sele2_heavy_index, dist, angle = current[-1]
            sele1 = u.atoms[sele1_index]
            sele2 = u.atoms[sele2_index]
            (s1_resname, s1_resid, s1_name) = (sele1.resname, sele1.resid, sele1.name)
            (s2_resname, s2_resid, s2_name) = (sele2.resname, sele2.resid, sele2.name)
            key = (s1_resname, s1_resid, s2_resname, s2_resid, len(current)-1)
            output[key] = 1
        result = [(('ALA', 1, 'ALA', 6, 0), 0.5),
                  (('ALA', 1, 'ALA', 6, 1), 0.25),
                  (('ALA', 1, 'ALA', 6, 2), 0.25)]
        assert_equal(sorted(wb_multiframe.count_by_type(analysis_func=analysis)), result)

    def test_count_by_time(self, wb_multiframe):
        '''
        This test tests if count_by_times() works.
        :return:
        '''
        assert_equal(wb_multiframe.count_by_time(), [(0, 1), (1, 1), (2, 1), (3, 1)])


    def test_count_by_time_weight(self, universe_AWA_AWWA):
        '''
        This test tests if modyfing the analysis_func allows the weight to be changed
        in count_by_type().
        :return:
        '''
        wb = WaterBridgeAnalysis(universe_AWA_AWWA, 'protein and (resid 1 or resid 5)',
                                 'protein and (resid 4 or resid 8)', order=2)
        wb.run(verbose=False)
        def analysis(current, output, u):
            sele1_index, sele1_heavy_index, atom2, heavy_atom2, dist, angle = current[0]
            atom1, heavy_atom1, sele2_index, sele2_heavy_index, dist, angle = current[-1]
            sele1 = u.atoms[sele1_index]
            sele2 = u.atoms[sele2_index]
            (s1_resname, s1_resid, s1_name) = (sele1.resname, sele1.resid, sele1.name)
            (s2_resname, s2_resid, s2_name) = (sele2.resname, sele2.resid, sele2.name)
            key = (s1_resname, s1_resid, s2_resname, s2_resid)
            output[key] += len(current)-1
        assert_equal(wb.count_by_time(analysis_func=analysis), [(0,3), ])

    def test_count_by_time_empty(self, universe_AWA_AWWA):
        '''
        See if count_by_time() can handle zero well.
        :return:
        '''
        wb = WaterBridgeAnalysis(universe_AWA_AWWA, 'protein and (resid 1 or resid 5)',
                                 'protein and (resid 4 or resid 8)', order=2)
        wb.run(verbose=False)
        def analysis(current, output, u):
            pass
        assert_equal(wb.count_by_time(analysis_func=analysis), [(0,0), ])

    def test_generate_table_hba(self, wb_multiframe):
        '''Test generate table using hydrogen bond analysis format'''
        table = wb_multiframe.generate_table(output_format='donor_acceptor')
        assert_array_equal(
            sorted(table.donor_resid),
            [1, 1, 2, 2, 2, 6, 6],
        )

    def test_generate_table_s1s2(self, wb_multiframe):
        '''Test generate table using hydrogen bond analysis format'''
        table = wb_multiframe.generate_table(output_format='sele1_sele2')
        assert_array_equal(
            sorted(table.sele1_resid),
            [1, 1, 1, 1, 2, 2, 3],
        )

    def test_timesteps_by_type(self, wb_multiframe):
        '''Test the timesteps_by_type function'''

        timesteps = sorted(wb_multiframe.timesteps_by_type())
        assert_array_equal(timesteps[3], [1, 12, 'ALA', 1, 'H', 'ALA', 6, 'O', 0, 2])

    def test_duplicate_water(self):
        '''A case #3119 where
        Acceptor···H−O···H-Donor
                     |
                     H···O-H
        will be recognised as 3rd order water bridge.
        '''
        grofile = '''Test gro file
    7
    1LEU      O    1   1.876   0.810   1.354
  117SOL    HW1    2   1.853   0.831   1.162
  117SOL     OW    3   1.877   0.890   1.081
  117SOL    HW2    4   1.908   0.828   1.007
  135SOL     OW    5   1.924   0.713   0.845
    1LEU      H    6   1.997   0.991   1.194
    1LEU      N    7   2.041   1.030   1.274
   2.22092   2.22092   2.22092'''
        u = MDAnalysis.Universe(StringIO(grofile), format='gro')
        wb = WaterBridgeAnalysis(u, 'resname LEU and name O',
                                 'resname LEU and name N H', order=4)
        wb.run()
        assert len(wb.results.timeseries[0]) == 2

    def test_warn_results_deprecated(self, universe_DA):
        wb = WaterBridgeAnalysis(universe_DA, 'protein and (resid 9)',
                                 'protein and (resid 10)', order=0)
        wb.run()

        wmsg = "The `network` attribute was deprecated in MDAnalysis 2.0.0"
        with pytest.warns(DeprecationWarning, match=wmsg):
            assert_equal(wb.network, wb.results.network)

        wmsg = "The `timeseries` attribute was deprecated in MDAnalysis 2.0.0"
        with pytest.warns(DeprecationWarning, match=wmsg):
            assert_equal(wb.timeseries, wb.results.timeseries)
