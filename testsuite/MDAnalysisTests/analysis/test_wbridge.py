import numpy as np
from collections import defaultdict

from numpy.testing import (
    assert_equal,
    assert_array_equal,
)
import pytest

from MDAnalysis.lib.util import NamedStream

import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.wbridge_analysis import (
    WaterBridgeAnalysis,
)


# Note: Currently the datafiles are added as a go-around for the
# incompatibility with StringIO with a fix in the future. Issue #5221
from MDAnalysisTests.datafiles import (
    WB_AD,
    WB_AWA,
    WB_AWA_AWWA,
    WB_AWD,
    WB_AWWA,
    WB_AWWWA,
    WB_AWWWWA,
    WB_BRANCH,
    WB_DA,
    WB_DA_PBC,
    WB_DWA,
    WB_DWD,
    WB_EMPTY,
    WB_LOOP,
    WB_DUPLICATE_WATER,
    WB_MULTIFRAME_GRO,
    WB_MULTIFRAME_DCD,
)


class TestWaterBridgeAnalysis(object):
    @staticmethod
    @pytest.fixture(scope="class")
    def wb_multiframe():
        """A water bridge object with multipley frames"""
        u = MDAnalysis.Universe(WB_MULTIFRAME_GRO, WB_MULTIFRAME_DCD)
        wb = WaterBridgeAnalysis(
            u, "protein and (resid 1)", "protein and (resid 4)", order=4
        )
        # Build an dummy WaterBridgeAnalysis object for testing
        wb.results.network = []
        wb.results.network.append({(1, 0, 12, None, 2.0, 180.0): None})
        wb.results.network.append({(0, None, 12, 13, 2.0, 180.0): None})
        wb.results.network.append(
            {(1, 0, 3, None, 2.0, 180.0): {(4, 2, 12, None, 2.0, 180.0): None}}
        )
        wb.results.network.append(
            {
                (0, None, 3, 2, 2.0, 180.0): {
                    (4, 2, 5, None, 2.0, 180.0): {
                        (5, None, 11, 12, 2.0, 180.0): None
                    }
                }
            }
        )
        wb.times = np.arange(len(wb.results.network))
        return wb

    def test_nodata(self):
        """Test if the funtions can run when there is no data.
        This is achieved by not runing the run() first."""
        universe_DA = MDAnalysis.Universe(WB_DA)
        wb = WaterBridgeAnalysis(
            universe_DA,
            "protein and (resid 1)",
            "protein and (resid 4)",
            order=0,
        )
        wb.generate_table()
        assert_equal(wb.timesteps_by_type(), None)
        assert_equal(wb.count_by_time(), None)
        assert_equal(wb.count_by_type(), None)

    def test_selection_type_error(self):
        """Test the case when the wrong selection1_type is given"""
        universe_DA = MDAnalysis.Universe(WB_DA)
        try:
            wb = WaterBridgeAnalysis(
                universe_DA,
                "protein and (resid 1)",
                "protein and (resid 4)",
                order=0,
                selection1_type="aaa",
            )
        except ValueError:
            pass
        else:
            raise pytest.fail("selection_type aaa should rasie error")

    def test_distance_type_error(self):
        """Test the case when the wrong selection1_type is given"""
        universe_DA = MDAnalysis.Universe(WB_DA)
        with pytest.raises(
            ValueError,
            match="Only 'hydrogen' and 'heavy' are allowed for option `distance_type'",
        ):
            WaterBridgeAnalysis(
                universe_DA,
                "protein and (resid 1)",
                "protein and (resid 4)",
                order=0,
                selection1_type="aaa",
                distance_type="aaa",
            )

    def test_selection2_type_error(self):
        """Test the case when the wrong selection1_type is given"""
        universe_DA = MDAnalysis.Universe(WB_DA)
        with pytest.raises(
            ValueError, match="`selection2_type` is not a keyword argument."
        ):
            WaterBridgeAnalysis(
                universe_DA,
                "protein and (resid 1)",
                "protein and (resid 4)",
                order=0,
                selection1_type="aaa",
                selection2_type="aaa",
            )

    def test_empty_selection(self, client_WaterBridgeAnalysis):
        """Test the case when selection yields empty result"""
        universe_DA = MDAnalysis.Universe(WB_DA)
        wb = WaterBridgeAnalysis(
            universe_DA,
            "protein and (resid 9)",
            "protein and (resid 10)",
            order=0,
        )
        wb.run(**client_WaterBridgeAnalysis)
        assert wb.results.network == [{}]

    def test_loop(self, client_WaterBridgeAnalysis):
        """Test if loop can be handled correctly"""
        universe_loop = MDAnalysis.Universe(WB_LOOP)
        wb = WaterBridgeAnalysis(
            universe_loop,
            "protein and (resid 1)",
            "protein and (resid 1 or resid 4)",
        )
        wb.run(**client_WaterBridgeAnalysis)
        assert_equal(len(wb.results.network[0].keys()), 2)

    @pytest.mark.parametrize("distance_type", ["hydrogen", "heavy"])
    def test_donor_accepter(self, distance_type, client_WaterBridgeAnalysis):
        """Test zeroth order donor to acceptor hydrogen bonding"""
        universe_DA = MDAnalysis.Universe(WB_DA)
        wb = WaterBridgeAnalysis(
            universe_DA,
            "protein and (resid 1)",
            "protein and (resid 4)",
            order=0,
            update_selection=True,
            debug=True,
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (1, 0, 2, None))

    @pytest.mark.parametrize("distance_type", ["hydrogen", "heavy"])
    def test_donor_accepter_pbc(
        self, distance_type, client_WaterBridgeAnalysis
    ):
        """Test zeroth order donor to acceptor hydrogen bonding in PBC conditions"""
        universe_DA_PBC = MDAnalysis.Universe(WB_DA_PBC)
        wb = WaterBridgeAnalysis(
            universe_DA_PBC,
            "protein and (resid 1)",
            "protein and (resid 4)",
            order=0,
            pbc=True,
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (1, 0, 2, None))

    @pytest.mark.parametrize("distance_type", ["hydrogen", "heavy"])
    def test_accepter_donor(self, distance_type, client_WaterBridgeAnalysis):
        """Test zeroth order acceptor to donor hydrogen bonding"""
        universe_AD = MDAnalysis.Universe(WB_AD)
        wb = WaterBridgeAnalysis(
            universe_AD,
            "protein and (resid 1)",
            "protein and (resid 4)",
            order=0,
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 1, 2))

    @pytest.mark.parametrize("distance_type", ["hydrogen", "heavy"])
    def test_acceptor_water_accepter(
        self, distance_type, client_WaterBridgeAnalysis
    ):
        """Test case where the hydrogen bond acceptor from selection 1 form
        water bridge with hydrogen bond acceptor from selection 2"""
        universe_AWA = MDAnalysis.Universe(WB_AWA)
        wb = WaterBridgeAnalysis(
            universe_AWA,
            "protein and (resid 1)",
            "protein and (resid 4)",
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        assert_equal(second[list(second.keys())[0]], None)

    @pytest.mark.parametrize("distance_type", ["hydrogen", "heavy"])
    def test_donor_water_accepter(
        self, distance_type, client_WaterBridgeAnalysis
    ):
        """Test case where the hydrogen bond donor from selection 1 form
        water bridge with hydrogen bond acceptor from selection 2"""
        universe_DWA = MDAnalysis.Universe(WB_DWA)
        wb = WaterBridgeAnalysis(
            universe_DWA,
            "protein and (resid 1)",
            "protein and (resid 4)",
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (1, 0, 2, None))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 2, 4, None))
        assert_equal(second[list(second.keys())[0]], None)

    @pytest.mark.parametrize("distance_type", ["hydrogen", "heavy"])
    def test_acceptor_water_donor(
        self, distance_type, client_WaterBridgeAnalysis
    ):
        """Test case where the hydrogen bond acceptor from selection 1 form
        water bridge with hydrogen bond donor from selection 2"""
        universe_AWD = MDAnalysis.Universe(WB_AWD)
        wb = WaterBridgeAnalysis(
            universe_AWD,
            "protein and (resid 1)",
            "protein and (resid 4)",
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (1, None, 3, 4))
        assert_equal(second[list(second.keys())[0]], None)

    @pytest.mark.parametrize("distance_type", ["hydrogen", "heavy"])
    def test_donor_water_donor(
        self, distance_type, client_WaterBridgeAnalysis
    ):
        """Test case where the hydrogen bond donor from selection 1 form
        water bridge with hydrogen bond donor from selection 2"""
        universe_DWD = MDAnalysis.Universe(WB_DWD)
        wb = WaterBridgeAnalysis(
            universe_DWD,
            "protein and (resid 1)",
            "protein and (resid 4)",
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (1, 0, 2, None))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (2, None, 3, 4))
        assert_equal(second[list(second.keys())[0]], None)

    def test_empty(self, client_WaterBridgeAnalysis):
        """Test case where no water bridge exists"""
        universe_empty = MDAnalysis.Universe(WB_EMPTY)
        wb = WaterBridgeAnalysis(universe_empty, "protein", "protein")
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        assert_equal(wb.results.network[0], defaultdict(dict))

    def test_same_selection(self, client_WaterBridgeAnalysis):
        """
        This test tests that if the selection 1 and selection 2 are both protein.
        However, the protein only forms one hydrogen bond with the water.
        This entry won't be included.
        """
        universe_DWA = MDAnalysis.Universe(WB_DWA)
        wb = WaterBridgeAnalysis(
            universe_DWA, "protein and resid 1", "protein and resid 1"
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        assert_equal(wb.results.network[0], defaultdict(dict))

    @pytest.mark.parametrize("distance_type", ["hydrogen", "heavy"])
    def test_acceptor_2water_accepter(
        self, distance_type, client_WaterBridgeAnalysis
    ):
        """Test case where the hydrogen bond acceptor from selection 1 form second order
        water bridge with hydrogen bond acceptor from selection 2"""
        # test first order
        universe_AWWA = MDAnalysis.Universe(WB_AWWA)
        wb = WaterBridgeAnalysis(
            universe_AWWA,
            "protein and (resid 1)",
            "protein and (resid 4)",
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        assert_equal(wb.results.network[0], defaultdict(dict))
        # test second order
        wb = WaterBridgeAnalysis(
            universe_AWWA,
            "protein and (resid 1)",
            "protein and (resid 4)",
            order=2,
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 4, 6, None))
        assert_equal(third[list(third.keys())[0]], None)
        # test third order
        wb = WaterBridgeAnalysis(
            universe_AWWA,
            "protein and (resid 1)",
            "protein and (resid 4)",
            order=3,
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 4, 6, None))
        assert_equal(third[list(third.keys())[0]], None)

    @pytest.mark.parametrize("distance_type", ["hydrogen", "heavy"])
    def test_acceptor_3water_accepter(
        self, distance_type, client_WaterBridgeAnalysis
    ):
        """Test case where the hydrogen bond acceptor from selection 1 form third order
        water bridge with hydrogen bond acceptor from selection 2"""
        universe_AWWWA = MDAnalysis.Universe(WB_AWWWA)
        wb = WaterBridgeAnalysis(
            universe_AWWWA,
            "protein and (resid 1)",
            "protein and (resid 5)",
            order=2,
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        assert_equal(wb.results.network[0], defaultdict(dict))

        wb = WaterBridgeAnalysis(
            universe_AWWWA,
            "protein and (resid 1)",
            "protein and (resid 5)",
            order=3,
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 4, 6, None))
        fourth = third[list(third.keys())[0]]
        assert_equal(list(fourth.keys())[0][:4], (7, 6, 8, None))
        assert_equal(fourth[list(fourth.keys())[0]], None)

        wb = WaterBridgeAnalysis(
            universe_AWWWA,
            "protein and (resid 1)",
            "protein and (resid 5)",
            order=4,
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        third = second[list(second.keys())[0]]
        assert_equal(list(third.keys())[0][:4], (5, 4, 6, None))
        fourth = third[list(third.keys())[0]]
        assert_equal(list(fourth.keys())[0][:4], (7, 6, 8, None))
        assert_equal(fourth[list(fourth.keys())[0]], None)

    @pytest.mark.parametrize("distance_type", ["hydrogen", "heavy"])
    def test_acceptor_4water_accepter(
        self, distance_type, client_WaterBridgeAnalysis
    ):
        """Test case where the hydrogen bond acceptor from selection 1 form fourth order
        water bridge with hydrogen bond acceptor from selection 2"""
        universe_AWWWWA = MDAnalysis.Universe(WB_AWWWWA)
        wb = WaterBridgeAnalysis(
            universe_AWWWWA,
            "protein and (resid 1)",
            "protein and (resid 6)",
            order=3,
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        assert_equal(wb.results.network[0], defaultdict(dict))

        wb = WaterBridgeAnalysis(
            universe_AWWWWA,
            "protein and (resid 1)",
            "protein and (resid 6)",
            order=4,
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
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

        wb = WaterBridgeAnalysis(
            universe_AWWWWA,
            "protein and (resid 1)",
            "protein and (resid 6)",
            order=5,
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
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

    @pytest.mark.parametrize("distance_type", ["hydrogen", "heavy"])
    def test_acceptor_22water_accepter(
        self, distance_type, client_WaterBridgeAnalysis
    ):
        """Test case where the hydrogen bond acceptor from selection 1 form a second order
        water bridge with hydrogen bond acceptor from selection 2
        and the last water is linked to two residues in selection 2"""
        universe_branch = MDAnalysis.Universe(WB_BRANCH)
        wb = WaterBridgeAnalysis(
            universe_branch,
            "protein and (resid 1)",
            "protein and (resid 4 or resid 5)",
            order=2,
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        third = second[list(second.keys())[0]]
        assert_equal(
            [(5, 4, 7, None), (6, 4, 8, None)],
            sorted([key[:4] for key in list(third.keys())]),
        )

    def test_timeseries_wba(self, client_WaterBridgeAnalysis):
        """Test if the time series data is correctly generated in water bridge analysis format"""
        universe_branch = MDAnalysis.Universe(WB_BRANCH)
        wb = WaterBridgeAnalysis(
            universe_branch,
            "protein and (resid 1)",
            "protein and (resid 4 or resid 5)",
            order=2,
        )
        wb.output_format = "sele1_sele2"
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        timeseries = sorted(wb.results.timeseries[0])

        assert_equal(
            timeseries[0][:4], (0, 2, ("ALA", 1, "O"), ("SOL", 2, "HW1"))
        )
        assert_equal(
            timeseries[1][:4], (3, 4, ("SOL", 2, "HW2"), ("SOL", 3, "OW"))
        )
        assert_equal(
            timeseries[2][:4], (5, 7, ("SOL", 3, "HW1"), ("ALA", 4, "O"))
        )
        assert_equal(
            timeseries[3][:4], (6, 8, ("SOL", 3, "HW2"), ("ALA", 5, "O"))
        )

    def test_timeseries_hba(self, client_WaterBridgeAnalysis):
        """Test if the time series data is correctly generated in hydrogen bond analysis format"""
        universe_branch = MDAnalysis.Universe(WB_BRANCH)
        wb = WaterBridgeAnalysis(
            universe_branch,
            "protein and (resid 1)",
            "protein and (resid 4 or resid 5)",
            order=2,
        )
        wb.output_format = "donor_acceptor"
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        timeseries = sorted(wb.results.timeseries[0])

        assert_equal(
            timeseries[0][:4], (2, 0, ("SOL", 2, "HW1"), ("ALA", 1, "O"))
        )
        assert_equal(
            timeseries[1][:4], (3, 4, ("SOL", 2, "HW2"), ("SOL", 3, "OW"))
        )
        assert_equal(
            timeseries[2][:4], (5, 7, ("SOL", 3, "HW1"), ("ALA", 4, "O"))
        )
        assert_equal(
            timeseries[3][:4], (6, 8, ("SOL", 3, "HW2"), ("ALA", 5, "O"))
        )

    @pytest.mark.parametrize("distance_type", ["hydrogen", "heavy"])
    def test_acceptor_12water_accepter(
        self, distance_type, client_WaterBridgeAnalysis
    ):
        """Test of independent first order and second can be recognised correctely"""
        universe_AWA_AWWA = MDAnalysis.Universe(WB_AWA_AWWA)
        wb = WaterBridgeAnalysis(
            universe_AWA_AWWA,
            "protein and (resid 1 or resid 5)",
            "protein and (resid 4 or resid 8)",
            order=1,
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        network = wb.results.network[0]
        assert_equal(list(network.keys())[0][:4], (0, None, 2, 1))
        second = network[list(network.keys())[0]]
        assert_equal(list(second.keys())[0][:4], (3, 1, 4, None))
        assert_equal(second[list(second.keys())[0]], None)
        network = wb.results.network[0]
        wb = WaterBridgeAnalysis(
            universe_AWA_AWWA,
            "protein and (resid 1 or resid 5)",
            "protein and (resid 4 or resid 8)",
            order=2,
            distance_type=distance_type,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        network = wb.results.network[0]
        assert_equal(
            [(0, None, 2, 1), (5, None, 7, 6)],
            sorted([key[:4] for key in list(network.keys())]),
        )

    def test_count_by_type_single_link(self, client_WaterBridgeAnalysis):
        """
        This test tests the simplest water bridge to see if count_by_type() works.
        """
        universe_DWA = MDAnalysis.Universe(WB_DWA)
        wb = WaterBridgeAnalysis(
            universe_DWA, "protein and (resid 1)", "protein and (resid 4)"
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        assert_equal(
            wb.count_by_type(), [(1, 4, "ALA", 1, "H", "ALA", 4, "O", 1.0)]
        )

    def test_count_by_type_multiple_link(self, client_WaterBridgeAnalysis):
        """
        This test tests if count_by_type() can give the correct result for more than 1 links.
        """
        universe_AWA_AWWA = MDAnalysis.Universe(WB_AWA_AWWA)
        wb = WaterBridgeAnalysis(
            universe_AWA_AWWA,
            "protein and (resid 1 or resid 5)",
            "protein and (resid 4 or resid 8)",
            order=2,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)
        assert_equal(
            sorted(wb.count_by_type()),
            [
                [0, 4, "ALA", 1, "O", "ALA", 4, "O", 1.0],
                [5, 11, "ALA", 5, "O", "ALA", 8, "O", 1.0],
            ],
        )

    def test_count_by_type_multiple_frame(self, wb_multiframe):
        """
        This test tests if count_by_type() works in multiply situations.
        :return:
        """
        result = [
            [0, 11, "ALA", 1, "O", "ALA", 6, "H", 0.25],
            [0, 12, "ALA", 1, "O", "ALA", 6, "O", 0.25],
            [1, 12, "ALA", 1, "H", "ALA", 6, "O", 0.5],
        ]

        assert_equal(sorted(wb_multiframe.count_by_type()), result)

    def test_count_by_type_filter(self, wb_multiframe):
        """
        This test tests if modifying analysis_func
        allows some results to be filtered out in count_by_type().
        :return:
        """

        def analysis(current, output, u):
            sele1_index, sele1_heavy_index, atom2, heavy_atom2, dist, angle = (
                current[0]
            )
            atom1, heavy_atom1, sele2_index, sele2_heavy_index, dist, angle = (
                current[-1]
            )
            sele1 = u.atoms[sele1_index]
            sele2 = u.atoms[sele2_index]
            (s1_resname, s1_resid, s1_name) = (
                sele1.resname,
                sele1.resid,
                sele1.name,
            )
            (s2_resname, s2_resid, s2_name) = (
                sele2.resname,
                sele2.resid,
                sele2.name,
            )

            key = (
                sele1_index,
                sele2_index,
                s1_resname,
                s1_resid,
                s1_name,
                s2_resname,
                s2_resid,
                s2_name,
            )
            if s2_name == "H":
                output[key] += 1

        result = [((0, 11, "ALA", 1, "O", "ALA", 6, "H"), 0.25)]
        assert_equal(
            sorted(wb_multiframe.count_by_type(analysis_func=analysis)), result
        )

    def test_count_by_type_merge(self, wb_multiframe):
        """
        This test tests if modifying analysis_func
        allows some same residue to be merged in count_by_type().
        """

        def analysis(current, output, u):
            sele1_index, sele1_heavy_index, atom2, heavy_atom2, dist, angle = (
                current[0]
            )
            atom1, heavy_atom1, sele2_index, sele2_heavy_index, dist, angle = (
                current[-1]
            )
            sele1 = u.atoms[sele1_index]
            sele2 = u.atoms[sele2_index]
            (s1_resname, s1_resid, s1_name) = (
                sele1.resname,
                sele1.resid,
                sele1.name,
            )
            (s2_resname, s2_resid, s2_name) = (
                sele2.resname,
                sele2.resid,
                sele2.name,
            )

            key = (s1_resname, s1_resid, s2_resname, s2_resid)
            output[key] = 1

        result = [(("ALA", 1, "ALA", 6), 1.0)]
        assert_equal(
            sorted(wb_multiframe.count_by_type(analysis_func=analysis)), result
        )

    def test_count_by_type_order(self, wb_multiframe):
        """
        This test tests if modifying analysis_func
        allows the order of water bridge to be separated in count_by_type().
        :return:
        """

        def analysis(current, output, u):
            sele1_index, sele1_heavy_index, atom2, heavy_atom2, dist, angle = (
                current[0]
            )
            atom1, heavy_atom1, sele2_index, sele2_heavy_index, dist, angle = (
                current[-1]
            )
            sele1 = u.atoms[sele1_index]
            sele2 = u.atoms[sele2_index]
            (s1_resname, s1_resid, s1_name) = (
                sele1.resname,
                sele1.resid,
                sele1.name,
            )
            (s2_resname, s2_resid, s2_name) = (
                sele2.resname,
                sele2.resid,
                sele2.name,
            )
            key = (
                s1_resname,
                s1_resid,
                s2_resname,
                s2_resid,
                len(current) - 1,
            )
            output[key] = 1

        result = [
            (("ALA", 1, "ALA", 6, 0), 0.5),
            (("ALA", 1, "ALA", 6, 1), 0.25),
            (("ALA", 1, "ALA", 6, 2), 0.25),
        ]
        assert_equal(
            sorted(wb_multiframe.count_by_type(analysis_func=analysis)), result
        )

    def test_count_by_time(self, wb_multiframe):
        """
        This test tests if count_by_times() works.
        :return:
        """
        assert_equal(
            wb_multiframe.count_by_time(), [(0, 1), (1, 1), (2, 1), (3, 1)]
        )

    def test_count_by_time_weight(self, client_WaterBridgeAnalysis):
        """
        This test tests if modyfing the analysis_func allows the weight to be changed
        in count_by_type().
        :return:
        """
        universe_AWA_AWWA = MDAnalysis.Universe(WB_AWA_AWWA)
        wb = WaterBridgeAnalysis(
            universe_AWA_AWWA,
            "protein and (resid 1 or resid 5)",
            "protein and (resid 4 or resid 8)",
            order=2,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)

        def analysis(current, output, u):
            sele1_index, sele1_heavy_index, atom2, heavy_atom2, dist, angle = (
                current[0]
            )
            atom1, heavy_atom1, sele2_index, sele2_heavy_index, dist, angle = (
                current[-1]
            )
            sele1 = u.atoms[sele1_index]
            sele2 = u.atoms[sele2_index]
            (s1_resname, s1_resid, s1_name) = (
                sele1.resname,
                sele1.resid,
                sele1.name,
            )
            (s2_resname, s2_resid, s2_name) = (
                sele2.resname,
                sele2.resid,
                sele2.name,
            )
            key = (s1_resname, s1_resid, s2_resname, s2_resid)
            output[key] += len(current) - 1

        assert_equal(
            wb.count_by_time(analysis_func=analysis),
            [
                (0, 3),
            ],
        )

    def test_count_by_time_empty(self, client_WaterBridgeAnalysis):
        """
        See if count_by_time() can handle zero well.
        :return:
        """
        universe_AWA_AWWA = MDAnalysis.Universe(WB_AWA_AWWA)
        wb = WaterBridgeAnalysis(
            universe_AWA_AWWA,
            "protein and (resid 1 or resid 5)",
            "protein and (resid 4 or resid 8)",
            order=2,
        )
        wb.run(**client_WaterBridgeAnalysis, verbose=False)

        def analysis(current, output, u):
            pass

        assert_equal(
            wb.count_by_time(analysis_func=analysis),
            [
                (0, 0),
            ],
        )

    def test_generate_table_hba(self, wb_multiframe):
        """Test generate table using hydrogen bond analysis format"""
        table = wb_multiframe.generate_table(output_format="donor_acceptor")
        assert_array_equal(
            sorted(table.donor_resid),
            [1, 1, 2, 2, 2, 6, 6],
        )

    def test_generate_table_s1s2(self, wb_multiframe):
        """Test generate table using hydrogen bond analysis format"""
        table = wb_multiframe.generate_table(output_format="sele1_sele2")
        assert_array_equal(
            sorted(table.sele1_resid),
            [1, 1, 1, 1, 2, 2, 3],
        )

    def test_timesteps_by_type(self, wb_multiframe):
        """Test the timesteps_by_type function"""

        timesteps = sorted(wb_multiframe.timesteps_by_type())
        assert_array_equal(
            timesteps[3], [1, 12, "ALA", 1, "H", "ALA", 6, "O", 0, 2]
        )

    def test_duplicate_water(self, client_WaterBridgeAnalysis):
        u = MDAnalysis.Universe(WB_DUPLICATE_WATER)
        wb = WaterBridgeAnalysis(
            u, "resname LEU and name O", "resname LEU and name N H", order=4
        )
        wb.run(**client_WaterBridgeAnalysis)
        assert len(wb.results.timeseries[0]) == 2

    def test_warn_results_deprecated(self, client_WaterBridgeAnalysis):
        universe_DA = MDAnalysis.Universe(WB_DA)
        wb = WaterBridgeAnalysis(
            universe_DA,
            "protein and (resid 9)",
            "protein and (resid 10)",
            order=0,
        )
        wb.run(**client_WaterBridgeAnalysis)

        wmsg = "The `network` attribute was deprecated in MDAnalysis 2.0.0"
        with pytest.warns(DeprecationWarning, match=wmsg):
            assert_equal(wb.network, wb.results.network)

        wmsg = "The `timeseries` attribute was deprecated in MDAnalysis 2.0.0"
        with pytest.warns(DeprecationWarning, match=wmsg):
            assert_equal(wb.timeseries, wb.results.timeseries)
