"""Test for MDAnalysis trajectory reader expectations
"""

import sys
import pytest
import pickle
from types import ModuleType

import numpy as np
from numpy.testing import (
    assert_almost_equal,
    assert_array_almost_equal,
    assert_equal,
    assert_allclose,
)

from MDAnalysis.transformations import translate
import MDAnalysis as mda
from MDAnalysis.coordinates.IMD import HAS_IMDCLIENT, parse_host_port

if HAS_IMDCLIENT:
    import imdclient
    from imdclient.tests.utils import (
        get_free_port,
        create_default_imdsinfo_v3,
    )
    from imdclient.tests.server import InThreadIMDServer

from MDAnalysis.coordinates.IMD import IMDReader

from MDAnalysisTests.datafiles import (
    COORDINATES_TOPOLOGY,
    COORDINATES_TRR,
    COORDINATES_H5MD,
)

from MDAnalysisTests.coordinates.base import (
    MultiframeReaderTest,
    BaseReference,
    assert_timestep_almost_equal,
)


def test_HAS_IMDCLIENT_too_old():
    # mock a version of imdclient that is too old
    module_name = "imdclient"

    sys.modules.pop(module_name, None)
    sys.modules.pop("MDAnalysis.coordinates.IMD", None)

    mocked_module = ModuleType(module_name)
    # too old version
    mocked_module.__version__ = "0.1.0"
    sys.modules[module_name] = mocked_module

    from MDAnalysis.coordinates.IMD import HAS_IMDCLIENT

    assert not HAS_IMDCLIENT


def test_HAS_IMDCLIENT_new_enough():
    module_name = "imdclient"
    sys.modules.pop(module_name, None)
    sys.modules.pop("MDAnalysis.coordinates.IMD", None)

    mocked_module = ModuleType(module_name)
    IMDClient_module = ModuleType(f"{module_name}.IMDClient")

    class MockIMDClient:
        pass

    IMDClient_module.IMDClient = MockIMDClient
    mocked_module.IMDClient = IMDClient_module
    # new enough version
    mocked_module.__version__ = "0.1.4"

    sys.modules[module_name] = mocked_module
    sys.modules[f"{module_name}.IMDClient"] = IMDClient_module

    from MDAnalysis.coordinates.IMD import HAS_IMDCLIENT

    assert HAS_IMDCLIENT


@pytest.mark.skipif(not HAS_IMDCLIENT, reason="IMDClient not installed")
class IMDReference(BaseReference):
    def __init__(self):
        super(IMDReference, self).__init__()
        self.port = get_free_port()
        # Serve TRR traj data via the server
        traj = mda.coordinates.TRR.TRRReader(COORDINATES_TRR)
        self.server = InThreadIMDServer(traj)
        self.server.set_imdsessioninfo(create_default_imdsinfo_v3())

        self.n_atoms = traj.n_atoms
        self.prec = 3

        self.trajectory = f"imd://localhost:{self.port}"
        self.topology = COORDINATES_TOPOLOGY
        self.changing_dimensions = True
        self.reader = IMDReader

        self.first_frame.velocities = self.first_frame.positions / 10
        self.first_frame.forces = self.first_frame.positions / 100

        self.second_frame.velocities = self.second_frame.positions / 10
        self.second_frame.forces = self.second_frame.positions / 100

        self.last_frame.velocities = self.last_frame.positions / 10
        self.last_frame.forces = self.last_frame.positions / 100

        self.jump_to_frame.velocities = self.jump_to_frame.positions / 10
        self.jump_to_frame.forces = self.jump_to_frame.positions / 100

    def iter_ts(self, i):
        ts = self.first_frame.copy()
        ts.positions = 2**i * self.first_frame.positions
        ts.velocities = ts.positions / 10
        ts.forces = ts.positions / 100
        ts.time = i
        ts.frame = i
        return ts


@pytest.mark.skipif(not HAS_IMDCLIENT, reason="IMDClient not installed")
class TestIMDReaderBaseAPI(MultiframeReaderTest):

    @pytest.fixture(scope="function")
    def ref(self):
        """Not a static method like in base class- need new server for each test"""
        return IMDReference()

    @staticmethod
    @pytest.fixture()
    def reader(ref):
        # This will start the test IMD Server, waiting for a connection
        # to then send handshake & first frame
        ref.server.handshake_sequence("localhost", ref.port)
        # This will connect to the test IMD Server and read the first frame
        reader = ref.reader(ref.trajectory, n_atoms=ref.n_atoms)
        # Send the rest of the frames- small enough to all fit in socket itself
        ref.server.send_frames(1, 5)

        reader.add_auxiliary(
            "lowf",
            ref.aux_lowf,
            dt=ref.aux_lowf_dt,
            initial_time=0,
            time_selector=None,
        )
        reader.add_auxiliary(
            "highf",
            ref.aux_highf,
            dt=ref.aux_highf_dt,
            initial_time=0,
            time_selector=None,
        )
        return reader

    @staticmethod
    @pytest.fixture()
    def transformed(ref):
        # This will start the test IMD Server, waiting for a connection
        # to then send handshake & first frame
        ref.server.handshake_sequence("localhost", ref.port)
        # This will connect to the test IMD Server and read the first frame
        transformed = ref.reader(ref.trajectory, n_atoms=ref.n_atoms)
        # Send the rest of the frames- small enough to all fit in socket itself
        ref.server.send_frames(1, 5)
        transformed.add_transformations(
            translate([1, 1, 1]), translate([0, 0, 0.33])
        )
        return transformed

    def test_n_frames(self, ref, reader):
        pytest.skip("`n_frames` is unknown for IMDReader")

    def test_first_frame(self, ref, reader):
        # don't rewind here as in inherited base test
        assert_timestep_almost_equal(
            reader.ts, ref.first_frame, decimal=ref.prec
        )

    def test_get_writer_1(self, ref, reader, tmpdir):
        pytest.skip("No Writer for IMDReader")

    def test_get_writer_2(self, ref, reader, tmpdir):
        pytest.skip("No Writer for IMDReader")

    def test_total_time(self, ref, reader):
        pytest.skip("`total_time` is unknown for IMDReader")

    def test_changing_dimensions(self, ref, reader):
        pytest.skip("IMDReader cannot be rewound")

    def test_iter(self, ref, reader):
        for i, ts in enumerate(reader):
            assert_timestep_almost_equal(ts, ref.iter_ts(i), decimal=ref.prec)

    def test_first_dimensions(self, ref, reader):
        # don't rewind here as in inherited base test
        if ref.dimensions is None:
            assert reader.ts.dimensions is None
        else:
            assert_allclose(
                reader.ts.dimensions,
                ref.dimensions,
                rtol=0,
                atol=1.5 * 10 ** (-ref.prec),
            )

    def test_transformed(self, ref, transformed):
        # see transformed fixture
        ref_trans = ref.first_frame.positions + 1
        ref_trans[:, 2]  += 0.33
        assert_allclose(transformed.ts.positions, ref_trans)

    def test_volume(self, ref, reader):
        # don't rewind here as in inherited base test
        vol = reader.ts.volume
        # Here we can only be sure about the numbers upto the decimal point due
        # to floating point impressions.
        assert_allclose(vol, ref.volume, rtol=0, atol=1.5e0)

    def test_reload_auxiliaries_from_description(self, ref, reader):
        pytest.skip("Cannot create two IMDReaders on the same stream")

    def test_stop_iter(self, reader):
        pytest.skip("IMDReader cannot be rewound")

    def test_iter_rewinds(self, reader):
        pytest.skip("IMDReader cannot be rewound")

    def test_timeseries_shape(self, reader):
        pytest.skip("IMDReader does not support timeseries")

    def test_timeseries_asel_shape(self, reader):
        pytest.skip("IMDReader does not support timeseries")

    def test_timeseries_values(self, reader):
        pytest.skip("IMDReader does not support timeseries")

    def test_transformations_2iter(self, ref, transformed):
        pytest.skip("IMDReader cannot be reopened")

    def test_transformations_slice(self, ref, transformed):
        pytest.skip("IMDReader cannot be reopened")

    def test_transformations_switch_frame(self, ref, transformed):
        pytest.skip("IMDReader cannot be reopened")

    def test_transformation_rewind(self, ref, transformed):
        pytest.skip("IMDReader cannot be reopened")

    def test_pickle_reader(self, reader):
        pytest.skip("IMDReader cannot be pickled")

    def test_pickle_next_ts_reader(self, reader):
        pytest.skip("IMDReader cannot be pickled")

    def test_pickle_last_ts_reader(self, reader):
        pytest.skip("IMDReader cannot be pickled")

    def test_transformations_copy(self, ref, transformed):
        pytest.skip("IMDReader cannot be copied")

    def test_timeseries_empty_asel(self, reader):
        pytest.skip("IMDReader does not support timeseries")

    def test_timeseries_empty_atomgroup(self, reader):
        pytest.skip("IMDReader does not support timeseries")

    def test_timeseries_asel_warns_deprecation(self, reader):
        pytest.skip("IMDReader does not support timeseries")

    def test_timeseries_atomgroup(self, reader):
        pytest.skip("IMDReader does not support timeseries")

    def test_timeseries_atomgroup_asel_mutex(self, reader):
        pytest.skip("IMDReader does not support timeseries")

    def test_last_frame(self, ref, reader):
        pytest.skip("IMDReader cannot be rewound")

    def test_go_over_last_frame(self, ref, reader):
        pytest.skip("IMDReader must be an indexed using a slice")

    def test_frame_jump(self, ref, reader):
        pytest.skip("IMDReader must be an indexed using a slice")

    def test_frame_jump_issue1942(self, ref, reader):
        pytest.skip("IMDReader must be an indexed using a slice")

    def test_next_gives_second_frame(self, ref, reader):
        # don't recreate reader here as in inherited base test
        ts = reader.next()
        assert_timestep_almost_equal(ts, ref.second_frame, decimal=ref.prec)

    def test_frame_collect_all_same(self, reader):
        pytest.skip("IMDReader has independent coordinates")


@pytest.mark.skipif(not HAS_IMDCLIENT, reason="IMDClient not installed")
class TestStreamIteration:

    @pytest.fixture
    def port(self):
        return get_free_port()

    @pytest.fixture
    def universe(self):
        return mda.Universe(COORDINATES_TOPOLOGY, COORDINATES_H5MD)

    @pytest.fixture
    def imdsinfo(self):
        return create_default_imdsinfo_v3()

    @pytest.fixture
    def reader(self, universe, imdsinfo, port):
        server = InThreadIMDServer(universe.trajectory)
        server.set_imdsessioninfo(imdsinfo)
        server.handshake_sequence("localhost", port, first_frame=True)
        reader = IMDReader(
            f"imd://localhost:{port}",
            n_atoms=universe.trajectory.n_atoms,
        )
        server.send_frames(1, 5)

        yield reader
        server.cleanup()

    def test_iterate_step(self, reader, universe):
        i = 0
        for ts in reader[::2]:
            assert ts.frame == i
            i += 2

    def test_iterate_twice_sliced_raises_error(self, reader):
        for ts in reader[::2]:
            pass
        with pytest.raises(RuntimeError, match="Cannot reopen stream"):
            for ts in reader[::2]:
                pass

    def test_iterate_twice_all_raises_error(self, reader):
        for ts in reader:
            pass
        with pytest.raises(RuntimeError, match="Cannot reopen stream"):
            for ts in reader:
                pass

    def test_iterate_twice_fi_all_raises_error(self, reader):
        for ts in reader[:]:
            pass
        with pytest.raises(RuntimeError, match="Cannot reopen stream"):
            for ts in reader[:]:
                pass

    def test_index_stream_raises_error(self, reader):
        with pytest.raises(TypeError, match="Streamed trajectories must be"):
            reader[0]

    def test_iterate_backwards_raises_error(self, reader):
        with pytest.raises(ValueError, match="Cannot go backwards"):
            for ts in reader[::-1]:
                pass

    def test_iterate_start_stop_raises_error(self, reader):
        with pytest.raises(ValueError, match="Cannot expect a start index"):
            for ts in reader[1:3]:
                pass

    def test_subslice_fi_all_after_iteration_raises_error(self, reader):
        sliced_reader = reader[:]
        for ts in sliced_reader:
            pass
        sub_sliced_reader = sliced_reader[::1]
        with pytest.raises(RuntimeError):
            for ts in sub_sliced_reader:
                pass

    def test_timeseries_raises(self, reader):
        with pytest.raises(
            RuntimeError,
            match="cannot access timeseries for streamed trajectories",
        ):
            reader.timeseries()


@pytest.mark.skipif(not HAS_IMDCLIENT, reason="IMDClient not installed")
def test_n_atoms_mismatch():
    universe = mda.Universe(COORDINATES_TOPOLOGY, COORDINATES_H5MD)
    port = get_free_port()
    server = InThreadIMDServer(universe.trajectory)
    server.set_imdsessioninfo(create_default_imdsinfo_v3())
    server.handshake_sequence("localhost", port, first_frame=True)
    with pytest.raises(
        EOFError,
        match="IMDProducer: Expected n_atoms value 6, got 5. Ensure you are using the correct topology file.",
    ):
        IMDReader(
            f"imd://localhost:{port}",
            n_atoms=universe.trajectory.n_atoms + 1,
        )



@pytest.mark.skipif(not HAS_IMDCLIENT, reason="IMDClient not installed")
def test_n_atoms_not_specified():
    universe = mda.Universe(COORDINATES_TOPOLOGY, COORDINATES_H5MD)
    port = get_free_port()
    server = InThreadIMDServer(universe.trajectory)
    server.set_imdsessioninfo(create_default_imdsinfo_v3())
    server.handshake_sequence("localhost", port, first_frame=True)
    with pytest.raises(
        ValueError,
        match="IMDReader: n_atoms must be specified",
    ):
        IMDReader(
            f"imd://localhost:{port}",

        )


def test_parse_host_port():
    # Test with a valid host and port
    host, port = parse_host_port("imd://localhost:8889")
    assert host == "localhost"
    assert port == 8889

    # Test with a valid host and invalid port
    with pytest.raises(ValueError, match="IMDReader: Port must be an integer"):
        host, port = parse_host_port("imd://localhost:abcd")


    with pytest.raises(ValueError, match="IMDReader: URL must be in the format 'imd://host:port'"):
        host, port = parse_host_port("imd://localhost:blah:bleh")

