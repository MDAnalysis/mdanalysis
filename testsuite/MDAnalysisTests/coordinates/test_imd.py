"""Test for MDAnalysis trajectory reader expectations
"""

import sys
import importlib
from weakref import ref
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
from MDAnalysis.coordinates.IMD import HAS_IMDCLIENT

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


@pytest.mark.skipif(not HAS_IMDCLIENT, reason="IMDClient not installed")
class IMDReference(BaseReference):
    def __init__(self):
        super(IMDReference, self).__init__()
        # Serve TRR traj data via the server
        traj = mda.coordinates.TRR.TRRReader(COORDINATES_TRR)
        self.server = InThreadIMDServer(traj)
        self.server.set_imdsessioninfo(create_default_imdsinfo_v3())

        self.n_atoms = traj.n_atoms
        self.prec = 3

        self.trajectory = "imd://localhost"
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
        reference = IMDReference()
        yield reference
        reference.server.cleanup()

    @staticmethod
    @pytest.fixture()
    def reader(ref):
        # This will start the test IMD Server, waiting for a connection
        # to then send handshake & first frame
        ref.server.handshake_sequence("localhost")
        # This will connect to the test IMD Server and read the first frame
        reader = ref.reader(
            ref.trajectory + str(ref.server.port), n_atoms=ref.n_atoms, buffer_size=1 * 1024 * 1024
        )
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
        yield reader
        reader.close()

    @staticmethod
    @pytest.fixture()
    def transformed(ref):
        # This will start the test IMD Server, waiting for a connection
        # to then send handshake & first frame
        ref.server.handshake_sequence("localhost")
        # This will connect to the test IMD Server and read the first frame
        transformed = ref.reader(
            ref.trajectory + str(ref.server.port), n_atoms=ref.n_atoms, buffer_size=1 * 1024 * 1024
        )
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
        ref_trans[:, 2] += 0.33
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


@pytest.fixture
def universe():
    return mda.Universe(COORDINATES_TOPOLOGY, COORDINATES_H5MD)


@pytest.fixture
def port():
    return get_free_port()


@pytest.fixture
def imdsinfo():
    return create_default_imdsinfo_v3()


@pytest.mark.skipif(not HAS_IMDCLIENT, reason="IMDClient not installed")
class TestStreamIteration:
    @pytest.fixture
    def reader(self, universe, imdsinfo):
        server = InThreadIMDServer(universe.trajectory)
        server.set_imdsessioninfo(imdsinfo)
        server.handshake_sequence("localhost", first_frame=True)
        reader = IMDReader(
            f"imd://localhost:{server.port}",
            n_atoms=universe.trajectory.n_atoms,
            buffer_size=1 * 1024 * 1024,
        )
        server.send_frames(1, 5)

        yield reader
        server.cleanup()
        reader.close()

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
def test_n_atoms_not_specified(universe, imdsinfo):
    server = InThreadIMDServer(universe.trajectory)
    server.set_imdsessioninfo(imdsinfo)
    server.handshake_sequence("localhost", first_frame=True)
    with pytest.raises(
        ValueError,
        match="IMDReader: n_atoms must be specified",
    ):
        IMDReader(
            f"imd://localhost:{server.port}",
        )
    server.cleanup()


@pytest.mark.skipif(not HAS_IMDCLIENT, reason="IMDClient not installed")
def test_imd_stream_empty(universe, imdsinfo):
    server = InThreadIMDServer(universe.trajectory)
    server.set_imdsessioninfo(imdsinfo)
    server.handshake_sequence("localhost", first_frame=False)
    with pytest.raises(
        RuntimeError,
        match="IMDReader: Read error",
    ):
        IMDReader(
            f"imd://localhost:{server.port}",
            n_atoms=universe.trajectory.n_atoms,
        )
    server.cleanup()


@pytest.mark.skipif(not HAS_IMDCLIENT, reason="IMDClient not installed")
def test_create_imd_universe(universe, imdsinfo):
    server = InThreadIMDServer(universe.trajectory)
    server.set_imdsessioninfo(imdsinfo)
    server.handshake_sequence("localhost", first_frame=True)
    u_imd = mda.Universe(
        COORDINATES_TOPOLOGY,
        f"imd://localhost:{server.port}",
        n_atoms=universe.trajectory.n_atoms,
    )
    assert type(u_imd.trajectory).__name__ == "IMDReader"
    with pytest.raises(ValueError, match="IMDReader: Invalid IMD URL"):
        u_imd = mda.Universe(
            COORDINATES_TOPOLOGY,
            f"imd://localhost:{port}/invalid",
            n_atoms=universe.trajectory.n_atoms,
        )
    server.cleanup()


def test_imd_format_hint():
    assert IMDReader._format_hint("imd://localhost:12345")
    assert IMDReader._format_hint("imd://localhost:12345/invalid")
    assert not IMDReader._format_hint("not_a_valid_imd_url")
    assert not IMDReader._format_hint(12345)
    assert not IMDReader._format_hint(None)
