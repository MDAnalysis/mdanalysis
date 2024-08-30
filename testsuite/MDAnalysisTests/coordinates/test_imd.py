from MDAnalysisTests.datafiles import (
    COORDINATES_TOPOLOGY,
    COORDINATES_TRR,
    COORDINATES_H5MD,
)
import MDAnalysis as mda
from imdclient.tests.utils import (
    get_free_port,
    create_default_imdsinfo_v3,
)

from imdclient.tests.server import InThreadIMDServer
from MDAnalysisTests.coordinates.base import (
    MultiframeReaderTest,
    BaseReference,
    BaseWriterTest,
    assert_timestep_almost_equal,
)
from numpy.testing import (
    assert_almost_equal,
    assert_array_almost_equal,
    assert_equal,
    assert_allclose,
)
import numpy as np
import logging
import pytest
from MDAnalysis.transformations import translate
import pickle


logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("test.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)


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

        self.trajectory = f"localhost:{self.port}"
        self.topology = COORDINATES_TOPOLOGY
        self.changing_dimensions = True
        self.reader = mda.coordinates.IMD.IMDReader

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


class TestIMDReaderBaseAPI(MultiframeReaderTest):

    @pytest.fixture()
    def ref(self):
        """Not a static method like in base class- need new server for each test"""
        return IMDReference()

    @pytest.fixture()
    def reader(self, ref):
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

    @pytest.fixture()
    def transformed(self, ref):
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

    @pytest.mark.skip(
        reason="Stream-based reader cannot determine n_frames until EOF"
    )
    def test_n_frames(self, reader, ref):
        assert_equal(
            self.universe.trajectory.n_frames,
            1,
            "wrong number of frames in pdb",
        )

    def test_first_frame(self, ref, reader):
        # don't rewind here as in inherited base test
        assert_timestep_almost_equal(
            reader.ts, ref.first_frame, decimal=ref.prec
        )

    @pytest.mark.skip(reason="IMD is not a writeable format")
    def test_get_writer_1(self, ref, reader, tmpdir):
        with tmpdir.as_cwd():
            outfile = "test-writer." + ref.ext
            with reader.Writer(outfile) as W:
                assert_equal(isinstance(W, ref.writer), True)
                assert_equal(W.n_atoms, reader.n_atoms)

    @pytest.mark.skip(reason="IMD is not a writeable format")
    def test_get_writer_2(self, ref, reader, tmpdir):
        with tmpdir.as_cwd():
            outfile = "test-writer." + ref.ext
            with reader.Writer(outfile, n_atoms=100) as W:
                assert_equal(isinstance(W, ref.writer), True)
                assert_equal(W.n_atoms, 100)

    @pytest.mark.skip(
        reason="Stream-based reader cannot determine total_time until EOF"
    )
    def test_total_time(self, reader, ref):
        assert_almost_equal(reader.totaltime, ref.totaltime, decimal=ref.prec)

    @pytest.mark.skip(reason="Stream-based reader can only be read iteratively")
    def test_changing_dimensions(self, ref, reader):
        if ref.changing_dimensions:
            reader.rewind()
            if ref.dimensions is None:
                assert reader.ts.dimensions is None
            else:
                assert_array_almost_equal(
                    reader.ts.dimensions, ref.dimensions, decimal=ref.prec
                )
            reader[1]
            if ref.dimensions_second_frame is None:
                assert reader.ts.dimensions is None
            else:
                assert_array_almost_equal(
                    reader.ts.dimensions,
                    ref.dimensions_second_frame,
                    decimal=ref.prec,
                )

    def test_iter(self, ref, reader):
        for i, ts in enumerate(reader):
            assert_timestep_almost_equal(ts, ref.iter_ts(i), decimal=ref.prec)

    def test_first_dimensions(self, ref, reader):
        # don't rewind here as in inherited base test
        if ref.dimensions is None:
            assert reader.ts.dimensions is None
        else:
            assert_array_almost_equal(
                reader.ts.dimensions, ref.dimensions, decimal=ref.prec
            )

    def test_volume(self, ref, reader):
        # don't rewind here as in inherited base test
        vol = reader.ts.volume
        # Here we can only be sure about the numbers upto the decimal point due
        # to floating point impressions.
        assert_almost_equal(vol, ref.volume, 0)

    @pytest.mark.skip(reason="Cannot create new reader from same stream")
    def test_reload_auxiliaries_from_description(self, ref, reader):
        # get auxiliary desscriptions form existing reader
        descriptions = reader.get_aux_descriptions()
        # load a new reader, without auxiliaries
        reader = ref.reader(ref.trajectory)
        # load auxiliaries into new reader, using description...
        for aux in descriptions:
            reader.add_auxiliary(**aux)
        # should have the same number of auxiliaries
        assert_equal(
            reader.aux_list,
            reader.aux_list,
            "Number of auxiliaries does not match",
        )
        # each auxiliary should be the same
        for auxname in reader.aux_list:
            assert_equal(
                reader._auxs[auxname],
                reader._auxs[auxname],
                "AuxReaders do not match",
            )

    @pytest.mark.skip(reason="Stream can only be read in for loop")
    def test_stop_iter(self, reader):
        # reset to 0
        reader.rewind()
        for ts in reader[:-1]:
            pass
        assert_equal(reader.frame, 0)

    @pytest.mark.skip(reason="Cannot rewind stream")
    def test_iter_rewinds(self, reader, accessor):
        for ts_indices in accessor(reader):
            pass
        assert_equal(ts_indices.frame, 0)

    @pytest.mark.skip(
        reason="Timeseries currently requires n_frames to be known"
    )
    @pytest.mark.parametrize(
        "order", ("fac", "fca", "afc", "acf", "caf", "cfa")
    )
    def test_timeseries_shape(self, reader, order):
        timeseries = reader.timeseries(order=order)
        a_index = order.index("a")
        # f_index = order.index("f")
        c_index = order.index("c")
        assert timeseries.shape[a_index] == reader.n_atoms
        # assert timeseries.shape[f_index] == len(reader)
        assert timeseries.shape[c_index] == 3

    @pytest.mark.skip(
        reason="Timeseries currently requires n_frames to be known"
    )
    @pytest.mark.parametrize("asel", ("index 1", "index 2", "index 1 to 3"))
    def test_timeseries_asel_shape(self, reader, asel):
        atoms = mda.Universe(reader.filename).select_atoms(asel)
        timeseries = reader.timeseries(atoms, order="fac")
        assert timeseries.shape[0] == len(reader)
        assert timeseries.shape[1] == len(atoms)
        assert timeseries.shape[2] == 3

    @pytest.mark.skip("Cannot slice stream")
    @pytest.mark.parametrize("slice", ([0, 2, 1], [0, 10, 2], [0, 10, 3]))
    def test_timeseries_values(self, reader, slice):
        ts_positions = []
        if isinstance(reader, mda.coordinates.memory.MemoryReader):
            pytest.xfail(
                "MemoryReader uses deprecated stop inclusive"
                " indexing, see Issue #3893"
            )
        if slice[1] > len(reader):
            pytest.skip("too few frames in reader")
        for i in range(slice[0], slice[1], slice[2]):
            ts = reader[i]
            ts_positions.append(ts.positions.copy())
        positions = np.asarray(ts_positions)
        timeseries = reader.timeseries(
            start=slice[0], stop=slice[1], step=slice[2], order="fac"
        )
        assert_allclose(timeseries, positions)

    @pytest.mark.skip(reason="Cannot rewind stream")
    def test_transformations_2iter(self, ref, transformed):
        # Are the transformations applied and
        # are the coordinates "overtransformed"?
        v1 = np.float32((1, 1, 1))
        v2 = np.float32((0, 0, 0.33))
        idealcoords = []
        for i, ts in enumerate(transformed):
            idealcoords.append(ref.iter_ts(i).positions + v1 + v2)
            assert_array_almost_equal(
                ts.positions, idealcoords[i], decimal=ref.prec
            )

        for i, ts in enumerate(transformed):
            assert_almost_equal(ts.positions, idealcoords[i], decimal=ref.prec)

    @pytest.mark.skip(reason="Cannot slice stream")
    def test_transformations_slice(self, ref, transformed):
        # Are the transformations applied when iterating over a slice of the trajectory?
        v1 = np.float32((1, 1, 1))
        v2 = np.float32((0, 0, 0.33))
        for i, ts in enumerate(transformed[2:3:1]):
            idealcoords = ref.iter_ts(ts.frame).positions + v1 + v2
            assert_array_almost_equal(
                ts.positions, idealcoords, decimal=ref.prec
            )

    @pytest.mark.skip(reason="Cannot slice stream")
    def test_transformations_switch_frame(self, ref, transformed):
        # This test checks if the transformations are applied and if the coordinates
        # "overtransformed" on different situations
        # Are the transformations applied when we switch to a different frame?
        v1 = np.float32((1, 1, 1))
        v2 = np.float32((0, 0, 0.33))
        first_ideal = ref.iter_ts(0).positions + v1 + v2
        if len(transformed) > 1:
            assert_array_almost_equal(
                transformed[0].positions, first_ideal, decimal=ref.prec
            )
            second_ideal = ref.iter_ts(1).positions + v1 + v2
            assert_array_almost_equal(
                transformed[1].positions, second_ideal, decimal=ref.prec
            )

            # What if we comeback to the previous frame?
            assert_array_almost_equal(
                transformed[0].positions, first_ideal, decimal=ref.prec
            )

            # How about we switch the frame to itself?
            assert_array_almost_equal(
                transformed[0].positions, first_ideal, decimal=ref.prec
            )
        else:
            assert_array_almost_equal(
                transformed[0].positions, first_ideal, decimal=ref.prec
            )

    @pytest.mark.skip(reason="Cannot rewind stream")
    def test_transformation_rewind(self, ref, transformed):
        # this test checks if the transformations are applied after rewinding the
        # trajectory
        v1 = np.float32((1, 1, 1))
        v2 = np.float32((0, 0, 0.33))
        ideal_coords = ref.iter_ts(0).positions + v1 + v2
        transformed.rewind()
        assert_array_almost_equal(
            transformed[0].positions, ideal_coords, decimal=ref.prec
        )

    @pytest.mark.skip(reason="Cannot make a copy of a stream")
    def test_copy(self, ref, transformed):
        # this test checks if transformations are carried over a copy and if the
        # coordinates of the copy are equal to the original's
        v1 = np.float32((1, 1, 1))
        v2 = np.float32((0, 0, 0.33))
        new = transformed.copy()
        assert_equal(
            transformed.transformations,
            new.transformations,
            "transformations are not equal",
        )
        for i, ts in enumerate(new):
            ideal_coords = ref.iter_ts(i).positions + v1 + v2
            assert_array_almost_equal(
                ts.positions, ideal_coords, decimal=ref.prec
            )

    @pytest.mark.skip(reason="Cannot pickle socket")
    def test_pickle_reader(self, reader):
        """It probably wouldn't be a good idea to pickle a
        reader that is connected to a server"""
        reader_p = pickle.loads(pickle.dumps(reader))
        assert_equal(len(reader), len(reader_p))
        assert_equal(
            reader.ts, reader_p.ts, "Timestep is changed after pickling"
        )

    @pytest.mark.skip(reason="Cannot pickle socket")
    def test_pickle_next_ts_reader(self, reader):
        reader_p = pickle.loads(pickle.dumps(reader))
        assert_equal(
            next(reader),
            next(reader_p),
            "Next timestep is changed after pickling",
        )

    @pytest.mark.skip(reason="Cannot pickle socket")
    def test_pickle_last_ts_reader(self, reader):
        #  move current ts to last frame.
        reader[-1]
        reader_p = pickle.loads(pickle.dumps(reader))
        assert_equal(
            len(reader),
            len(reader_p),
            "Last timestep is changed after pickling",
        )
        assert_equal(
            reader.ts, reader_p.ts, "Last timestep is changed after pickling"
        )

    @pytest.mark.skip(reason="Cannot copy stream")
    def test_transformations_copy(self, ref, transformed):
        # this test checks if transformations are carried over a copy and if the
        # coordinates of the copy are equal to the original's
        v1 = np.float32((1, 1, 1))
        v2 = np.float32((0, 0, 0.33))
        new = transformed.copy()
        assert_equal(
            transformed.transformations,
            new.transformations,
            "transformations are not equal",
        )
        for i, ts in enumerate(new):
            ideal_coords = ref.iter_ts(i).positions + v1 + v2
            assert_array_almost_equal(
                ts.positions, ideal_coords, decimal=ref.prec
            )

    @pytest.mark.skip(
        reason="Timeseries currently requires n_frames to be known"
    )
    def test_timeseries_empty_asel(self, reader):
        with pytest.warns(
            UserWarning,
            match="Empty string to select atoms, empty group returned.",
        ):
            atoms = mda.Universe(reader.filename).select_atoms(None)
        with pytest.raises(ValueError, match="Timeseries requires at least"):
            reader.timeseries(asel=atoms)

    @pytest.mark.skip(
        reason="Timeseries currently requires n_frames to be known"
    )
    def test_timeseries_empty_atomgroup(self, reader):
        with pytest.warns(
            UserWarning,
            match="Empty string to select atoms, empty group returned.",
        ):
            atoms = mda.Universe(reader.filename).select_atoms(None)
        with pytest.raises(ValueError, match="Timeseries requires at least"):
            reader.timeseries(atomgroup=atoms)

    @pytest.mark.skip(
        reason="Timeseries currently requires n_frames to be known"
    )
    def test_timeseries_asel_warns_deprecation(self, reader):
        atoms = mda.Universe(reader.filename).select_atoms("index 1")
        with pytest.warns(DeprecationWarning, match="asel argument to"):
            timeseries = reader.timeseries(asel=atoms, order="fac")

    @pytest.mark.skip(
        reason="Timeseries currently requires n_frames to be known"
    )
    def test_timeseries_atomgroup(self, reader):
        atoms = mda.Universe(reader.filename).select_atoms("index 1")
        timeseries = reader.timeseries(atomgroup=atoms, order="fac")

    @pytest.mark.skip(
        reason="Timeseries currently requires n_frames to be known"
    )
    def test_timeseries_atomgroup_asel_mutex(self, reader):
        atoms = mda.Universe(reader.filename).select_atoms("index 1")
        with pytest.raises(ValueError, match="Cannot provide both"):
            timeseries = reader.timeseries(
                atomgroup=atoms, asel=atoms, order="fac"
            )

    @pytest.mark.skip("Cannot slice stream")
    def test_last_frame(self, ref, reader):
        ts = reader[-1]
        assert_timestep_almost_equal(ts, ref.last_frame, decimal=ref.prec)

    @pytest.mark.skip("Cannot slice stream")
    def test_go_over_last_frame(self, ref, reader):
        with pytest.raises(IndexError):
            reader[ref.n_frames + 1]

    @pytest.mark.skip("Cannot slice stream")
    def test_frame_jump(self, ref, reader):
        ts = reader[ref.jump_to_frame.frame]
        assert_timestep_almost_equal(ts, ref.jump_to_frame, decimal=ref.prec)

    @pytest.mark.skip("Cannot slice stream")
    def test_frame_jump_issue1942(self, ref, reader):
        """Test for issue 1942 (especially XDR on macOS)"""
        reader.rewind()
        try:
            for ii in range(ref.n_frames + 2):
                reader[0]
        except StopIteration:
            pytest.fail("Frame-seeking wrongly iterated (#1942)")

    def test_next_gives_second_frame(self, ref, reader):
        # don't recreate reader here as in inherited base test
        ts = reader.next()
        assert_timestep_almost_equal(ts, ref.second_frame, decimal=ref.prec)

    @pytest.mark.skip(
        reason="Stream isn't rewound after iteration- base reference is the same but it is the last frame"
    )
    def test_frame_collect_all_same(self, reader):
        # check that the timestep resets so that the base reference is the same
        # for all timesteps in a collection with the exception of memoryreader
        # and DCDReader
        if isinstance(reader, mda.coordinates.memory.MemoryReader):
            pytest.xfail("memoryreader allows independent coordinates")
        if isinstance(reader, mda.coordinates.DCD.DCDReader):
            pytest.xfail(
                "DCDReader allows independent coordinates."
                "This behaviour is deprecated and will be changed"
                "in 3.0"
            )
        collected_ts = []
        for i, ts in enumerate(reader):
            collected_ts.append(ts.positions)
        for array in collected_ts:
            assert_allclose(array, collected_ts[0])
