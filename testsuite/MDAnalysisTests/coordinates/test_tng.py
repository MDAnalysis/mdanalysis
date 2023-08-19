# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysis.coordinates.TNG import HAS_PYTNG
from MDAnalysis.lib.mdamath import triclinic_box
from numpy.testing import assert_allclose, assert_equal

from MDAnalysisTests.coordinates.base import (
    BaseReference,
    MultiframeReaderTest,
)
from MDAnalysisTests.datafiles import (
    COORDINATES_TNG,
    COORDINATES_TOPOLOGY,
    TNG_traj,
    TNG_traj_gro,
    TNG_traj_uneven_blocks,
    TNG_traj_vels_forces,
)


@pytest.mark.skipif(HAS_PYTNG, reason="pytng present")
def test_pytng_not_present_raises():
    with pytest.raises(ImportError, match="please install pytng"):
        _ = mda.Universe(TNG_traj_gro, TNG_traj)


@pytest.mark.skipif(HAS_PYTNG, reason="pytng present")
def test_parse_n_atoms_no_pytng():
    with pytest.raises(ImportError, match="please install pytng"):
        mda.coordinates.TNG.TNGReader.parse_n_atoms(TNG_traj)


@pytest.mark.skipif(not HAS_PYTNG, reason="pytng not installed")
class TNGReference(BaseReference):
    """Reference synthetic trajectory that was
    copied from test_xdr.TRReference"""

    def __init__(self):
        super(TNGReference, self).__init__()
        self.trajectory = COORDINATES_TNG
        self.topology = COORDINATES_TOPOLOGY
        self.reader = mda.coordinates.TNG.TNGReader
        self.writer = mda.coordinates.TNG.TNGReader.Writer
        self.ext = "tng"
        self.changing_dimensions = True
        self.prec = 4

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


@pytest.mark.filterwarnings("ignore:Empty string")
@pytest.mark.filterwarnings("ignore:Stride of block")
@pytest.mark.skipif(not HAS_PYTNG, reason="pytng not installed")
class TestTNGCoordinatesTraj(MultiframeReaderTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return TNGReference()

    def test_get_writer_1(self, reader):
        with pytest.raises(
            NotImplementedError,
            match="There is currently no writer for TNG files",
        ):
            reader.Writer()

    def test_get_writer_2(self, reader):
        with pytest.raises(
            NotImplementedError,
            match="There is currently no writer for TNG files",
        ):
            reader.Writer()


@pytest.mark.skipif(not HAS_PYTNG, reason="pytng not installed")
def test_tng_traj_uneven_blocks():
    with pytest.raises(IOError, match="Strides of TNG special blocks"):
        _ = mda.Universe(TNG_traj_gro, TNG_traj_uneven_blocks)


@pytest.mark.filterwarnings("ignore:Failed read for block")
@pytest.mark.skipif(not HAS_PYTNG, reason="pytng not installed")
class TestTNGTraj(object):

    _n_atoms = 1000
    _n_frames = 101
    _stride = 5000

    prec = 5

    # these values all taken from GMX dump / 10 to get MDA units
    _pos_frame_0_first_3_atoms = (
        np.array(
            [
                [2.53300e00, 1.24400e00, 3.50600e00],
                [8.30000e-01, 2.54400e00, 3.44800e00],
                [1.09100e00, 1.10000e-01, 3.12900e00],
            ]
        )
        * 10
    )

    _pos_frame_100_first_3_atoms = (
        np.array(
            [
                [4.40000e-01, 3.89000e-01, 1.37400e00],
                [1.43200e00, 1.64900e00, 2.93900e00],
                [2.01500e00, 2.10300e00, 2.65700e00],
            ]
        )
        * 10
    )

    _box_frame_0 = (
        np.array(
            [
                3.60140e00,
                0.00000e00,
                0.00000e00,
                0.00000e00,
                3.60140e00,
                0.00000e00,
                0.00000e00,
                0.00000e00,
                3.60140e00,
            ]
        ).reshape(3, 3)
        * 10
    )

    _box_frame_100 = (
        np.array(
            [
                3.60140e00,
                0.00000e00,
                0.00000e00,
                0.00000e00,
                3.60140e00,
                0.00000e00,
                0.00000e00,
                0.00000e00,
                3.60140e00,
            ]
        ).reshape(3, 3)
        * 10
    )

    _box_frame_100 = (
        np.array(
            [
                3.58965e00,
                0.00000e00,
                0.00000e00,
                0.00000e00,
                3.58965e00,
                0.00000e00,
                0.00000e00,
                0.00000e00,
                3.58965e00,
            ]
        ).reshape(3, 3)
        * 10
    )

    @pytest.fixture(scope="class")
    def universe(self):
        return mda.Universe(TNG_traj_gro, TNG_traj)

    def test_n_atoms(self, universe):
        assert_equal(universe.trajectory.n_atoms, self._n_atoms)

    def test_n_frames(self, universe):
        assert_equal(
            universe.trajectory.n_frames,
            self._n_frames,
            "wrong number of frames in TNG file",
        )

    def test_block_names(self, universe):
        """Check the block names in the file"""
        assert "TNG_TRAJ_BOX_SHAPE" in universe.trajectory.blocks
        assert "TNG_TRAJ_POSITIONS" in universe.trajectory.blocks
        assert "TNG_GMX_LAMBDA" in universe.trajectory.blocks

    def test_special_blocks(self, universe):
        """Check the position and box special blocks are present"""
        assert "TNG_TRAJ_BOX_SHAPE" in universe.trajectory.special_blocks
        assert "TNG_TRAJ_POSITIONS" in universe.trajectory.special_blocks

    def test_additional_blocks(self, universe):
        """Check the lambda special block is present"""
        assert "TNG_GMX_LAMBDA" in universe.trajectory.additional_blocks

    def check_strides(self, universe):
        """Check the stride of trajectory frames in integrator steps"""
        assert universe.trajectory._global_stride == self._stride

    def test_initial_frame_is_0(self, universe):
        assert_equal(
            universe.trajectory.ts.frame,
            0,
            "initial frame is not 0 but {0}".format(universe.trajectory.ts.frame),
        )

    def test_starts_with_first_frame(self, universe):
        """Test that coordinate arrays are filled as soon as the trajectory
        has been opened."""
        assert np.any(universe.atoms.positions > 0)

    def test_rewind(self, universe):
        trj = universe.trajectory
        trj.next()
        trj.next()  # for readers that do not support indexing
        assert_equal(trj.ts.frame, 2, "failed to forward to frame 2 (frameindex 2)")
        trj.rewind()
        assert_equal(trj.ts.frame, 0, "failed to rewind to first frame")
        assert np.any(universe.atoms.positions > 0)

    def test_full_slice(self, universe):
        trj_iter = universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(universe.trajectory.n_frames))

    def test_random_access(self, universe):
        pos1 = universe.atoms[0].position
        universe.trajectory.next()
        universe.trajectory.next()
        pos3 = universe.atoms[0].position

        universe.trajectory[0]

        assert_equal(universe.atoms[0].position, pos1)

        universe.trajectory[2]

        assert_equal(universe.atoms[0].position, pos3)

    def test_frame_overrun(self, universe):
        with pytest.raises(IndexError, match="exceeds length of trajectory"):
            universe.trajectory[101]

    def test_positions_first_frame(self, universe):
        pos = universe.trajectory[0].positions
        assert_allclose(
            pos[0:3, :], self._pos_frame_0_first_3_atoms, rtol=10**-self.prec
        )

    def test_box_first_frame(self, universe):
        dims = universe.trajectory[0].dimensions
        assert_allclose(dims, triclinic_box(*self._box_frame_0), rtol=10**-self.prec)

    def test_positions_last_frame(self, universe):
        pos = universe.trajectory[100].positions
        assert_allclose(
            pos[0:3, :],
            self._pos_frame_100_first_3_atoms,
            rtol=10**-self.prec,
        )

    def test_box_last_frame(self, universe):
        dims = universe.trajectory[100].dimensions
        assert_allclose(
            dims, triclinic_box(*self._box_frame_100), rtol=10**-self.prec
        )

    @pytest.mark.parametrize("frame", [0, 20, 50, 100])
    def test_step(self, universe, frame):
        ts = universe.trajectory[frame]
        step = ts.data["step"]
        assert step == ts.frame * universe.trajectory._global_stride

    def test_lambda_in_ts(self, universe):
        ts = universe.trajectory[10]
        assert "TNG_GMX_LAMBDA" in ts.data.keys()
        assert isinstance(ts.data["TNG_GMX_LAMBDA"], np.ndarray)
        assert_equal(ts.data["TNG_GMX_LAMBDA"], np.asarray([[0]], dtype=np.float32))

    def test_read_box_fail_strange_step(self, universe):
        stepnum = 123  # step number with no data
        iterator_step = universe.trajectory._file_iterator.read_step(stepnum)
        with pytest.raises(IOError, match="Failed to read box from TNG file"):
            universe.trajectory._frame_to_ts(iterator_step, universe.trajectory.ts)

    def test_read_pos_fail_strange_step(self, universe):
        stepnum = 123  # step number with no data
        iterator_step = universe.trajectory._file_iterator.read_step(stepnum)
        # set _has_box to False to trigger position reading error
        universe.trajectory._has_box = False
        with pytest.raises(IOError, match="Failed to read positions from TNG file"):
            universe.trajectory._frame_to_ts(iterator_step, universe.trajectory.ts)

    def test_additional_block_read_fails(self, universe):
        stepnum = 123  # step number with no data
        iterator_step = universe.trajectory._file_iterator.read_step(stepnum)
        # set has_box, has_pos, to false to trigger GMX_LAMBDA reading error
        universe.trajectory._has_box = False
        universe.trajectory._has_positions = False
        # doesn't have velocities or forces
        with pytest.raises(
            IOError, match="Failed to read additional block TNG_GMX_LAMBDA"
        ):
            universe.trajectory._frame_to_ts(iterator_step, universe.trajectory.ts)

    def test_parse_n_atoms(self, universe):
        assert universe.trajectory.parse_n_atoms(TNG_traj) == self._n_atoms


@pytest.mark.filterwarnings("ignore:Off stride read for block")
@pytest.mark.skipif(not HAS_PYTNG, reason="pytng not installed")
class TestTNGTraj_vels_forces(object):

    _n_atoms = 1000
    _n_frames = 51
    _stride = 10

    prec = 5

    @pytest.fixture(scope="class")
    def universe(self):
        return mda.Universe(TNG_traj_gro, TNG_traj_vels_forces)

    def test_n_atoms(self, universe):
        assert_equal(universe.trajectory.n_atoms, self._n_atoms)

    def test_n_frames(self, universe):
        assert_equal(
            universe.trajectory.n_frames,
            self._n_frames,
            "wrong number of frames in TNG file",
        )

    def test_block_names(self, universe):
        """Check the block names in the file"""
        assert "TNG_TRAJ_BOX_SHAPE" in universe.trajectory.blocks
        assert "TNG_TRAJ_POSITIONS" in universe.trajectory.blocks
        assert "TNG_TRAJ_VELOCITIES" in universe.trajectory.blocks
        assert "TNG_TRAJ_FORCES" in universe.trajectory.blocks
        assert "TNG_GMX_LAMBDA" in universe.trajectory.blocks

    def test_special_blocks(self, universe):
        """Check the position and box special blocks are present"""
        assert "TNG_TRAJ_BOX_SHAPE" in universe.trajectory.special_blocks
        assert "TNG_TRAJ_POSITIONS" in universe.trajectory.special_blocks
        assert "TNG_TRAJ_VELOCITIES" in universe.trajectory.special_blocks
        assert "TNG_TRAJ_FORCES" in universe.trajectory.special_blocks

    def check_strides(self, universe):
        """Check the stride of trajectory frames in integrator steps"""
        assert universe.trajectory._global_stride == self._stride

    def test_read_vels_fail_strange_step(self, universe):
        stepnum = 123  # step number with no data
        iterator_step = universe.trajectory._file_iterator.read_step(stepnum)
        # set _has_* attrs to False to trigger velocities reading error
        universe.trajectory._has_box = False
        universe.trajectory._has_positions = False
        with pytest.raises(IOError, match="Failed to read velocities from TNG file"):
            universe.trajectory._frame_to_ts(iterator_step, universe.trajectory.ts)

    def test_read_force_fail_strange_step(self, universe):
        stepnum = 123  # step number with no data
        iterator_step = universe.trajectory._file_iterator.read_step(stepnum)
        # set _has_* attrs to False to trigger forces reading error
        universe.trajectory._has_box = False
        universe.trajectory._has_positions = False
        universe.trajectory._has_velocities = False
        with pytest.raises(IOError, match="Failed to read forces from TNG file"):
            universe.trajectory._frame_to_ts(iterator_step, universe.trajectory.ts)


@pytest.mark.skipif(not HAS_PYTNG, reason="pytng not installed")
def test_writer_raises_notimpl():
    u = mda.Universe(TNG_traj_gro, TNG_traj)
    with pytest.raises(
        NotImplementedError,
        match="There is currently no writer for TNG files",
    ):
        u.trajectory.Writer()
