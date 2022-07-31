# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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

"""
TNG Trajectory IO --- :mod:`MDAnalysis.coordinates.TNG`
=======================================================

Some stuff
"""

import numpy as np
import MDAnalysis as mda
import warnings
from . import base
from ..due import due, Doi
from ..lib.mdamath import  triclinic_box


try:
    import pytng
except ImportError:
    HAS_PYTNG = False
else:
    HAS_PYTNG = True


class TNGReader(base.ReaderBase):
    r"""Reader for the TNG format"""

    format = 'TNG'
    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps',
             'force': 'kJ/(mol*nm)'}

    _box_blockname = "TNG_TRAJ_BOX_SHAPE"
    _positions_blockname = "TNG_TRAJ_POSITIONS"
    _velocities_blockname = "TNG_TRAJ_VELOCITIES"
    _forces_blockname = "TNG_TRAJ_FORCES"
    _special_blocks = [_box_blockname, _positions_blockname,
                       _velocities_blockname, _forces_blockname]

    def __init__(self, filename, **kwargs):

        if not HAS_PYTNG:
            raise RuntimeError("To read TNG files please install pytng")

        super(TNGReader, self).__init__(filename, **kwargs)

        self.filename = filename

        self._file_iterator = pytng.TNGFileIterator(self.filename, 'r')
        self.n_atoms = self._file_iterator.n_atoms

        self._block_names = self._file_iterator.block_ids.keys()
        self._block_dictionary = self._file_iterator.block_ids
        self._block_strides = self._file_iterator.block_strides
        self._data_frames = self._file_iterator.n_data_frames
        self._special_block_present = {
            k: False for k in self._special_blocks}

        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        self._has_box = self._box_blockname in self._block_names
        if self._has_box:
            self._special_block_present[self._box_blockname] = True

        self._has_positions = self._positions_blockname in self._block_names
        if self._has_positions:
            self._special_block_present[self._positions_blockname] = True

        self._has_velocities = self._velocities_blockname in self._block_names
        if self._has_velocities:
            self._special_block_present[self._velocities_blockname] = True

        self._has_forces = self._forces_blockname in self._block_names
        if self._has_forces:
            self._special_block_present[self._forces_blockname] = True

        self._additional_blocks = [
            block for block in self._block_names if block not in self._special_blocks]
        self._check_strides_and_frames()
        self._frame = 0

    def _check_strides_and_frames(self):
        strides = []
        n_data_frames = []

        []
        for block in self.special_blocks:
            stride = self._block_strides[block]
            strides.append(stride)
            n_data_frame = self._data_frames[block]
            n_data_frames.append(n_data_frame)

        if not all(element == strides[0] for element in strides):
            raise IOError("Strides of TNG special blocks not equal,"
                          " file cannot be read")

        if not all(element == n_data_frames[0] for element in n_data_frames):
            raise IOError("Number of data containing frames for TNG special"
                          " blocks not equal, file cannot be read")

        self._global_stride = strides[0]
        self._n_frames = n_data_frames[0]

        self._additional_blocks_to_read = []
        for block in self._additional_blocks:
            stride_add = self._block_strides[block]
            n_data_frame_add = self._data_frames[block]
            if (stride_add != self._global_stride) or (n_data_frame_add != self.n_frames):
                warnings.warn("TNG additional block {block} does not match"
                              " strides of other blocks and will not be read")
            else:
                self._additional_blocks_to_read.append(block)

    def close(self):
        """close reader"""
        self._file_iterator._close()

    @staticmethod
    def parse_n_atoms(filename, **kwargs):
        with pytng.TNGFileIterator(filename, 'r') as tng:
            n_atoms = tng.n_atoms
        return n_atoms

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        return self._n_frames

    @property
    def special_blocks(self):
        "list of the special blocks that are in the file"
        return [k for k,v in self._special_block_present.items() if v]

    def _reopen(self):
        """reopen trajectory"""
        self.ts.frame = 0
        self._frame = -1
        self._file_iterator._close()
        self._file_iterator._open(self.filename, 'r')

    def _frame_to_step(self, frame):
        return frame * self._global_stride

    def _read_frame(self, i):
        """read frame i"""
        self._frame = i - 1
        ts = self._read_next_frame()
        return ts


    def _read_next_timestep(self, ts=None):
        """Read next frame into a timestep"""
        if self._frame == self.n_frames - 1:
            raise IOError('trying to go over trajectory limit')
        if ts is None:
            ts = self.ts
        # convert from frames to integrator steps
        step = self._frame_to_step(self._frame)
        iterator_step = self._file_iterator.read_step(step)
        self._frame += 1
        ts = self._frame_to_ts(iterator_step, ts)
        return ts

    def _frame_to_ts(self, curr_step, ts):
        """convert a TNGIteratorStep to an MDA Timestep"""

        ts.frame = self._frame
        ts.time = curr_step.get_time()
        ts.data['step'] = curr_step.step

        if self._has_box:
            box = self._file_iterator.make_ndarray_for_block_from_name(self._box_blockname)
            curr_step.get_box(box)
            ts.dimensions = triclinic_box(box)
            if not curr_step.read_success():
                raise IOError("Failed to read box from TNG file")
        if self._has_positions:
            curr_step.get_positions(ts.positions)
            if not curr_step.read_success():
                raise IOError("Failed to read positions from TNG file")
        if self._has_velocities:
            curr_step.get_velocities(ts.velocities)
            if not curr_step.read_velocities():
                raise IOError("Failed to read velocities from TNG file")
        if self._has_forces:
            curr_step.get_forces(ts.forces)
            if not curr_step.read_success():
                raise IOError("Failed to read forces from TNG file")

        for block in self._additional_blocks_to_read:
            block_data = self._file_iterator.make_ndarray_for_block_from_name(
                block)
            ts.data[block] = curr_step.get_blockid(
                self._block_dictionary[block], block_data)
            if not curr_step.read_success():
                raise IOError(
                    f"Failed to read additional block {block} from TNG file")

        return ts

    def Writer(self):
        raise NotImplementedError("There is currently no writer for TNG files")
