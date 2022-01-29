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

:Authors: Hugo MacDermott-Opeskin
:Year: 2021
:Copyright: GNU Public License v2
"""

import numpy as np
import MDAnalysis as mda
from . import base, core
from ..exceptions import NoDataError
from ..due import due, Doi

try:
    import pytng
except ImportError:
    HAS_PYTNG = False
else:
    HAS_PYTNG = True


class TNGReader(base.ReaderBase):
    r"""Reader for the TNG format."""
    format = 'TRR'
    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps',
             'force': 'kJ/(mol*nm)'}

    _special_blocks = ["TNG_TRAJ_BOX_SHAPE", "TNG_TRAJ_POSITIONS",
                       "TNG_TRAJ_VELOCITIES", "TNG_TRAJ_FORCES"]

    def __init__(self, filename, **kwargs):

        if not HAS_PYTNG:
            raise RuntimeError("Please install pytng")

        super(TNGReader, self).__init__(filename, **kwargs)

        self.filename = filename


        self._file_iterator = pytng.TNGFileIterator(self.filename,'r')
        self.n_atoms = self._file_iterator.n_atoms

        self._block_names = self._file_iterator.block_ids.keys()
        self._block_ids = self._file_iterator.block_ids.values()
        self._block_strides = self._file_iterator.block_strides
        self._data_frames   = self._file_iterator.n_data_frames

        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        self._has_box = "TNG_TRAJ_BOX_SHAPE" in self._block_names
        self._box = None

        self._has_positions = "TNG_TRAJ_POSITIONS" in self._block_names
        self._positions = None 
        
        self._has_velocities = "TNG_TRAJ_VELOCITIES" in self._block_names
        self._velocities = None

        self._has_forces = "TNG_TRAJ_FORCES" in self._block_names
        self._forces = None

        self._additional_blocks = [block for block in self._block_names if block not in self._special_blocks]
        self._additional_block_data = {block: None for block in self._additional_blocks}
        self._make_ndarrays()
        self._check_strides()
        self._frame = 0 
    
    def _make_ndarrays(self):
        if self._has_box:
            self._box = self._file_iterator.make_ndarray_for_block_from_name("TNG_TRAJ_BOX_SHAPE")
        if self._has_positions:
            self._positions = self._file_iterator.make_ndarray_for_block_from_name("TNG_TRAJ_POSITIONS")
        if self._has_velocities:
            self._velocities = self._file_iterator.make_ndarray_for_block_from_name("TNG_TRAJ_VELOCITIES")
        if self._has_forces:
            self._forces = self._file_iterator.make_ndarray_for_block_from_name("TNG_TRAJ_FORCES")
        
        for block in self._additional_blocks:
            self._additional_block_data[block] = self._file_iterator.make_ndarray_for_block_from_name(block)
    
    def _check_strides(self):
        if self._has_box and self._has_positions:
            if 
        
    def close(self):
        """close reader"""
        self._file_iterator._close()

    @staticmethod
    def parse_n_atoms(filename, **kwargs):
        with pytng.TNGFileIterator(filename 'r') as tng:
            n_atoms = tng.n_atoms
        return n_atoms

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        return len(self._file_iterator)

    def _reopen(self):
        """reopen trajectory"""
        self.ts.frame = 0
        self._frame = -1
        self._file_iterator._close()
        self._file_iterator._open(self.filename, 'r')

    def _read_frame(self, i):
        """read frame i"""
        raise NotImplementedError
        


    def _read_next_timestep(self, ts=None):
        """copy next frame into timestep"""
        raise NotImplementedError
        

    def Writer(self):
        return None

    
