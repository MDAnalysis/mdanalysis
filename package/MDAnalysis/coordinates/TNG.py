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

    def __init__(self, filename, **kwargs):

        if not HAS_PYTNG:
            raise RuntimeError("Please install pytng")

        super(TNGReader, self).__init__(filename, **kwargs)
        self.filename = filename
        self._frame = 0 
        self._file_iterator = pytng.TNGFileIterator(self.filename,'r')
    


   def close(self):
        """close reader"""
        self._file_iterator._close()

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


    def _read_next_timestep(self, ts=None):

        

    def Writer(self):
        return None

    
