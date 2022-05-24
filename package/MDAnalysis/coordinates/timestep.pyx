# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
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
#

import weakref
import copy
import numpy as np
cimport numpy as cnp
cnp.import_array()

from libc.stdint cimport uint64_t
from libcpp cimport bool

from . import core
from .. import NoDataError
from .. import (
    _READERS, _READER_HINTS,
    _SINGLEFRAME_WRITERS,
    _MULTIFRAME_WRITERS,
    _CONVERTERS
)
from .. import units
from ..auxiliary.base import AuxReader
from ..auxiliary.core import auxreader
from ..lib.util import asiterable, Namespace

cdef class Timestep:

    order = 'C'

    cdef uint64_t n_atoms
    cdef uint64_t frame 

    cdef bool _has_dimensions

    cdef bool _has_positions

    cdef bool _has_velocities

    cdef bool _has_forces

    cdef cnp.ndarray _dimensions
    cdef cnp.ndarray _positions
    cdef cnp.ndarray _velocities
    cdef cnp.ndarray _forces

    cdef object _dtype
    cdef public dict data
    cdef public object aux




    def __cinit__(self, uint64_t n_atoms, dtype=np.float32, **kwargs):
        # c++ level objects
        self.n_atoms =  n_atoms
        self.frame = -1
        self._has_dimensions = False
        self._has_positions = False
        self._has_velocities = False
        self._has_forces = False

        # match original can this be removed?
        self._dimensions = np.zeros(6, dtype=np.float32)

    def __init__(self, uint64_t n_atoms, dtype=np.float32, **kwargs):
        #python objects
        self._dtype = dtype
        self.data = {}

        for att in ('dt', 'time_offset'):
            try:
                self.data[att] = kwargs[att]
            except KeyError:
                pass
        try:
            # do I have a hook back to the Reader?
            self._reader = weakref.ref(kwargs['reader'])
        except KeyError:
            pass

        # set up aux namespace for adding auxiliary data
        self.aux = Namespace()
        
    
    def __dealloc__(self):
            pass

    @property
    def dtype(self):
        return self._dtype


    @property
    def has_positions(self):
        return self._has_positions

    @property
    def has_dimensions(self):
        return self._has_dimensions
    
    @property
    def has_velocities(self):
        return self._has_velocities
    
    @property
    def has_forces(self):
        return self._has_forces


    @property
    def positions(self):
        if self._has_positions:
            return self._positions
        else:
            raise ValueError("This Timestep has no position information")

 
    @positions.setter
    def positions(self,  cnp.ndarray new_positions):
        # force C contig memory order
        self._positions = np.ascontiguousarray(new_positions).copy()
        self._has_positions = True



    @property
    def dimensions(self):
        if self._has_dimensions:
           return self._dimensions
        else:
            raise ValueError("This Timestep has no dimension information")

    
    @dimensions.setter

    def dimensions(self, cnp.ndarray new_dimensions):
        # force C contig memory order
        self._dimensions = np.ascontiguousarray(new_dimensions).copy()
        self._has_dimensions = True


    @property
    def velocities(self):
        if self._has_velocities:
            return self._velocities
        else:
            raise ValueError("This Timestep has no velocities information")


 
    @velocities.setter

    def velocities(self,  cnp.ndarray new_velocities):
        # force C contig memory order
        self._velocities = np.ascontiguousarray(new_velocities).copy()
        self._has_velocities = True



    @property
    def forces(self):
        if self._has_forces:
          return self._forces
        else:
            raise ValueError("This Timestep has no force information")


 
    @forces.setter
    def forces(self,  cnp.ndarray new_forces):
        # force C contig memory order
        self._forces = np.ascontiguousarray(new_forces).copy()
        self._has_forces = True


    @classmethod
    def from_timestep(cls, other, **kwargs):
        """Create a copy of another Timestep, in the format of this Timestep

        .. versionadded:: 0.11.0
        """
        ts = cls(other.n_atoms,
                 positions=other.has_positions,
                 velocities=other.has_velocities,
                 forces=other.has_forces,
                 **kwargs)
        ts.frame = other.frame
        ts.dimensions = other.dimensions
        try:
            ts.positions = other.positions.copy(order=cls.order)
        except NoDataError:
            pass
        try:
            ts.velocities = other.velocities.copy(order=cls.order)
        except NoDataError:
            pass
        try:
            ts.forces = other.forces.copy(order=cls.order)
        except NoDataError:
            pass

        # Optional attributes that don't live in .data
        # should probably iron out these last kinks
        for att in ('_frame',):
            try:
                setattr(ts, att, getattr(other, att))
            except AttributeError:
                pass

        if hasattr(ts, '_reader'):
            other._reader = weakref.ref(ts._reader())

        ts.data = copy.deepcopy(other.data)

        return ts