# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""\
=========================================================================
Reading trajectories from memory --- :mod:`MDAnalysis.coordinates.memory`
=========================================================================

:Author: Wouter Boomsma
:Year: 2016
:Copyright: GNU Public License v2
:Maintainer: Wouter Boomsma <wb@di.ku.dk>, wouterboomsma on github


.. versionadded:: 0.16.0

The module contains a trajectory reader that operates on an array in
memory, rather than reading from file. This makes it possible to use
operate on raw coordinate using existing MDAnalysis tools. In
addition, it allows the user to make changes to the coordinates in a
trajectory (e.g. through
:attr:`MDAnalysis.core.groups.AtomGroup.positions`) without having
to write the entire state to file.


How to use the :class:`MemoryReader`
====================================

The :class:`MemoryReader` can be used to either directly generate a
trajectory as a numpy array or by transferring an existing trajectory
to memory.

In-memory representation of arbitrary trajectories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If sufficient memory is available to hold a whole trajectory in memory
then analysis can be sped up substantially by transferring the
trajectory to memory.

The most straightforward use of the :class:`MemoryReader` is to simply
use the ``in_memory=True`` flag for the
:class:`~MDAnalysis.core.universe.Universe` class, which
automatically transfers a trajectory to memory::

 import MDAnalysis as mda
 from MDAnalysisTests.datafiles import TPR, XTC

 universe = mda.Universe(TPR, XTC, in_memory=True)

Of course, sufficient memory has to be available to hold the whole
trajectory.


Switching a trajectory to an in-memory representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The decision to transfer the trajectory to memory can be made at any
time with the
:meth:`~MDAnalysis.core.universe.Universe.transfer_to_memory` method
of a :class:`~MDAnalysis.core.universe.Universe`::

    universe = mda.Universe(TPR, XTC)
    universe.transfer_to_memory()

This operation may take a while (with `verbose=True` a progress bar is
displayed) but then subsequent operations on the trajectory directly
operate on the in-memory array and will be very fast.


Constructing a Reader from a numpy array
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`MemoryReader` provides great flexibility because it
becomes possible to create a :class:`~MDAnalysis.core.universe.Universe` directly
from a numpy array.

A simple example consists of a new universe created from the array
extracted from a DCD
:meth:`~MDAnalysis.coordinates.DCD.DCDReader.timeseries`::

    import MDAnalysis as mda
    from MDAnalysisTests.datafiles import DCD, PSF
    from MDAnalysis.coordinates.memory import MemoryReader

    universe = mda.Universe(PSF, DCD)

    coordinates = universe.trajectory.timeseries(universe.atoms)
    universe2 = mda.Universe(PSF, coordinates, format=MemoryReader, order='afc')


.. _create-in-memory-trajectory-with-AnalysisFromFunction:

.. rubric:: Creating an in-memory trajectory with
            :func:`~MDAnalysis.analysis.base.AnalysisFromFunction`

The :meth:`~MDAnalysis.coordinates.DCD.DCDReader.timeseries` is
currently only implemented for the
:class:`~MDAnalysis.coordinates.DCD.DCDReader`. However, the
:func:`MDAnalysis.analysis.base.AnalysisFromFunction` can provide the
same functionality for any supported trajectory format::

  import MDAnalysis as mda
  from MDAnalysis.tests.datafiles import PDB, XTC

  from MDAnalysis.coordinates.memory import MemoryReader
  from MDAnalysis.analysis.base import AnalysisFromFunction

  u = mda.Universe(PDB, XTC)

  coordinates = AnalysisFromFunction(lambda ag: ag.positions.copy(),
                                     u.atoms).run().results
  u2 = mda.Universe(PDB, coordinates, format=MemoryReader)

.. _creating-in-memory-trajectory-label:

Creating an in-memory trajectory of a sub-system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Creating a trajectory for just a selection of an existing trajectory
requires the transfer of the appropriate coordinates as well as
creation of a topology of the sub-system. For the latter one can use
the :func:`~MDAnalysis.core.universe.Merge` function, for the former
the :meth:`~MDAnalysis.core.universe.Universe.load_new` method of a
:class:`~MDAnalysis.core.universe.Universe` together with the
:class:`MemoryReader`. In the following, an in-memory trajectory of
only the protein is created::

  import MDAnalysis as mda
  from MDAnalysis.tests.datafiles import PDB, XTC

  from MDAnalysis.coordinates.memory import MemoryReader
  from MDAnalysis.analysis.base import AnalysisFromFunction

  u = mda.Universe(PDB, XTC)
  protein = u.select_atoms("protein")

  coordinates = AnalysisFromFunction(lambda ag: ag.positions.copy(),
                                     protein).run().results
  u2 = mda.Merge(protein)            # create the protein-only Universe
  u2.load_new(coordinates, format=MemoryReader)

The protein coordinates are extracted into ``coordinates`` and then
the in-memory trajectory is loaded from these coordinates. In
principle, this could have all be done in one line::

  u2 = mda.Merge(protein).load_new(
           AnalysisFromFunction(lambda ag: ag.positions.copy(),
                                protein).run().results,
           format=MemoryReader)

The new :class:`~MDAnalysis.core.universe.Universe` ``u2`` can be used
to, for instance, write out a new trajectory or perform fast analysis
on the sub-system.


Classes
=======

.. autoclass:: Timestep
   :members:
   :inherited-members:

.. autoclass:: MemoryReader
   :members:
   :inherited-members:

"""
from __future__ import absolute_import
import logging
import errno
import numpy as np

from . import base


class Timestep(base.Timestep):
    """
    Timestep for the :class:`MemoryReader`

    Overrides the positions property in
    :class:`MDAnalysis.coordinates.base.Timestep` to use avoid
    duplication of the array.

    """

    @property
    def positions(self):
        return base.Timestep.positions.fget(self)

    @positions.setter
    def positions(self, new):
        self.has_positions = True
        # Use reference to original rather than a copy
        self._pos = new


class MemoryReader(base.ProtoReader):
    """
    MemoryReader works with trajectories represented as numpy arrays.

    A trajectory reader interface to a numpy array of the coordinates.
    For compatibility with the timeseries interface, support is provided for
    specifying the order of columns through the format option.

    .. versionadded:: 0.16.0

    """

    format = 'MEMORY'
    _Timestep = Timestep

    def __init__(self, coordinate_array, order='fac',
                 dimensions=None, dt=1, filename=None, **kwargs):
        """
        Parameters
        ----------
        coordinate_array : numpy.ndarray
            The underlying array of coordinates
        order : {"afc", "acf", "caf", "fac", "fca", "cfa"} (optional)
            the order/shape of the return data array, corresponding
            to (a)tom, (f)rame, (c)oordinates all six combinations
            of 'a', 'f', 'c' are allowed ie "fac" - return array
            where the shape is (frame, number of atoms,
            coordinates).
        dimensions: [A, B, C, alpha, beta, gamma] (optional)
            unitcell dimensions (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)
            lengths *A*, *B*, *C* are in the MDAnalysis length unit (Å), and
            angles are in degrees.
        dt: float (optional)
            The time difference between frames (ps).  If :attr:`time`
            is set, then `dt` will be ignored.
        filename: string (optional)
            The name of the file from which this instance is created. Set to ``None``
            when created from an array

        Note
        ----
        At the moment, only a fixed `dimension` is supported, i.e., the same
        unit cell for all frames in `coordinate_array`. See issue `#1041`_.

        .. _`#1041`: https://github.com/MDAnalysis/mdanalysis/issues/1041

        """

        super(MemoryReader, self).__init__()

        self.filename = filename
        self.stored_order = order
        self.set_array(np.asarray(coordinate_array), order)
        self.n_frames = \
            self.coordinate_array.shape[self.stored_order.find('f')]
        self.n_atoms = \
            self.coordinate_array.shape[self.stored_order.find('a')]

        provided_n_atoms = kwargs.pop("n_atoms", None)
        if (provided_n_atoms is not None and
            provided_n_atoms != self.n_atoms):
                raise ValueError("The provided value for n_atoms ({}) "
                                 "does not match the shape of the coordinate "
                                 "array ({})"
                                 .format(provided_n_atoms, self.n_atoms))

        self.ts = self._Timestep(self.n_atoms, **kwargs)
        self.ts.dt = dt
        if dimensions is not None:
            self.ts.dimensions = dimensions
        self.ts.frame = -1
        self.ts.time = -1
        self._read_next_timestep()

    def set_array(self, coordinate_array, order='fac'):
        """
        Set underlying array in desired column order.

        Parameters
        ----------
        coordinate_array : :class:`~numpy.ndarray` object
            The underlying array of coordinates
        order : {"afc", "acf", "caf", "fac", "fca", "cfa"} (optional)
            the order/shape of the return data array, corresponding
            to (a)tom, (f)rame, (c)oordinates all six combinations
            of 'a', 'f', 'c' are allowed ie "fac" - return array
            where the shape is (frame, number of atoms,
            coordinates).
        """
        # Only make copy if not already in float32 format
        self.coordinate_array = coordinate_array.astype('float32', copy=False)
        self.stored_format = order

    def get_array(self):
        """
        Return underlying array.
        """
        return self.coordinate_array

    def _reopen(self):
        """Reset iteration to first frame"""
        self.ts.frame = -1
        self.ts.time = -1

    def timeseries(self, asel=None, start=0, stop=-1, step=1, format='afc'):
        """Return a subset of coordinate data for an AtomGroup in desired
        column order/format. If no selection is given, it will return a view of
        the underlying array, while a copy is returned otherwise.

        Parameters
        ---------
        asel : AtomGroup (optional)
            Atom selection. Defaults to ``None``, in which case the full set of
            coordinate data is returned. Note that in this case, a view
            of the underlying numpy array is returned, while a copy of the
            data is returned whenever `asel` is different from ``None``.
        start : int (optional)
        stop : int (optional)
        step : int (optional)
            range of trajectory to access, `start` and `stop` are *inclusive*
        format : {"afc", "acf", "caf", "fac", "fca", "cfa"} (optional)
            the order/shape of the return data array, corresponding
            to (a)tom, (f)rame, (c)oordinates all six combinations
            of 'a', 'f', 'c' are allowed ie "fac" - return array
            where the shape is (frame, number of atoms,
            coordinates).

        Note
        ----
        The `format` parameter name is used to mimic the
        :class:`MDAnalysis.coordinates.DCD.timeseries` interface. It is
        identical to the `order` parameter for :class:`MemoryReader`. In a
        future version, `format` will be renamed to `order`.
        """
        # Renaming 'format' to 'order' here for internal consistency in this class
        order = format

        array = self.get_array()
        if order == self.stored_order:
            pass
        elif order[0] == self.stored_order[0]:
            array = np.swapaxes(array, 1, 2)
        elif order[1] == self.stored_order[1]:
            array = np.swapaxes(array, 0, 2)
        elif order[2] == self.stored_order[2]:
            array = np.swapaxes(array, 0, 1)
        elif order[0] == self.stored_order[1]:
            array = np.swapaxes(array, 1, 0)
            array = np.swapaxes(array, 1, 2)
        elif order[0] == self.stored_order[2]:
            array = np.swapaxes(array, 2, 0)
            array = np.swapaxes(array, 1, 2)

        a_index = order.find('a')
        f_index = order.find('f')
        stop_index = stop+1
        if stop_index == 0:
            stop_index = None
        basic_slice = ([slice(None)] * f_index +
                       [slice(start, stop_index, step)] +
                       [slice(None)] * (2-f_index))

        # Return a view if either:
        #   1) asel is None
        #   2) asel corresponds to the selection of all atoms.
        array = array[basic_slice]
        if (asel is None or asel is asel.universe.atoms):
            return array
        else:
            # If selection is specified, return a copy
            return array.take(asel.indices, a_index)

    def _read_next_timestep(self, ts=None):
        """copy next frame into timestep"""

        if self.ts.frame >= self.n_frames-1:
            raise IOError(errno.EIO, 'trying to go over trajectory limit')
        if ts is None:
            ts = self.ts
        ts.frame += 1
        f_index = self.stored_order.find('f')
        basic_slice = ([slice(None)]*(f_index) +
                       [self.ts.frame] +
                       [slice(None)]*(2-f_index))
        ts.positions = self.coordinate_array[basic_slice]

        ts.time = self.ts.frame*self.dt
        return ts

    def _read_frame(self, i):
        """read frame i"""
        # Frame number is incremented to zero by _read_next_timestep()
        self.ts.frame = i - 1
        return self._read_next_timestep()

    def __repr__(self):
        """String representation"""
        return ("<{cls} with {nframes} frames of {natoms} atoms>"
                "".format(
                    cls=self.__class__.__name__,
                    nframes=self.n_frames,
                    natoms=self.n_atoms
                ))
