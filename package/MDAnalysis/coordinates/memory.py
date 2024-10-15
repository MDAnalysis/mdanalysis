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
memory, rather than reading from file. This makes it possible to
operate on raw coordinates using existing MDAnalysis tools. In
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
                                     u.atoms).run().results['timeseries']
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
                                     protein).run().results['timeseries']
  u2 = mda.Merge(protein)            # create the protein-only Universe
  u2.load_new(coordinates, format=MemoryReader)

The protein coordinates are extracted into ``coordinates`` and then
the in-memory trajectory is loaded from these coordinates. In
principle, this could have all be done in one line::

  u2 = mda.Merge(protein).load_new(
           AnalysisFromFunction(lambda ag: ag.positions.copy(),
                                protein).run().results['timeseries'],
           format=MemoryReader)

The new :class:`~MDAnalysis.core.universe.Universe` ``u2`` can be used
to, for instance, write out a new trajectory or perform fast analysis
on the sub-system.


Classes
=======


.. autoclass:: MemoryReader
   :members:
   :inherited-members:

"""
import logging
import errno
import numpy as np
import warnings
import copy

from . import base
from .timestep import Timestep


# These methods all pass in an existing *view* onto a larger array
def _replace_positions_array(ts, new):
    """Replace the array of positions

    Replaces the array of positions by another array.


    Note
    ----
    The behavior of :meth:`_replace_positions_array` is different from the
    behavior of the :attr:`position` property that replaces the **content**
    of the array. The :meth:`_replace_positions_array` method should only be
    used to set the positions to a different frame in
    :meth:`MemoryReader._read_next_timestep`; there, the memory reader sets
    the positions to a view of the correct frame.  Modifying the positions
    for a given frame should be done with the :attr:`positions` attribute
    that does not break the link between the array of positions in the time
    step and the :attr:`MemoryReader.coordinate_array`.


    .. versionadded:: 0.19.0
    .. versionchanged:: 2.0.0
       This function, and the _repalace helper functions for velocities,
       forces, and dimensions, have been moved out of the now removed
       custom timestep object for :class:`MemoryReader`.
    """
    ts.has_positions = True
    ts._pos = new


def _replace_velocities_array(ts, new):
    ts.has_velocities = True
    ts._velocities = new


def _replace_forces_array(ts, new):
    ts.has_forces = True
    ts._forces = new


def _replace_dimensions(ts, new):
    ts._unitcell = new


class MemoryReader(base.ProtoReader):
    """
    MemoryReader works with trajectories represented as numpy arrays.

    A trajectory reader interface to a numpy array of the coordinates.
    For compatibility with the timeseries interface, support is provided for
    specifying the order of columns through the `order` keyword.

    .. versionadded:: 0.16.0
    .. versionchanged:: 1.0.0
       Support for the deprecated `format` keyword for
       :meth:`MemoryReader.timeseries` has now been removed.
    """

    format = 'MEMORY'

    def __init__(self, coordinate_array, order='fac',
                 dimensions=None, dt=1, filename=None,
                 velocities=None, forces=None,
                 **kwargs):
        """
        Parameters
        ----------
        coordinate_array : numpy.ndarray
            The underlying array of coordinates. The MemoryReader now
            necessarily requires a np.ndarray
        order : {"afc", "acf", "caf", "fac", "fca", "cfa"} (optional)
            the order/shape of the return data array, corresponding
            to (a)tom, (f)rame, (c)oordinates all six combinations
            of 'a', 'f', 'c' are allowed ie "fac" - return array
            where the shape is (frame, number of atoms,
            coordinates).
        dimensions: [A, B, C, alpha, beta, gamma] (optional)
            unitcell dimensions (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)
            lengths *A*, *B*, *C* are in the MDAnalysis length unit (Å), and
            angles are in degrees. An array of dimensions can be given,
            which must then be shape (nframes, 6)
        dt: float (optional)
            The time difference between frames (ps).  If :attr:`time`
            is set, then `dt` will be ignored.
        filename: string (optional)
            The name of the file from which this instance is created. Set to ``None``
            when created from an array
        velocities : numpy.ndarray (optional)
            Atom velocities.  Must match shape of coordinate_array.  Will share order
            with coordinates.
        forces : numpy.ndarray (optional)
            Atom forces.  Must match shape of coordinate_array  Will share order
            with coordinates

        Raises
        ------
        TypeError if the coordinate array passed is not a np.ndarray

        Note
        ----
        At the moment, only a fixed `dimension` is supported, i.e., the same
        unit cell for all frames in `coordinate_array`. See issue `#1041`_.


        .. _`#1041`: https://github.com/MDAnalysis/mdanalysis/issues/1041

        .. versionchanged:: 0.19.0
            The input to the MemoryReader now must be a np.ndarray
            Added optional velocities and forces
        .. versionchanged:: 2.2.0
            Input kwargs are now stored under the :attr:`_kwargs` attribute,
            and are passed on class creation in :meth:`copy`.
        """

        super(MemoryReader, self).__init__()
        self.filename = filename
        self.stored_order = order
        self._kwargs = kwargs

        # See Issue #1685. The block below checks if the coordinate array
        # passed is of shape (N, 3) and if it is, the coordiante array is
        # reshaped to (1, N, 3)
        try:
            if coordinate_array.ndim == 2 and coordinate_array.shape[1] == 3:
                coordinate_array = coordinate_array[np.newaxis, :, :]
        except AttributeError:
            errmsg = ("The input has to be a numpy.ndarray that corresponds "
                      "to the layout specified by the 'order' keyword.")
            raise TypeError(errmsg) from None

        self.set_array(coordinate_array, order)
        self.n_frames = \
            self.coordinate_array.shape[self.stored_order.find('f')]
        self.n_atoms = \
            self.coordinate_array.shape[self.stored_order.find('a')]

        if velocities is not None:
            try:
                velocities = np.asarray(velocities, dtype=np.float32)
            except ValueError:
                errmsg = (f"'velocities' must be array-like got "
                          f"{type(velocities)}")
                raise TypeError(errmsg) from None
            # if single frame, make into array of 1 frame
            if velocities.ndim == 2:
                velocities = velocities[np.newaxis, :, :]
            if not velocities.shape == self.coordinate_array.shape:
                raise ValueError('Velocities has wrong shape {} '
                                 'to match coordinates {}'
                                 ''.format(velocities.shape,
                                           self.coordinate_array.shape))
            self.velocity_array = velocities.astype(np.float32, copy=False)
        else:
            self.velocity_array = None

        if forces is not None:
            try:
                forces = np.asarray(forces, dtype=np.float32)
            except ValueError:
                errmsg = f"'forces' must be array like got {type(forces)}"
                raise TypeError(errmsg) from None
            if forces.ndim == 2:
                forces = forces[np.newaxis, :, :]
            if not forces.shape == self.coordinate_array.shape:
                raise ValueError('Forces has wrong shape {} '
                                 'to match coordinates {}'
                                 ''.format(forces.shape,
                                           self.coordinate_array.shape))
            self.force_array = forces.astype(np.float32, copy=False)
        else:
            self.force_array = None

        provided_n_atoms = kwargs.pop("n_atoms", None)
        if (provided_n_atoms is not None and
            provided_n_atoms != self.n_atoms
        ):
            raise ValueError(
                "The provided value for n_atoms ({}) "
                "does not match the shape of the coordinate "
                "array ({})".format(provided_n_atoms, self.n_atoms)
            )

        self.ts = self._Timestep(self.n_atoms, **kwargs)
        self.ts.dt = dt

        if dimensions is None:
            self.dimensions_array = np.zeros((self.n_frames, 6), dtype=np.float32)
        else:
            try:
                dimensions = np.asarray(dimensions, dtype=np.float32)
            except ValueError:
                errmsg = (f"'dimensions' must be array-like got "
                          f"{type(dimensions)}")
                raise TypeError(errmsg) from None
            if dimensions.shape == (6,):
                # single box, tile this to trajectory length
                # allows modifying the box of some frames
                dimensions = np.tile(dimensions, (self.n_frames, 1))
            elif dimensions.shape != (self.n_frames, 6):
                raise ValueError("Provided dimensions array has shape {}. "
                                 "This must be a array of shape (6,) or "
                                 "(n_frames, 6)".format(dimensions.shape))
            self.dimensions_array = dimensions

        self.ts.frame = -1
        self.ts.time = -1
        self._read_next_timestep()

    @staticmethod
    def _format_hint(thing):
        """For internal use: Check if MemoryReader can operate on *thing*

        .. versionadded:: 1.0.0
        """
        return isinstance(thing, np.ndarray)

    @staticmethod
    def parse_n_atoms(filename, order='fac', **kwargs):
        """Deduce number of atoms in a given array of coordinates

        Parameters
        ----------
        filename : numpy.ndarray
          data which will be used later in MemoryReader
        order : {"afc", "acf", "caf", "fac", "fca", "cfa"} (optional)
            the order/shape of the return data array, corresponding
            to (a)tom, (f)rame, (c)oordinates all six combinations
            of 'a', 'f', 'c' are allowed ie "fac" - return array
            where the shape is (frame, number of atoms,
            coordinates).

        Returns
        -------
        n_atoms : int
          number of atoms in system
        """
        # assume filename is a numpy array
        return filename.shape[order.find('a')]

    def copy(self):
        """Return a copy of this Memory Reader"""
        vels = (self.velocity_array.copy()
                if self.velocity_array is not None else None)
        fors = (self.force_array.copy()
                if self.force_array is not None else None)
        dims = self.dimensions_array.copy()

        new = self.__class__(
            self.coordinate_array.copy(),
            order=self.stored_order,
            dimensions=dims,
            velocities=vels,
            forces=fors,
            dt=self.ts.dt,
            filename=self.filename,
            **self._kwargs
        )
        new[self.ts.frame]

        for auxname, auxread in self._auxs.items():
            new.add_auxiliary(auxname, auxread.copy())
        # since transformations are already applied to the whole trajectory
        # simply copy the property
        new.transformations = self.transformations

        return new

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

    def timeseries(self, asel=None, atomgroup=None, start=0, stop=-1, step=1, order='afc'):
        """Return a subset of coordinate data for an AtomGroup in desired
        column order. If no selection is given, it will return a view of the
        underlying array, while a copy is returned otherwise.

        Parameters
        ---------
        asel : AtomGroup (optional)
            Atom selection. Defaults to ``None``, in which case the full set of
            coordinate data is returned. Note that in this case, a view
            of the underlying numpy array is returned, while a copy of the
            data is returned whenever `asel` is different from ``None``.

            .. deprecated:: 2.7.0
               asel argument will be renamed to atomgroup in 3.0.0

        atomgroup: AtomGroup (optional)
            Same as `asel`, will replace `asel` in 3.0.0
        start : int (optional)
            the start trajectory frame
        stop : int (optional)
            the end trajectory frame

            .. deprecated:: 2.4.0
               Note that `stop` is currently *inclusive* but will be
               changed in favour of being *exclusive* in version 3.0.  

        step : int (optional)
            the number of trajectory frames to skip
        order : {"afc", "acf", "caf", "fac", "fca", "cfa"} (optional)
            the order/shape of the return data array, corresponding
            to (a)tom, (f)rame, (c)oordinates all six combinations
            of 'a', 'f', 'c' are allowed ie "fac" - return array
            where the shape is (frame, number of atoms,
            coordinates).


        .. versionchanged:: 1.0.0
           Deprecated `format` keyword has been removed. Use `order` instead.
        .. versionchanged:: 2.4.0
            ValueError now raised instead of NoDataError for empty input
            AtomGroup
        """
        if asel is not None:
            warnings.warn(
                "asel argument to timeseries will be renamed to"
                "'atomgroup' in 3.0, see #3911",
                category=DeprecationWarning)
            if atomgroup:
                raise ValueError("Cannot provide both asel and atomgroup kwargs")
            atomgroup = asel


        if stop != -1:
            warnings.warn("MemoryReader.timeseries inclusive `stop` "
                      "indexing will be removed in 3.0 in favour of exclusive "
                      "indexing", category=DeprecationWarning)

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
        #   1) atomgroup is None
        #   2) atomgroup corresponds to the selection of all atoms.
        array = array[tuple(basic_slice)]
        if (atomgroup is None or atomgroup is atomgroup.universe.atoms):
            return array
        else:
            if len(atomgroup) == 0:
                raise ValueError("Timeseries requires at least one atom "
                                  "to analyze")
            # If selection is specified, return a copy
            return array.take(asel.indices, a_index)

    def _read_next_timestep(self, ts=None):
        """copy next frame into timestep"""
        if ts:
            warnings.warn("ts argument to _read_next_timestep is deprecated as of 2.7.0 and will be removed in 3.0.0, see #3928")

        if self.ts.frame >= self.n_frames-1:
            raise IOError(errno.EIO, 'trying to go over trajectory limit')
        if ts is None:
            ts = self.ts
        ts.frame += 1
        f_index = self.stored_order.find('f')
        basic_slice = ([slice(None)]*(f_index) +
                       [self.ts.frame] +
                       [slice(None)]*(2-f_index))
        _replace_positions_array(ts, self.coordinate_array[tuple(basic_slice)])
        _replace_dimensions(ts, self.dimensions_array[self.ts.frame])
        if self.velocity_array is not None:
            _replace_velocities_array(ts, self.velocity_array[tuple(basic_slice)])
        if self.force_array is not None:
            _replace_forces_array(ts, self.force_array[tuple(basic_slice)])

        ts.time = self.ts.frame * self.dt
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

    def add_transformations(self, *transformations):
        """ Add all transformations to be applied to the trajectory.

        This function take as list of transformations as an argument. These
        transformations are functions that will be called by the Reader and given
        a :class:`Timestep` object as argument, which will be transformed and returned
        to the Reader.
        The transformations can be part of the :mod:`~MDAnalysis.transformations`
        module, or created by the user, and are stored as a list `transformations`.
        This list can only be modified once, and further calls of this function will
        raise an exception.

        .. code-block:: python

          u = MDAnalysis.Universe(topology, coordinates)
          workflow = [some_transform, another_transform, this_transform]
          u.trajectory.add_transformations(*workflow)

        Parameters
        ----------
        transform_list : list
            list of all the transformations that will be applied to the coordinates

        See Also
        --------
        :mod:`MDAnalysis.transformations`
        """
        #Overrides :meth:`~MDAnalysis.coordinates.base.ProtoReader.add_transformations`
        #to avoid unintended behaviour where the coordinates of each frame are transformed
        #multiple times when iterating over the trajectory.
        #In this method, the trajectory is modified all at once and once only.

        super(MemoryReader, self).add_transformations(*transformations)
        for i, ts in enumerate(self):
            for transform in self.transformations:
                ts = transform(ts)

    def _apply_transformations(self, ts):
        """ Applies the transformations to the timestep."""
        # Overrides :meth:`~MDAnalysis.coordinates.base.ProtoReader.add_transformations`
        # to avoid applying the same transformations multiple times on each frame

        return ts
