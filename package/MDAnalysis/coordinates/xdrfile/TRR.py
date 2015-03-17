# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
Gromacs TRR trajectory I/O --- :mod:`MDAnalysis.coordinates.xdrfile.TRR`
========================================================================

Classes for reading and writing of `Gromacs TRR trajectories`_
together with supporting code.

.. note:: Users should access classes from :mod:`MDAnalysis.coordinates.TRR`.

.. _Gromacs TRR trajectories: http://www.gromacs.org/Documentation/File_Formats/.trr_File
.. _Gromacs: http://www.gromacs.org


.. SeeAlso:: :mod:`MDAnalysis.coordinates.xdrfile.libxdrfile2` for low-level
   bindings to the Gromacs trajectory file formats

Classes
-------

.. autoclass:: Timestep
   :members:
   :inherited-members:
.. autoclass:: TRRReader
   :members:
   :inherited-members:
.. autoclass:: TRRWriter
   :members:
   :inherited-members:

"""

import numpy

import core
import libxdrfile2
from MDAnalysis import NoDataError


class Timestep(core.Timestep):
    """Timestep for a Gromacs_ TRR trajectory.

    The Timestep can be initialized with *arg* being

      1. an integer (the number of atoms)
      2. another :class:`Timestep` instance, in which case a copy is made; attention:
         loss of attributes that do not exist within the TRR :class:`Timestep` may occur;
      3. a :class:`numpy.ndarray` of shape ``(numatoms, 3)`` (for positions only) or
         ``(numatoms, 9)`` (for positions, velocities, and forces): ``positions = arg[:,:3]``,
         ``velocities = arg[:,3:6]``, and ``forces = arg[:,6:]``.

    The constructor also takes the named arguments *has_x*, *has_v*, and
    *has_f*, which are used to set the :class:`Timestep` flags
    :attr:`~Timestep.has_x`, :attr:`~Timestep.has_v`, and
    :attr:`~Timestep.has_f`, described below.  Depending on the *arg* use-case
    above, the defaults set for these flags will vary:

      1. when *arg* is an integer :attr:`~Timestep.has_x` defaults to ``True``
         and :attr:`~Timestep.has_v` and :attr:`~Timestep.has_f` to ``False``.

      2. when *arg* is another :class:`Timestep` instance the flags will
         default to being copied from the passed :class:`Timestep`. If that
         instance has no 'has_*' flags the behavior is to assign them to
         ``True`` depending on the existence of :attr:`~Timestep._velocities`
         and :attr:`~Timestep._forces` (:attr:`~Timestep._pos` is assumed to
         always be there, so in this case :attr:`~Timestep.has_x` defaults to
         ``True``).

      3. when *arg* is a numpy array, the default flags will reflect what
         information is passed in the array.

    TRR :class:`Timestep` objects are now fully aware of the existence or not of
    coordinate/velocity/force information in frames, reflected in the
    :attr:`~Timestep.has_x`, :attr:`~Timestep.has_v`, and :attr:`~Timestep.has_f` flags.
    Accessing either kind of information while the corresponding flag is set to ``False``
    wil raise a :exc:`NoDataError`. Internally, however, the arrays are always populated,
    even when the flags are ``False``; upon creation of a :class:`Timestep` they are
    zero-filled, but this might not always be the case later on for properties flagged as
    ``False`` if the same :class:`Timestep` instance is used to read from a TRR frame.

    When doing low-level writing to :attr:`~Timestep._pos`, :attr:`~Timestep._velocities`,
    or :attr:`~Timestep._forces:attr:, the corresponding flags must be set beforehand. The
    TRR :class:`Timestep` constructor allows for the named boolean arguments *has_x*,
    *has_v*, and *has_f* to be passed for automatic setting of the corresponding flag.
    An exception to this is assignment to the full property array thus::

        ts = MDAnalysis.coordinates.TRR.Timestep(N)     # N being the number of atoms
        ts._velocities = vel_array   # Where vel_array is an existing array of shape (N, DIM)
                                     #  This will also automatically set 'has_v' to True.

    Attempting to populate the array instead will, however, raise a NoDataError exception::

        ts = MDAnalysis.coordinates.TRR.Timestep(N)     # N being the number of atoms
        ts._velocities[:] = vel_array   #  This will fail if 'has_v' hasn't been set to True.

    .. versionchanged:: 0.8.0
       TRR :class:`Timestep` objects are now fully aware of the existence or
       not of coordinate/velocity/force information in frames.

    .. _Gromacs: http://www.gromacs.org
    """
    # The exception error
    _nodataerr = "You are accessing the %s of a Timestep but there are none in this frame. \
It might be the case that your trajectory doesn't have %s, or not every frame. In the latter case you can \
(1) get rid of %s-less frames before passing the trajectory to MDAnalysis, or \
(2) reflow your code to not access %s when they're not there, by making use of the '%s' flag of Timestep objects."

    def __init__(self, arg, **kwargs):
        DIM = libxdrfile2.DIM  # compiled-in dimension (most likely 3)
        # These flags control which attributes were actually initialized in this Timestep. Due to the possibility of TRR
        #  frames with no coords/vels/forces, and the way reading is handled, the _tpos, _tvelocities and _tforces
        # arrays
        #  may end up containing data from different frames when the relevant attribute is missing in the trajectory.
        #  This is flagged and exceptions can be raised whenever there is an attempt to read flagged data.
        if numpy.dtype(type(arg)) == numpy.dtype(int):
            self.has_x = kwargs.pop('has_x', True)
            self.has_v = kwargs.pop('has_v', False)
            self.has_f = kwargs.pop('has_f', False)
            self.frame = 0
            self.numatoms = arg
            # C floats and C-order for arrays (see libxdrfile2.i)
            self._tpos = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')
            self._tvelocities = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')
            self._tforces = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')
            self._unitcell = numpy.zeros((DIM, DIM), dtype=numpy.float32)
            # additional data for xtc
            self.status = libxdrfile2.exdrOK
            self.step = 0
            self.time = 0
            self.lmbda = 0
        elif isinstance(arg, Timestep):  # Copy constructor
            # This makes a deepcopy of the timestep
            self.frame = arg.frame
            self.numatoms = arg.numatoms
            self._unitcell = numpy.array(arg._unitcell)
            # The 'has_' flags are set to True if the passed ts has them set to True OR if there is no such flag
            #  and the attributes are there. This is overriden by the has_ named arguments passed to the constructor.
            # has_x is an exception because we can assume that, if the flag doesn't exist, the _pos attribute will
            # always be there.
            self.has_x = kwargs.pop('has_x', arg.__dict__.get("has_x", True))
            self.has_v = kwargs.pop('has_v', arg.__dict__.get("has_v", hasattr(arg, "_velocities")))
            self.has_f = kwargs.pop('has_f', arg.__dict__.get("has_f", hasattr(arg, "_forces")))
            # We now either copy the properties or initialize a zeros array, depending on whether the passed ts
            #  has the appropriate 'has_' flags set or has the attribute. If the attribute is there but the 'has_'
            #  flag is False we initialize to zeros (it's just as slow and it makes no sense to copy data from
            #  undefined behavior).
            #
            # Copies are done using the casted numpy.array function, which serves to early identify invalid input data.
            # COORDINATES:
            if self.has_x:
                try:
                    self._tpos = numpy.array(arg._pos, dtype=numpy.float32)
                except ValueError as err:
                    raise ValueError("Attempted to create a Timestep with invalid coordinate data: " + err.message)
            else:
                self._tpos = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')
            # VELOCITIES:
            if self.has_v:
                try:
                    self._tvelocities = numpy.array(arg._velocities, dtype=numpy.float32)
                except ValueError as err:
                    raise ValueError("Attempted to create a Timestep with invalid velocity data: " + err.message)
                except AttributeError:
                    raise AttributeError(
                        "Attempted to create a new Timestep from inconsistent Timestep data (the passed Timestep "
                        "object has the 'has_v' flag set, but the '_velocities' attribute is missing.)")
            else:
                self._tvelocities = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')
            # FORCES:
            if self.has_f:
                try:
                    self._tforces = numpy.array(arg._forces, dtype=numpy.float32)
                except ValueError as err:
                    raise ValueError("Attempted to create a Timestep with invalid force data: " + err.message)
                except AttributeError:
                    raise AttributeError(
                        "Attempted to create a new Timestep from inconsistent Timestep data (the passed Timestep "
                        "object has the 'has_f' flag set, but the '_forces' attribute is missing.)")
            else:
                self._tforces = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')
            ##
            for attr in ('status', 'step', 'time', 'lmbda', 'prec'):
                if hasattr(arg, attr):
                    self.__setattr__(attr, arg.__getattribute__(attr))
        elif isinstance(arg, numpy.ndarray):
            # provide packed array shape == (natoms, 3*DIM)
            # which contains pos = arg[:,0:3], v = arg[:,3:6], f = arg[:, 6:9]
            # or just positions: pos = arg[:,0:3] == arg
            #
            # The 'has_' flags are set from the named arguments passed to the constructor,
            # but default to what fields are available from the input array.
            self.has_x = kwargs.pop('has_x', True)
            self.has_v = kwargs.pop('has_v', arg.shape[1] == 3 * DIM)
            self.has_f = kwargs.pop('has_f', arg.shape[1] == 3 * DIM)
            if len(arg.shape) != 2:
                raise ValueError("packed numpy array (x,v,f) can only have 2 dimensions")
            self._unitcell = numpy.zeros((DIM, DIM), dtype=numpy.float32)
            self.frame = 0
            if (arg.shape[0] == 3 * DIM and arg.shape[1] != 3 * DIM) or \
                    (arg.shape[0] == DIM and arg.shape[1] != DIM):
                # wrong order (but need to exclude case where natoms == DIM or natoms == 3*DIM!)
                raise ValueError("TRR timestep is to be initialized from an array with dimensions (natoms, "
                                 "3*%d) or (natoms, %d). \
If you wish to skip a property set the corresponding 'has_' flag to False on construction." % (DIM, DIM))
            self.numatoms = arg.shape[0]
            if arg.shape[1] != DIM and arg.shape[1] != 3 * DIM:
                raise ValueError(
                    "TRR timestep doesn't have second dimension %d or 3*%d: shape=%r" % (DIM, DIM, arg.shape))
            # COORDINATES
            self._tpos = arg[:, 0:DIM].copy('C')  # C-order
            # Velocities and forces, if the array seems to have them.
            # VELOCITIES
            if arg.shape[1] == DIM:
                self._tvelocities = arg[:, DIM:2 * DIM].copy('C')  # C-order
                # FORCES
                self._tforces = arg[:, 2 * DIM:3 * DIM].copy('C')  # C-order
            else:
                self._tvelocities = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')
                self._tforces = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')

            # additional data for trr
            self.status = libxdrfile2.exdrOK
            self.step = 0
            self.time = 0
            self.lmbda = 0
        else:
            raise ValueError("Cannot create an empty Timestep")

        # The entire coordinate/velocity/force information is now hidden behind
        #  existence-checking decorated functions.

    # COORDINATES
    @property
    def _pos(self):
        if self.has_x:
            return self._tpos
        else:
            raise NoDataError(self._nodataerr % ("coordinates", "coordinates", "coordinate", "coordinates", "has_x"))

    @_pos.setter
    def _pos(self, x):
        if x.shape == (self.numatoms, libxdrfile2.DIM):
            self._tpos = x
            self.has_x = True
        else:
            raise ValueError("You are attempting to set the positions array of a Timestep with an array \
that doesn't have the same number of atoms or the same number of dimensions. The Timestep \
has number-of-atoms,dimensions %r, and you supplied an array of shape %r." %
                             ((self.numatoms, libxdrfile2.DIM), x.shape))

    #
    @property
    def _x(self):
        return self._pos[:, 0]

    @_x.setter
    def _x(self, x):
        self._tpos[:, 0] = x

    @property
    def _y(self):
        return self._pos[:, 1]

    @_y.setter
    def _y(self, y):
        self._tpos[:, 1] = y

    @property
    def _z(self):
        return self._pos[:, 2]

    @_z.setter
    def _z(self, z):
        self._tpos[:, 2] = z

    # VELOCITIES
    @property
    def _velocities(self):
        if self.has_v:
            return self._tvelocities
        else:
            raise NoDataError(self._nodataerr % ("velocities", "velocities", "velocity", "velocities", "has_v"))

    @_velocities.setter
    def _velocities(self, v):
        if v.shape == (self.numatoms, libxdrfile2.DIM):
            self._tvelocities = v
            self.has_v = True
        else:
            raise ValueError("You are attempting to set the velocities array of a Timestep with an array \
that doesn't have the same number of atoms or the same number of dimensions. The Timestep \
has number-of-atoms,dimensions %r, and you supplied an array of shape %r." %
                             ((self.numatoms, libxdrfile2.DIM), v.shape))

    #FORCES
    @property
    def _forces(self):
        if self.has_f:
            return self._tforces
        else:
            raise NoDataError(self._nodataerr % ("forces", "forces", "force", "forces", "has_f"))

    @_forces.setter
    def _forces(self, f):
        if f.shape == (self.numatoms, libxdrfile2.DIM):
            self._tforces = f
            self.has_f = True
        else:
            raise ValueError("You are attempting to set the forces array of a Timestep with an array \
that doesn't have the same number of atoms or the same number of dimensions. The Timestep \
has number-of-atoms,dimensions %r, and you supplied an array of shape %r." %
                             ((self.numatoms, libxdrfile2.DIM), f.shape))


class TRRWriter(core.TrjWriter):
    """Write a Gromacs_ TRR trajectory.

    .. _Gromacs: http://www.gromacs.org
    """
    format = "TRR"
    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps', 'force': 'kJ/(mol*nm)'}


class TRRReader(core.TrjReader):
    """Read a Gromacs_ TRR trajectory.

    .. versionchanged:: 0.8.0
       :class:`Timestep` objects returned from TRR files now have
       :attr:`~Timestep.has_x`, :attr:`~Timestep.has_v`, and :attr:`~Timestep.has_f`
       flags reflecting whether coordinates/velocities/forces were read.
       Attempting to access such data when the corresponding flag is set to ``False``
       will raise a :exc:`NoDataError`.

    .. _Gromacs: http://www.gromacs.org
    """
    format = "TRR"
    _Timestep = Timestep
    _Writer = TRRWriter
    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps', 'force': 'kJ/(mol*nm)'}
