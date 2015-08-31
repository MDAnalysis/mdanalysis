# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
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
Common high-level Gromacs XDR functionality --- :mod:`MDAnalysis.coordinates.xdrfile.core`
==========================================================================================

The :mod:`MDAnalysis.coordinates.xdrfile.core` module contains generic
classes to access Gromacs_ XDR-encoded trajectory formats such as TRR
and XTC.

A generic Gromacs_ trajectory is simply called "trj" within this
module.

.. SeeAlso:: :mod:`MDAnalysis.coordinates.base` for the generic MDAnalysis base
             classes and :mod:`MDAnalysis.coordinates.xdrfile.libxdrfile2` for
             the low-level bindings to the XDR trajectories.

.. _Gromacs: http://www.gromacs.org

Generic xdr trj classes
-----------------------

The generic classes are subclassed to generate the specific classes
for the XTC and TRR format.

.. versionchanged:: 0.8.0
   The XTC/TRR I/O interface now uses
   :mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile2`, which has seeking and
   indexing capabilities. Note that unlike
   :mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile` before it,
   :mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile2` is distributed under the
   GNU GENERAL PUBLIC LICENSE, version 2 (or higher).

.. versionchanged:: 0.9.0
   TrjReader now stores the offsets used for frame seeking automatically as a
   hidden file in the same directory as the source trajectory. These offsets
   are automatically retrieved upon TrjReader instantiation, resulting in
   substantially quicker initialization times for long trajectories. The
   ctime and filesize of the trajectory are stored with the offsets, and these
   are checked against the trajectory on load to ensure the offsets aren't
   stale. The offsets are automatically regenerated if they are stale or
   missing.

.. versionchanged:: 0.11.0
   Frames now 0-based instead of 1-based

.. autoclass:: Timestep
   :members:
.. autoclass:: TrjReader
   :members:
.. autoclass:: TrjWriter
   :members:

"""

import os
import errno
import numpy as np
import sys
import cPickle
import warnings
import weakref

from . import libxdrfile2
from MDAnalysis.coordinates import base
from MDAnalysis.coordinates.core import triclinic_box, triclinic_vectors
import MDAnalysis.core
from ...lib.util import cached


# This is the XTC class. The TRR overrides with it's own.
class Timestep(base.Timestep):
    """Timestep for a Gromacs trajectory.

    .. versionchanged:: 0.11.0
       Attributes status, lmbda, prec all stored in the :attr:`data` dictionary
       Native frame number now stored as `_frame`, was `step`
    """
    order = 'C'

    def __init__(self, n_atoms, **kwargs):
        super(Timestep, self).__init__(n_atoms, **kwargs)
        self.data['status'] = libxdrfile2.exdrOK
        self._frame = 0
        self.data['prec'] = 0

    def _init_unitcell(self):
        return np.zeros((3, 3), dtype=np.float32)

    @property
    def dimensions(self):
        """unitcell dimensions (A, B, C, alpha, beta, gamma)

        - A, B, C are the lengths of the primitive cell vectors e1, e2, e3
        - alpha = angle(e1, e2)
        - beta = angle(e1, e3)
        - gamma = angle(e2, e3)
        """
        # Layout of unitcell is [X, Y, Z] with the primitive cell vectors
        x = self._unitcell[0]
        y = self._unitcell[1]
        z = self._unitcell[2]
        return triclinic_box(x, y, z)

    @dimensions.setter
    def dimensions(self, box):
        self._unitcell = triclinic_vectors(box)


class TrjWriter(base.Writer):
    """Writes to a Gromacs trajectory file

    (Base class)
    """
    #: units of time (ps) and length (nm) in Gromacs
    units = {'time': 'ps', 'length': 'nm'}
    #: override to define trajectory format of the reader (XTC or TRR)
    format = None

    def __init__(self, filename, n_atoms, start=0, step=1, dt=None, precision=1000.0, remarks=None,
                 convert_units=None):
        """ Create a new TrjWriter

        :Arguments:
          *filename*
             name of output file
          *n_atoms*
             number of atoms in trajectory file

        :Keywords:
          *start*
             starting timestep frame; only used when *dt* is set.
          *step*
             skip in frames between subsequent timesteps; only used when *dt* is set.
          *dt*
             time between frames to use. If set will override any time information
             contained in the passed :class:`Timestep` objects, which will otherwise
             be used.
             The :attr:`~Timestep.time` attribute defaults to a timestep of
             to setting the trajectory time at 1 ps per step if there is no
             time information.
          *precision*
              accuracy for lossy XTC format as a power of 10 (ignored
              for TRR) [1000.0]
          *convert_units*
             ``True``: units are converted to the MDAnalysis base format; ``None`` selects
             the value of :data:`MDAnalysis.core.flags` ['convert_lengths'].
             (see :ref:`flags-label`)

        .. versionchanged:: 0.8.0
           The TRR writer is now able to write TRRs without coordinates/velocities/forces,
           depending on the properties available in the :class:`Timestep` objects passed to
           :meth:`~TRRWriter.write`.
        .. versionchanged:: 0.11.0
           Keyword "delta" renamed to "dt"
        """
        if n_atoms == 0:
            raise ValueError("TrjWriter: no atoms in output trajectory")
        self.filename = filename
        # Convert filename to ascii because of SWIG bug.
        # See: http://sourceforge.net/p/swig/feature-requests/75
        # Only needed for Python < 3
        if sys.version_info[0] < 3:
            if isinstance(filename, unicode):
                self.filename = filename.encode("UTF-8")

        if convert_units is None:
            convert_units = MDAnalysis.core.flags['convert_lengths']
        self.convert_units = convert_units  # convert length and time to base units on the fly?
        self.n_atoms = n_atoms

        self.frames_written = 0
        self.start = start
        self.step = step
        self.dt = dt
        self.remarks = remarks
        self.precision = precision  # only for XTC
        self.xdrfile = libxdrfile2.xdrfile_open(self.filename, 'w')

        self.ts = None
        # To flag empty properties to be skipped when writing a TRR it suffices to pass an empty 2D array with shape(
        # natoms,0)
        if self.format == 'TRR':
            self._emptyarr = np.array([], dtype=np.float32).reshape(self.n_atoms, 0)

    def write_next_timestep(self, ts=None):
        """ write a new timestep to the trj file

        *ts* is a :class:`Timestep` instance containing coordinates to
        be written to trajectory file
        """
        if self.xdrfile is None:
            raise IOError("Attempted to write to closed file %r", self.filename)
        if ts is None:
            if not hasattr(self, "ts"):
                raise IOError("TrjWriter: no coordinate data to write to trajectory file")
            else:
                ts = self.ts
        elif not ts.n_atoms == self.n_atoms:
            # Check to make sure Timestep has the correct number of atoms
            raise IOError("TrjWriter: Timestep does not have the correct number of atoms")

        status = self._write_next_timestep(ts)

        if status != libxdrfile2.exdrOK:
            raise IOError(errno.EIO, "Error writing %s file (status %d)" % (self.format, status), self.filename)
        self.frames_written += 1

    def _write_next_timestep(self, ts):
        """Generic writer for XTC and TRR with minimum intelligence; override if necessary."""

        # (1) data common to XTC and TRR

        # Time-writing logic: if the writer was created with a dt parameter,
        #  use dt*(start+step*frames_written)
        #  otherwise use the provided Timestep obj time attribute
        if self.dt is None:
            time = ts.time
        else:
            time = (self.start + self.step * self.frames_written) * self.dt

        if self.convert_units:
            time = self.convert_time_to_native(time, inplace=False)

        try:
            step = int(ts._frame)
        except AttributeError:
            # bogus, should be actual MD step number, i.e. frame * dt/delta
            step = ts.frame

        unitcell = self.convert_dimensions_to_unitcell(ts).astype(np.float32)  # must be float32 (!)

        # make a copy of the scaled positions so that the in-memory
        # timestep is not changed (would have lead to wrong results if
        # analysed *after* writing a time step to disk). The new
        # implementation could lead to memory problems and/or slow-down for
        # very big systems because we temporarily create a new array pos
        # for each frame written
        #
        # For TRR only go through the trouble if the frame actually has valid
        #  coords/vels/forces; otherwise they won't be written anyway (pointers
        #  set to an empty array that libxdrfile2.py knows it should set to NULL).
        #
        # (2) have to treat XTC and TRR somewhat differently
        if self.format == 'XTC':
            if self.convert_units:
                pos = self.convert_pos_to_native(ts._pos, inplace=False)
            else:
                pos = ts._pos
            status = libxdrfile2.write_xtc(self.xdrfile, step, float(time), unitcell, pos, self.precision)
        elif self.format == 'TRR':
            try:
                lmbda = ts.data['lmbda']
            except KeyError:
                lmbda = 1.0

            # COORDINATES
            if ts.has_positions:
                if self.convert_units:
                    pos = self.convert_pos_to_native(ts._pos, inplace=False)
                else:
                    pos = ts._pos
            else:
                pos = self._emptyarr
            #VELOCITIES
            if ts.has_velocities:
                if self.convert_units:
                    velocities = self.convert_velocities_to_native(ts._velocities, inplace=False)
                else:
                    velocities = ts._velocities
            else:
                velocities = self._emptyarr
            # FORCES
            if ts.has_forces:
                if self.convert_units:
                    forces = self.convert_forces_to_native(ts._forces, inplace=False)
                else:
                    forces = ts._forces
            else:
                forces = self._emptyarr

            status = libxdrfile2.write_trr(self.xdrfile, step, float(time), lmbda, unitcell,
                                           pos, velocities, forces)

        return status

    def close(self):
        status = libxdrfile2.exdrCLOSE
        if not self.xdrfile is None:
            status = libxdrfile2.xdrfile_close(self.xdrfile)
            self.xdrfile = None
        return status

    def convert_dimensions_to_unitcell(self, ts):
        """Read dimensions from timestep *ts* and return Gromacs box vectors"""
        return self.convert_pos_to_native(triclinic_vectors(ts.dimensions))


class TrjReader(base.Reader):
    """Generic base class for reading Gromacs trajectories inside MDAnalysis.

    Derive classes and set :attr:`TrjReader.format`, :attr:`TrjReader._read_trj` and :attr:`TrjReader._read_trj_atoms`.

    Example::
       reader = TrjReader("file.trj")
       for ts in reader:
          print ts

    """
    #: units of time (ps) and length (nm) in Gromacs
    units = {'time': 'ps', 'length': 'nm'}
    #: override to define trajectory format of the reader (XTC or TRR)
    format = None
    #: supply the appropriate Timestep class, e.g.
    #: :class:`MDAnalysis.coordinates.xdrfile.XTC.Timestep` for XTC
    _Timestep = Timestep
    #: writer class that matches this reader (override appropriately)
    _Writer = TrjWriter

    def __init__(self, filename, sub=None, **kwargs):
        """
        :Arguments:
            *filename*
                the name of the trr file.

        :Keywords:
            *sub*
                an numpy integer array of what subset of trajectory atoms to load into
                the timestep. Intended to work similarly to the 'sub' argument to Gromacs_' trjconv.

                This is usefull when one has a Universe loaded with only an unsolvated protein, and
                wants to read a solvated trajectory.

                The length of this array must be <= to the actual number of atoms in the trajectory, and
                equal to number of atoms in the Universe.
            *refresh_offsets*
                if ``True``, do not retrieve stored offsets, but instead generate
                new ones; if ``False``, use retrieved offsets if available [``False``]

        .. versionchanged:: 0.9.0
           New keyword *refresh_offsets*
        .. versionchanged:: 0.11.0
           Renamed "delta" attribute to "dt"
           Now passes weakref of self to ts (as "_reader")
        """
        super(TrjReader, self).__init__(filename, **kwargs)
        self._cache = dict()
        # Convert filename to ascii because of SWIG bug.
        # See: http://sourceforge.net/p/swig/feature-requests/75
        # Only needed for Python < 3
        if sys.version_info[0] < 3:
            if isinstance(filename, unicode):
                self.filename = filename.encode("UTF-8")

        self.xdrfile = None

        self._n_frames = None  # takes a long time, avoid accessing self.n_frames
        self._dt = None  # compute from time in first two frames!
        self._offsets = None  # storage of offsets in the file

        # actual number of atoms in the trr file
        # first time file is opened, exception should be thrown if bad file
        self._trr_n_atoms = self._read_trj_natoms(self.filename)

        # logic for handling sub sections of trr:
        # this class has some tmp buffers into which the libxdrfile2 functions read the
        # entire frame, then this class copies the relevant sub section into the timestep.
        # the sub section logic is contained entierly within this class.

        # check to make sure sub is valid (if given)
        if sub is not None:
            # only valid type
            if not isinstance(sub, np.ndarray) or len(sub.shape) != 1 or sub.dtype.kind != 'i':
                raise TypeError("sub MUST be a single dimensional numpy array of integers")
            if len(sub) > self._trr_n_atoms:
                raise ValueError("sub MUST be less than or equal to the number of actual trr atoms,"
                                 " {0} in this case".format(self._trr_n_atoms))
            if np.max(sub) >= self._trr_n_atoms or np.min(sub) < 0:
                raise IndexError("sub contains out-of-range elements for the given trajectory")
            # sub appears to be valid
            self._sub = sub

            # make tmp buffers
            # C floats and C-order for arrays (see libxdrfile2.i)
            DIM = libxdrfile2.DIM  # compiled-in dimension (most likely 3)
            # XTC and TRR allocate different things, so call this
            self._allocate_sub(DIM)
        else:
            self._sub = None
            self._pos_buf = None
            self._velocities_buf = None
            self._forces_buf = None

        # make the timestep, this is ALWAYS the used the public number of atoms
        # (same as the calling Universe)
        # at this time, _trr_n_atoms and _sub are set, so self.n_atoms has all it needs
        # to determine number of atoms.
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self.ts._reader = weakref.ref(self)

        # Read in the first timestep
        self._read_next_timestep()

        # try retrieving stored offsets
        if not kwargs.pop('refresh_offsets', False):
            self._retrieve_offsets()

    @property
    def n_atoms(self):
        """The number of publically available atoms that this reader will store in the timestep.

        If 'sub' was not given in the ctor, then this value will just be the actual
        number of atoms in the underlying trajectory file. If however 'sub' was given, then
        this value is the number specified by the 'sub' sub-selection.

        If for any reason the trajectory cannot be read then a negative value is returned.
        """
        return len(self._sub) if self._sub is not None else self._trr_n_atoms

    @property
    def n_frames(self):
        """Read the number of frames from the trajectory.

        The result is cached. If for any reason the trajectory cannot
        be read then 0 is returned.

        This takes a long time because the frames are counted by
        iterating through the whole trajectory. If the trajectory
        was previously loaded and saved offsets exist, then
        loading will be significantly faster.

        .. SeeAlso:: :meth:`TrjReader.load_offsets` and :meth:`TrjReader.save_offsets`
        """
        if not self._n_frames is None:  # return cached value
            return self._n_frames
        try:
            self._read_trj_n_frames(self.filename)
        except IOError:
            self._n_frames = 0
            return 0
        else:
            return self._n_frames

    @property
    def offsets(self):
        if self._offsets is not None:
            return self._offsets
        try:
            self._read_trj_n_frames(self.filename)
        except IOError:
            self._offsets = []
            return 0
        else:
            return self._offsets

    def _get_dt(self):
        """Time step length in ps.

        The result is computed from the trajectory and cached. If for
        any reason the trajectory cannot be read then 0 is returned.
        """
        curr = self.ts.frame
        # no need for conversion: it's alread in our base unit ps
        try:
            t0 = self.ts.time
            self.next()
            t1 = self.ts.time
            dt = t1 - t0
            return dt
        except IOError:
            return 0
        finally:
            self[curr]

    def _offset_filename(self):
        head, tail = os.path.split(self.filename)
        return os.path.join(head, '.{0}_offsets.pkl'.format(tail))

    def _store_offsets(self):
        """Stores offsets for trajectory as a hidden file in the same directory
            as the trajectory itself.

        .. versionadded: 0.9.0
        """
        # try to store offsets; if fails (due perhaps to permissions), then
        # don't bother
        try:
            self.save_offsets(self._offset_filename())
        except IOError:
            warnings.warn("Offsets could not be stored; they will rebuilt when needed next.")

    def _retrieve_offsets(self):
        """Attempts to retrieve previously autosaved offsets for trajectory.

        .. versionadded: 0.9.0
        """
        try:
            self.load_offsets(self._offset_filename(), check=True)
        except IOError:
            warnings.warn("Offsets could not be retrieved; they will be rebuilt instead.")

    def save_offsets(self, filename):
        """Saves current trajectory offsets into *filename*, as a pickled object.

        Along with the offsets themselves, the ctime and file size of the
        trajectory file are also saved.  These are used upon load as a check to
        ensure the offsets still match the trajectory they are being applied
        to.

        The offset file is a pickled dictionary with keys/values::
          *ctime*
             the ctime of the trajectory file
          *size*
             the size of the trajectory file
          *offsets*
             a numpy array of the offsets themselves

        :Arguments:
          *filename*
              filename in which to save the frame offsets

        .. versionadded: 0.8.0
        .. versionchanged: 0.9.0
           Format of the offsets file has changed. It is no longer a pickled numpy
           array, but now a pickled dictionary. See details above. Old offset files
           can no longer be loaded.
            
        """
        if self._offsets is None:
            self._read_trj_n_frames(self.filename)

        output = {'ctime': os.path.getctime(self.filename),
                  'size': os.path.getsize(self.filename),
                  'offsets': self._offsets}

        with open(filename, 'wb') as f:
            cPickle.dump(output, f)

    def load_offsets(self, filename, check=False):
        """Loads current trajectory offsets from pickled *filename*. 

        Checks if ctime and size of trajectory file matches that stored in
        pickled *filename*.  If either one does not match (and *check* == ``True``)
        then the offsets are not loaded.  This is intended to conservatively
        avoid loading out-of-date offsets.

        The offset file is expected to be a pickled dictionary with keys/values::
          *ctime*
             the ctime of the trajectory file
          *size*
             the size of the trajectory file
          *offsets*
             a numpy array of the offsets themselves

        :Arguments:
          *filename*
              filename of pickle file saved with :meth:`~TrjReader.save_offsets` 
              with the frame offsets for the loaded trajectory

        :Keywords:
          *check*
              if False, ignore ctime and size check of trajectory file
        
        :Raises: :exc:`IOError` if the file cannot be read (see :func:`open`).

        .. versionadded: 0.8.0
        .. versionchanged: 0.9.0
           Format of the offsets file has changed. It is no longer a pickled numpy
           array, but now a pickled dictionary. See details above. Old offset files
           can no longer be loaded.

        """
        if not os.path.isfile(filename):  # Return silently if the offset file is not present
            return

        with open(filename, 'rb') as f:
            offsets = cPickle.load(f)

        if check:
            conditions = False
            try:
                ## ensure all conditions are met
                # ctime of file must match that stored
                key = 'ctime'
                conditions = (os.path.getctime(self.filename) == offsets[key])

                # file size must also match
                key = 'size'
                conditions = (os.path.getsize(self.filename) == offsets[key]) and conditions
            except KeyError:
                warnings.warn("Offsets in file '{0}' not suitable;"
                              " missing {1}.".format(filename, key))
                return

            # if conditions not met, abort immediately
            if not conditions:
                warnings.warn("Aborted loading offsets from file; ctime or size did not match.")
                return

        # try to load offsets
        try:
            self._offsets = offsets['offsets']
        except KeyError:
            warnings.warn("Missing key 'offsets' in file '{0}';"
                          " aborting load of offsets.".format(filename))
            return
        self._n_frames = len(self._offsets)

        # finally, check that loaded offsets appear to work by trying
        # to load last frame; otherwise, dump them so they get regenerated
        # on next call to ``self.n_frames``

        #store current frame
        frame = self.frame
        try:
            self.__getitem__(-1)
            # ensure we return to the frame we started with
            self.__getitem__(frame)
        except (IndexError, IOError):
            warnings.warn("Could not access last frame with loaded offsets;"
                          " will rebuild offsets instead.")
            self._offsets = None
            self._n_frames = None

    def open_trajectory(self):
        """Open xdr trajectory file.

        :Returns: pointer to XDRFILE (and sets self.xdrfile)
        :Raises:  :exc:`IOError` with code EALREADY if file was already opened or
                  ENOENT if the file cannot be found
        """
        if not self.xdrfile is None:
            raise IOError(errno.EALREADY, 'XDR file already opened', self.filename)
        if not os.path.exists(self.filename):
            # must check; otherwise might segmentation fault
            raise IOError(errno.ENOENT, 'XDR file not found', self.filename)
        self.xdrfile = libxdrfile2.xdrfile_open(self.filename, 'r')
        # reset ts
        ts = self.ts
        ts.data['status'] = libxdrfile2.exdrOK
        ts.frame = -1 
        ts._frame = 0

        # additional data for XTC
        ts.data['prec'] = 0
        # additional data for TRR
        ts.data['lmbda'] = 0
        return self.xdrfile

    def close(self):
        """Close xdr trajectory file if it was open."""
        if self.xdrfile is None:
            return
        libxdrfile2.xdrfile_close(self.xdrfile)
        self.xdrfile = None  # guard against  crashing with a double-free pointer

    def Writer(self, filename, **kwargs):
        """Returns a Gromacs TrjWriter for *filename* with the same parameters as this trajectory.

        All values can be changed through keyword arguments.

        :Arguments:
          *filename*
              filename of the output trajectory

        :Keywords:
          *n_atoms*
              number of atoms
          *dt*
              Time interval between frames.
          *precision*
              accuracy for lossy XTC format as a power of 10 (ignored
              for TRR) [1000.0]

        :Returns: appropriate :class:`TrjWriter`

        .. versionchanged:: 0.11.0
           Changed "delta" keyword to "dt"
        """
        n_atoms = kwargs.pop('n_atoms', self.n_atoms)

        kwargs.setdefault('dt', self.dt)
        try:
            kwargs['start'] = self[0].time / kwargs['dt']
        except (AttributeError, ZeroDivisionError):
            kwargs['start'] = 0
        try:
            kwargs.setdefault('precision', self.precision)
        except AttributeError:
            pass  # not needed for TRR
        return self._Writer(filename, n_atoms, **kwargs)

    def __iter__(self):
        self._reopen()
        while True:
            try:
                ts = self._read_next_timestep()
            except IOError as err:
                if err.errno == errno.EIO:
                    break
                else:
                    self.close()
                    raise
            except:
                self.close()
                raise
            else:
                yield ts

    def _reopen(self):
        self.close()
        self.open_trajectory()

    def _read_frame(self, frame):
        """ Fast, index-based, random frame access

        """
        if self._offsets is None:
            self._read_trj_n_frames(self.filename)
        self._seek(self._offsets[frame])
        self.ts.frame = frame - 1 # frame gets +1'd in _read_next_timestep   
        self._read_next_timestep()
        return self.ts

    # Renamed this once upon a time.
    _goto_frame = _read_frame

    def timeseries(self, asel, start=0, stop=-1, skip=1, format='afc'):
        raise NotImplementedError("timeseries not available for Gromacs trajectories")

    def correl(self, timeseries, start=0, stop=-1, skip=1):
        raise NotImplementedError("correl not available for Gromacs trajectories")

    def _seek(self, pos, rel=False):
        """Traj seeker"""
        if rel:
            status = libxdrfile2.xdr_seek(self.xdrfile, long(pos), libxdrfile2.SEEK_CUR)
        else:
            status = libxdrfile2.xdr_seek(self.xdrfile, long(pos), libxdrfile2.SEEK_SET)
        if status != libxdrfile2.exdrOK:
            raise IOError(errno.EIO,
                          "Problem seeking to offset %d (relative to current = %s) on file %s, status %s.\n"
                          "Perhaps you are trying to read a file >2GB and your system does not have large file"
                          "support?" % (pos, rel, self.filename, status))

    def _tell(self):
        """Traj pos getter"""
        return libxdrfile2.xdr_tell(self.xdrfile)

