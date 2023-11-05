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
# TRZ Reader written by Richard J. Gowers (2013)

"""TRZ trajectory I/O  --- :mod:`MDAnalysis.coordinates.TRZ`
============================================================

Classes to read `IBIsCO`_ / `YASP`_ TRZ binary trajectories, including
coordinates, velocities and more (see attributes of the :class:`Timestep`).

Data are read and written in binary representation but because this depends on
the machine hardware architecture, MDAnalysis *always* reads and writes TRZ
trajectories in *little-endian* byte order.

.. _IBIsCO: http://www.theo.chemie.tu-darmstadt.de/ibisco/IBISCO.html
.. _YASP: http://www.theo.chemie.tu-darmstadt.de/group/services/yaspdoc/yaspdoc.html

Classes
-------

.. autoclass:: TRZReader
   :members:

.. autoclass:: TRZWriter
   :members:

"""
import warnings
import numpy as np
import os
import errno

from . import base
from ..exceptions import NoDataError
from .timestep import Timestep
from ..lib import util
from ..lib.util import cached, store_init_arguments
from .core import triclinic_box, triclinic_vectors


class TRZReader(base.ReaderBase):
    """Reads an IBIsCO or YASP trajectory file

    Attributes
    ----------
    ts : timestep.Timestep
         :class:`~MDAnalysis.coordinates.timestep.Timestep` object containing
         coordinates of current frame

    Note
    ----
    Binary TRZ trajectories are *always* assumed to be written in
    *little-endian* byte order and are read as such.


    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based.
       Extra data (Temperature, Energies, Pressures, etc) now read
       into ts.data dictionary.
       Now passes a weakref of self to ts (ts._reader).
    .. versionchanged:: 1.0.1
       Now checks for the correct `n_atoms` on reading
       and can raise :exc:`ValueError`.
    .. versionchanged:: 2.1.0
       TRZReader now returns a default :attr:`dt` of 1.0 when it cannot be
       obtained from the difference between two frames.
    .. versionchanged:: 2.3.0
       _frame attribute moved to `ts.data` dictionary.
    .. deprecated:: 2.7.0
       The TRZ Reader and Writer are deprecated as of version 2.7.0
       and will be removed in version 3.0.0.
    """

    format = "TRZ"

    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps'}

    @store_init_arguments
    def __init__(self, trzfilename, n_atoms=None, **kwargs):
        """Creates a TRZ Reader

        Parameters
        ----------
        trzfilename : str
            name of input file
        n_atoms : int
            number of atoms in trajectory, must be taken from topology file!
        convert_units : bool (optional)
            converts units to MDAnalysis defaults

        Raises
        ------
        ValueError
           If `n_atoms` or the number of atoms in the topology file do not
           match the number of atoms in the trajectory.
        """
        wmsg = ("The TRZ reader is deprecated and will be removed in "
                "MDAnalysis version 3.0.0")
        warnings.warn(wmsg, DeprecationWarning)

        super(TRZReader, self).__init__(trzfilename,  **kwargs)

        if n_atoms is None:
            raise ValueError('TRZReader requires the n_atoms keyword')

        self.trzfile = util.anyopen(self.filename, 'rb')
        self._cache = dict()
        self._n_atoms = n_atoms

        self._read_trz_header()
        self.ts = Timestep(self.n_atoms,
                           velocities=True,
                           forces=self.has_force,
                           reader=self,
                           **self._ts_kwargs)

        # structured dtype of a single trajectory frame
        readarg = str(n_atoms) + '<f4'
        frame_contents = [
            ('p1', '<i4'),
            ('nframe', '<i4'),
            ('ntrj', '<i4'),
            ('natoms', '<i4'),
            ('treal', '<f8'),
            ('p2', '<2i4'),
            ('box', '<9f8'),
            ('p3', '<2i4'),
            ('pressure', '<f8'),
            ('ptensor', '<6f8'),
            ('p4', '<3i4'),
            ('etot', '<f8'),
            ('ptot', '<f8'),
            ('ek', '<f8'),
            ('T', '<f8'),
            ('p5', '<6i4'),
            ('rx', readarg),
            ('pad2', '<2i4'),
            ('ry', readarg),
            ('pad3', '<2i4'),
            ('rz', readarg),
            ('pad4', '<2i4'),
            ('vx', readarg),
            ('pad5', '<2i4'),
            ('vy', readarg),
            ('pad6', '<2i4'),
            ('vz', readarg)]
        if not self.has_force:
            frame_contents += [('pad7', '<i4')]
        else:
            frame_contents += [
                ('pad7', '<2i4'),
                ('fx', readarg),
                ('pad8', '<2i4'),
                ('fy', readarg),
                ('pad9', '<2i4'),
                ('fz', readarg),
                ('pad10', '<i4')]
        self._dtype = np.dtype(frame_contents)

        self._read_next_timestep()

    def _read_trz_header(self):
        """Reads the header of the trz trajectory"""
        self._headerdtype = np.dtype([
            ('p1', '<i4'),
            ('title', '80c'),
            ('p2', '<2i4'),
            ('force', '<i4'),
            ('p3', '<i4')])
        data = np.fromfile(self.trzfile, dtype=self._headerdtype, count=1)
        self.title = ''.join(c.decode('utf-8') for c in data['title'][0]).strip()
        if data['force'] == 10:
            self.has_force = False
        elif data['force'] == 20:
            self.has_force = True
        else:
            raise IOError

    def _read_next_timestep(self, ts=None):
        if ts is None:
            ts = self.ts

        try:
            data = np.fromfile(self.trzfile, dtype=self._dtype, count=1)
            if data['natoms'][0] != self.n_atoms:
                raise ValueError("Supplied n_atoms {} is incompatible "
                                 "with provided trajectory file. "
                                 "Maybe `topology` is wrong?".format(
                                                             self.n_atoms))
            ts.frame = data['nframe'][0] - 1  # 0 based for MDA
            ts.data['frame'] = data['ntrj'][0] # moved from attr to data
            ts.time = data['treal'][0]
            ts.dimensions = triclinic_box(*(data['box'].reshape(3, 3)))
            ts.data['pressure'] = data['pressure']
            ts.data['pressure_tensor'] = data['ptensor']
            ts.data['total_energy'] = data['etot']
            ts.data['potential_energy'] = data['ptot']
            ts.data['kinetic_energy'] = data['ek']
            ts.data['temperature'] = data['T']
            ts._x[:] = data['rx']
            ts._y[:] = data['ry']
            ts._z[:] = data['rz']
            ts._velocities[:, 0] = data['vx']
            ts._velocities[:, 1] = data['vy']
            ts._velocities[:, 2] = data['vz']
            if self.has_force:
                ts._forces[:, 0] = data['fx']
                ts._forces[:, 1] = data['fy']
                ts._forces[:, 2] = data['fz']
        except IndexError: # Raises indexerror if data has no data (EOF)
            raise IOError from None
        else:
            # Convert things read into MDAnalysis' native formats (nm -> angstroms)
            if self.convert_units:
                self.convert_pos_from_native(self.ts._pos)
                if self.ts.dimensions is not None:
                    self.convert_pos_from_native(self.ts.dimensions[:3])
                self.convert_velocities_from_native(self.ts._velocities)

            return ts

    @property
    def n_atoms(self):
        """Number of atoms in a frame"""
        return self._n_atoms

    @property
    @cached('n_frames')
    def n_frames(self):
        """Total number of frames in a trajectory"""
        try:
            return self._read_trz_n_frames(self.trzfile)
        except IOError:
            return 0

    def _read_trz_n_frames(self, trzfile):
        """Uses size of file and dtype information to determine how many frames exist

        .. versionchanged:: 0.9.0
           Now is based on filesize rather than reading entire file
        """
        fsize = os.fstat(trzfile.fileno()).st_size  # size of file in bytes

        if not (fsize - self._headerdtype.itemsize) % self._dtype.itemsize == 0:
            raise IOError("Trajectory has incomplete frames")  # check that division is sane

        nframes = int((fsize - self._headerdtype.itemsize) / self._dtype.itemsize)  # returns long int otherwise

        return nframes

    def _get_dt(self):
        """The amount of time between frames in ps

        Assumes that this step is constant (ie. 2 trajectories with different
        steps haven't been stitched together).
        Returns ``AttributeError`` in case of ``StopIteration``
        (which makes :attr:`dt` return 1.0).

        .. versionchanged:: 2.1.0
           Now returns an ``AttributeError`` if dt can't be obtained from the
           time difference between two frames.
        """
        curr_frame = self.ts.frame
        try:
            t0 = self.ts.time
            self.next()
            t1 = self.ts.time
            dt = t1 - t0
        except StopIteration:
            raise AttributeError
        else:
            return dt
        finally:
            self._read_frame(curr_frame)

    @property
    @cached('delta')
    def delta(self):
        """MD integration timestep"""
        return self.dt / self.skip_timestep

    @property
    @cached('skip_timestep')
    def skip_timestep(self):
        """Timesteps between trajectory frames"""
        curr_frame = self.ts.frame
        try:
            t0 = self.ts.data['frame']
            self.next()
            t1 = self.ts.data['frame']
            skip_timestep = t1 - t0
        except StopIteration:
            return 0
        else:
            return skip_timestep
        finally:
            self._read_frame(curr_frame)

    def _read_frame(self, frame):
        """Move to *frame* and fill timestep with data.

        .. versionchanged:: 0.11.0
           Frames now 0-based instead of 1-based
        """
        move = frame - self.ts.frame

        self._seek(move - 1)
        self._read_next_timestep()
        return self.ts

    def _seek(self, nframes):
        """Move *nframes* in the trajectory

        Note that this doens't read the trajectory (ts remains unchanged)

        .. versionadded:: 0.9.0
        """
        framesize = self._dtype.itemsize
        seeksize = framesize * nframes
        maxi_l = seeksize + 1

        if seeksize > maxi_l:
            # Workaround for seek not liking long ints
            # On python 3 this branch will never be used as we defined maxi_l
            # greater than seeksize.
            framesize = long(framesize)
            seeksize = framesize * nframes

            nreps = int(seeksize / maxi_l)  # number of max seeks we'll have to do
            rem = int(seeksize % maxi_l)  # amount leftover to do once max seeks done

            for _ in range(nreps):
                self.trzfile.seek(maxint, 1)
            self.trzfile.seek(rem, 1)
        else:
            seeksize = int(seeksize)

            self.trzfile.seek(seeksize, 1)

    def _reopen(self):
        self.close()
        self.open_trajectory()
        self._read_trz_header()  # Moves to start of first frame

    def open_trajectory(self):
        """Open the trajectory file"""
        if self.trzfile is not None:
            raise IOError(errno.EALREADY, 'TRZ file already opened', self.filename)
        if not os.path.exists(self.filename):
            raise IOError(errno.ENOENT, 'TRZ file not found', self.filename)

        self.trzfile = util.anyopen(self.filename, 'rb')

        #Reset ts
        ts = self.ts
        ts.frame = -1

        return self.trzfile

    def Writer(self, filename, n_atoms=None):
        if n_atoms is None:
            # guess that they want to write the whole timestep unless told otherwise?
            n_atoms = self.ts.n_atoms
        return TRZWriter(filename, n_atoms)

    def close(self):
        """Close trz file if it was open"""
        if self.trzfile is not None:
            self.trzfile.close()
            self.trzfile = None


class TRZWriter(base.WriterBase):
    """Writes a TRZ format trajectory.

    Note
    ----
    Binary TRZ trajectories are *always* written in *little-endian* byte order.


    .. deprecated:: 2.7.0
       The TRZ Reader and Writer are deprecated as of version 2.7.0
       and will be removed in version 3.0.0.
    """

    format = 'TRZ'
    multiframe = True

    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps'}

    def __init__(self, filename, n_atoms, title='TRZ', convert_units=True):
        """Create a TRZWriter

        Parameters
        ----------
        filename : str
            name of output file
        n_atoms : int
            number of atoms in trajectory
        title : str (optional)
            title of the trajectory; the title must be 80 characters or
            shorter, a longer title raises a ValueError exception.
        convert_units : bool (optional)
            units are converted to the MDAnalysis base format; [``True``]
        """
        wmsg = ("The TRZ writer is deprecated and will be removed in "
                "MDAnalysis version 3.0.0")
        warnings.warn(wmsg, DeprecationWarning)

        self.filename = filename
        if n_atoms is None:
            raise ValueError("TRZWriter requires the n_atoms keyword")
        if n_atoms == 0:
            raise ValueError("TRZWriter: no atoms in output trajectory")
        self.n_atoms = n_atoms

        if len(title) > 80:
            raise ValueError("TRZWriter: 'title' must be 80 characters of shorter")

        self.convert_units = convert_units

        self.trzfile = util.anyopen(self.filename, 'wb')

        self._writeheader(title)

        floatsize = str(n_atoms) + '<f4'
        self.frameDtype = np.dtype([
            ('p1a', '<i4'),
            ('nframe', '<i4'),
            ('ntrj', '<i4'),
            ('natoms', '<i4'),
            ('treal', '<f8'),
            ('p1b', '<i4'),
            ('p2a', '<i4'),
            ('box', '<9f8'),
            ('p2b', '<i4'),
            ('p3a', '<i4'),
            ('pressure', '<f8'),
            ('ptensor', '<6f8'),
            ('p3b', '<i4'),
            ('p4a', '<i4'),
            ('six', '<i4'),
            ('etot', '<f8'),
            ('ptot', '<f8'),
            ('ek', '<f8'),
            ('T', '<f8'),
            ('blanks', '<2f8'),
            ('p4b', '<i4'),
            ('p5a', '<i4'),
            ('rx', floatsize),
            ('p5b', '<i4'),
            ('p6a', '<i4'),
            ('ry', floatsize),
            ('p6b', '<i4'),
            ('p7a', '<i4'),
            ('rz', floatsize),
            ('p7b', '<i4'),
            ('p8a', '<i4'),
            ('vx', floatsize),
            ('p8b', '<i4'),
            ('p9a', '<i4'),
            ('vy', floatsize),
            ('p9b', '<i4'),
            ('p10a', '<i4'),
            ('vz', floatsize),
            ('p10b', '<i4')])

    def _writeheader(self, title):
        hdt = np.dtype([
            ('pad1', '<i4'), ('title', '80c'), ('pad2', '<i4'),
            ('pad3', '<i4'), ('nrec', '<i4'), ('pad4', '<i4')])
        out = np.zeros((), dtype=hdt)
        out['pad1'], out['pad2'] = 80, 80
        out['title'] = title + ' ' * (80 - len(title))
        out['pad3'], out['pad4'] = 4, 4
        out['nrec'] = 10
        out.tofile(self.trzfile)

    def _write_next_frame(self, obj):
        """Write information associated with ``obj`` at current frame into trajectory

        Parameters
        ----------
        ag : AtomGroup or Universe


        .. versionchanged:: 1.0.0
           Renamed from `write_next_timestep` to `_write_next_frame`.
        .. versionchanged:: 2.0.0
           Deprecated support for Timestep argument has now been removed.
           Use AtomGroup or Universe as an input instead.
        """
        # Check size of ts is same as initial
        try:  # atomgroup?
            ts = obj.ts
        except AttributeError:  # universe?
            try:
                ts = obj.trajectory.ts
            except AttributeError:
                errmsg = "Input obj is neither an AtomGroup or Universe"
                raise TypeError(errmsg) from None

        if not obj.atoms.n_atoms == self.n_atoms:
            errmsg = "Number of atoms in ts different to initialisation"
            raise ValueError(errmsg)

        # Gather data, faking it when unavailable
        data = {}
        faked_attrs = []
        for att in ['pressure', 'pressure_tensor', 'total_energy',
                    'potential_energy', 'kinetic_energy', 'temperature']:
            try:
                data[att] = ts.data[att]
            except KeyError:
                if att == 'pressure_tensor':
                    data[att] = np.zeros(6, dtype=np.float64)
                else:
                    data[att] = 0.0
                faked_attrs.append(att)
        try:
            data['step'] = ts.data['frame']
        except KeyError:
            data['step'] = ts.frame
            faked_attrs.append('step')
        try:
            data['time'] = ts.time
        except AttributeError:
            data['time'] = float(ts.frame)
            faked_attrs.append('time')
        if faked_attrs:
            warnings.warn("Timestep didn't have the following attributes: '{0}', "
                          "these will be set to 0 in the output trajectory"
                          "".format(", ".join(faked_attrs)))

        # Convert other stuff into our format
        if ts.dimensions is not None:
            unitcell = triclinic_vectors(ts.dimensions).reshape(9)
        else:
            warnings.warn("Timestep didn't have dimensions information, "
                          "box will be written as all zero values")
            unitcell = np.zeros(9, dtype=np.float32)

        try:
            vels = ts.velocities
        except NoDataError:
            vels = np.zeros((self.n_atoms, 3), dtype=np.float32, order='F')
            warnings.warn("Timestep didn't have velocity information, "
                          "this will be set to zero in output trajectory. ")

        out = np.zeros((), dtype=self.frameDtype)
        out['p1a'], out['p1b'] = 20, 20
        out['nframe'] = ts.frame + 1  # TRZ wants 1 based
        out['ntrj'] = data['step']
        out['natoms'] = self.n_atoms
        out['treal'] = data['time']
        out['p2a'], out['p2b'] = 72, 72
        out['box'] = self.convert_pos_to_native(unitcell, inplace=False)
        out['p3a'], out['p3b'] = 56, 56
        out['pressure'] = data['pressure']
        out['ptensor'] = data['pressure_tensor']
        out['p4a'], out['p4b'] = 60, 60
        out['six'] = 6
        out['etot'] = data['total_energy']
        out['ptot'] = data['potential_energy']
        out['ek'] = data['kinetic_energy']
        out['T'] = data['temperature']
        out['blanks'] = 0.0, 0.0
        size = ts.n_atoms * 4  # size of float for vels & coords
        out['p5a'], out['p5b'] = size, size
        out['rx'] = self.convert_pos_to_native(ts._x, inplace=False)
        out['p6a'], out['p6b'] = size, size
        out['ry'] = self.convert_pos_to_native(ts._y, inplace=False)
        out['p7a'], out['p7b'] = size, size
        out['rz'] = self.convert_pos_to_native(ts._z, inplace=False)
        out['p8a'], out['p8b'] = size, size
        out['vx'] = self.convert_velocities_to_native(vels[:, 0], inplace=False)
        out['p9a'], out['p9b'] = size, size
        out['vy'] = self.convert_velocities_to_native(vels[:, 1], inplace=False)
        out['p10a'], out['p10b'] = size, size
        out['vz'] = self.convert_velocities_to_native(vels[:, 2], inplace=False)
        out.tofile(self.trzfile)

    def close(self):
        """Close if it was open"""
        if self.trzfile is None:
            return
        self.trzfile.close()
        self.trzfile = None
