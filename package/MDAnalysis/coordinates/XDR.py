# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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
XDR based trajectory files --- :mod:`MDAnalysis.coordinates.XDR`
================================================================

This module contains helper function and classes to read the XTC and TRR file
formats.

See Also
--------
MDAnalysis.coordinates.XTC: Read and write GROMACS XTC trajectory files.
MDAnalysis.coordinates.TRR: Read and write GROMACS TRR trajectory files.
MDAnalysis.lib.formats.libmdaxdr: Low level xdr format reader
"""

import errno
import numpy as np
from os.path import getctime, getsize, isfile, split, join
import warnings
from filelock import FileLock

from . import base
from ..lib.mdamath import triclinic_box
from ..lib.util import store_init_arguments


def offsets_filename(filename, ending="npz"):
    """Return offset or its lock filename for XDR files.
    For this the filename is appended
    with `_offsets.{ending}`.

    Parameters
    ----------
    filename : str
        filename of trajectory
    ending : str (optional)
        fileending of offsets file

    Returns
    -------
    offset_filename : str

    """
    head, tail = split(filename)
    return join(head, f".{tail}_offsets.{ending}")


def read_numpy_offsets(filename):
    """read offsets into dictionary.

    This assume offsets have been saved using numpy

    Parameters
    ----------
    filename : str
        filename of offsets

    Returns
    -------
    offsets : dict
        dictionary of offsets information

    """
    try:
        return {k: v for k, v in np.load(filename).items()}

    #  `ValueError` is encountered when the offset file is corrupted.
    except (ValueError, IOError):
        warnings.warn(f"Failed to load offsets file {filename}\n")
        return False


class XDRBaseReader(base.ReaderBase):
    """Base class for libmdaxdr file formats xtc and trr

    This class handles integration of XDR based formats into MDAnalysis. The
    XTC and TRR classes only implement `_write_next_frame` and
    `_frame_to_ts`.

    .. _offsets-label:

    Notes
    -----
    XDR based readers store persistent offsets on disk. The offsets are used to
    enable access to random frames efficiently. These offsets will be generated
    automatically the  first time the  trajectory is opened.  Generally offsets
    are stored  in hidden  `*_offsets.npz` files.  Afterwards opening  the same
    file again is fast. It sometimes can happen that the stored offsets get out
    off sync with the trajectory they refer to. For this the offsets also store
    the number of atoms, size of the file and last modification time. If any of
    them change  the offsets are recalculated.  Writing of the offset  file can
    fail when the  directory where the trajectory file resides  is not writable
    or if the  disk is full. In this  case a warning message will  be shown but
    the offsets will nevertheless be used during the lifetime of the trajectory
    Reader. However, the  next time the trajectory is opened,  the offsets will
    have to be rebuilt again.

    .. versionchanged:: 1.0.0
       XDR offsets read from trajectory if offsets file read-in fails
    .. versionchanged:: 2.0.0
       Add a InterProcessLock when generating offsets
    .. versionchanged:: 2.4.0
       Use a direct read into ts attributes
    .. versionchanged:: 2.9.0
       Changed fasteners.InterProcessLock() to filelock.FileLock
    """

    @store_init_arguments
    def __init__(
        self,
        filename,
        convert_units=True,
        sub=None,
        refresh_offsets=False,
        dt=None,
        **kwargs,
    ):
        """
        Parameters
        ----------
        filename : str
            trajectory filename
        convert_units : bool (optional)
            convert units to MDAnalysis units
        sub : array_like (optional)
            `sub` is an array of indices to pick out the corresponding
            coordinates and load only them; this requires that the topology
            itself is that of the sub system.
        refresh_offsets : bool (optional)
            force refresh of offsets
        dt : float (optional)
            timestep in MDAnalysis units to load trajectory with;
            if `dt` is ``None``, the time is taken from the xdr file;
            else, the time is set to `dt` * frame
        **kwargs : dict
            General reader arguments.

        """
        super(XDRBaseReader, self).__init__(
            filename, convert_units=convert_units, **kwargs
        )
        self._xdr = self._file(self.filename)

        self._sub = sub
        if self._sub is not None:
            self.n_atoms = len(self._sub)
        else:
            self.n_atoms = self._xdr.n_atoms

        if not refresh_offsets:
            self._load_offsets()
        else:
            self._read_offsets(store=True)
        frame = self._xdr.read()
        try:
            xdr_frame = self._xdr.read()
            if dt is None:
                dt = xdr_frame.time - frame.time
            else:
                self._ts_kwargs["dt"] = dt
            self._xdr.seek(1)
        except StopIteration:
            dt = 0

        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self._frame = 0
        self._frame_to_ts(frame, self.ts)
        # these should only be initialized once
        self.ts.dt = dt
        self.ts.dimensions = triclinic_box(*frame.box)
        if self.convert_units:
            if self.ts.dimensions is not None:
                self.convert_pos_from_native(self.ts.dimensions[:3])

    @classmethod
    def parse_n_atoms(cls, filename, **kwargs):
        with cls._file(filename) as f:
            n_atoms = f.n_atoms
        return n_atoms

    def close(self):
        """close reader"""
        self._xdr.close()

    def _load_offsets(self):
        """load frame offsets from file, reread them from the trajectory if that
        fails. To prevent the competition of generating the same offset file
        from multiple processes, an `InterProcessLock` is used."""
        fname = offsets_filename(self.filename)
        lock_name = offsets_filename(self.filename, ending="lock")

        #  check if the location of the lock is writable.
        try:
            with FileLock(lock_name) as filelock:
                pass
        except OSError as e:
            if isinstance(e, PermissionError) or e.errno == errno.EROFS:
                warnings.warn(
                    f"Cannot write lock/offset file in same location as "
                    f"{self.filename}. Using slow offset calculation."
                )
                self._read_offsets(store=False)
                return
            else:
                raise

        with FileLock(lock_name) as filelock:
            if not isfile(fname):
                self._read_offsets(store=True)
                return

            # if offsets file read correctly, data will be a dictionary of offsets
            # if not, data will be False
            # if False, offsets should be read from the trajectory
            # this warning can be avoided by loading Universe like:
            # u = mda.Universe(data.TPR, data.XTC, refresh_offsets=True)
            # refer to Issue #1893
            data = read_numpy_offsets(fname)
            if not data:
                warnings.warn(
                    f"Reading offsets from {fname} failed, "
                    "reading offsets from trajectory instead.\n"
                    "Consider setting 'refresh_offsets=True' "
                    "when loading your Universe."
                )
                self._read_offsets(store=True)
                return

            ctime_ok = size_ok = n_atoms_ok = False

            try:
                ctime_ok = getctime(self.filename) == data["ctime"]
                size_ok = getsize(self.filename) == data["size"]
                n_atoms_ok = self._xdr.n_atoms == data["n_atoms"]
            except KeyError:
                # we tripped over some old offset formated file
                pass

            if not (ctime_ok and size_ok and n_atoms_ok):
                warnings.warn(
                    "Reload offsets from trajectory\n "
                    "ctime or size or n_atoms did not match"
                )
                self._read_offsets(store=True)
            else:
                self._xdr.set_offsets(data["offsets"])

    def _read_offsets(self, store=False):
        """read frame offsets from trajectory"""
        fname = offsets_filename(self.filename)
        offsets = self._xdr.offsets
        if store:
            ctime = getctime(self.filename)
            size = getsize(self.filename)
            try:
                np.savez(
                    fname,
                    offsets=offsets,
                    size=size,
                    ctime=ctime,
                    n_atoms=self._xdr.n_atoms,
                )
            except Exception as e:
                warnings.warn(f"Couldn't save offsets because: {e}")

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        return len(self._xdr)

    def _reopen(self):
        """reopen trajectory"""
        self.ts.frame = 0
        self._frame = -1
        offsets = self._xdr.offsets.copy()
        self._xdr.close()
        self._xdr.open(self.filename.encode("utf-8"), "r")
        # only restore in case we actually had offsets
        if len(offsets) != 0:
            self._xdr.set_offsets(offsets)

    def _read_frame(self, i):
        """read frame i"""
        self._frame = i - 1
        try:
            self._xdr.seek(i)
            timestep = self._read_next_timestep()
        except IOError:
            warnings.warn("seek failed, recalculating offsets and retrying")
            offsets = self._xdr.calc_offsets()
            self._xdr.set_offsets(offsets)
            self._read_offsets(store=True)
            self._xdr.seek(i)
            timestep = self._read_next_timestep()
        return timestep

    def Writer(self, filename, n_atoms=None, **kwargs):
        """Return writer for trajectory format"""
        if n_atoms is None:
            n_atoms = self.n_atoms
        return self._writer(filename, n_atoms=n_atoms, **kwargs)


class XDRBaseWriter(base.WriterBase):
    """Base class for libmdaxdr file formats xtc and trr"""

    def __init__(
        self, filename, n_atoms, convert_units=True, dt=None, **kwargs
    ):
        """
        Parameters
        ----------
        filename : str
            filename of trajectory
        n_atoms : int
            number of atoms to be written
        convert_units : bool (optional)
            convert from MDAnalysis units to format specific units
        dt : float (optional)
            timestep in MDAnalysis units to write trajectory with;
            if `dt` is ``None``, time for a frame is set from the timestep;
            else, the time for a frame is `dt` * frame
        **kwargs : dict
            General writer arguments
        """
        self.filename = filename
        self._convert_units = convert_units
        self.n_atoms = n_atoms
        self._xdr = self._file(self.filename, "w")
        self._dt = dt

    def close(self):
        """close trajectory"""
        self._xdr.close()

    def __del__(self):
        self.close()
