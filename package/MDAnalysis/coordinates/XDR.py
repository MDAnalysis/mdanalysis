# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.


# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import six

import errno
import numpy as np
from os.path import getctime, getsize, isfile, split, join
import warnings

from . import base
from ..lib.mdamath import triclinic_box


def offsets_filename(filename, ending='npz'):
    """Return offset filename

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
    return join(head, '.{tail}_offsets.{ending}'.format(tail=tail,
                                                        ending=ending))


def read_numpy_offsets(filename):
    """read offsets into a dictionary

    Parameters
    ----------
    filename : str
        filename of offsets

    Returns
    -------
    offsets : dict
        dictionary of offsets information
    """
    return {k: v for k, v in six.iteritems(np.load(filename))}


class XDRBaseReader(base.Reader):
    """Base class for libmdaxdr file formats xtc and trr"""
    def __init__(self, filename, convert_units=True, sub=None,
                 refresh_offsets=False, **kwargs):
        super(XDRBaseReader, self).__init__(filename,
                                            convert_units=convert_units,
                                            **kwargs)
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
            dt = xdr_frame.time - frame.time
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
            self.convert_pos_from_native(self.ts.dimensions[:3])

    def close(self):
        """close reader"""
        self._xdr.close()

    def _load_offsets(self):
        """load frame offsets from file, reread them from the trajectory if that
        fails"""
        fname = offsets_filename(self.filename)

        if not isfile(fname):
            self._read_offsets(store=True)
            return

        data = read_numpy_offsets(fname)
        ctime_ok = size_ok = n_atoms_ok = False

        try:
            ctime_ok = getctime(self.filename) == data['ctime']
            size_ok = getsize(self.filename) == data['size']
            n_atoms_ok = self._xdr.n_atoms == data['n_atoms']
        except KeyError:
            # we tripped over some old offset formated file
            pass

        if not (ctime_ok and size_ok and n_atoms_ok):
            warnings.warn("Reload offsets from trajectory\n "
                          "ctime or size or n_atoms did not match")
            self._read_offsets(store=True)
        else:
            self._xdr.set_offsets(data['offsets'])

    def _read_offsets(self, store=False):
        """read frame offsets from trajectory"""
        offsets = self._xdr.offsets
        if store:
            ctime = getctime(self.filename)
            size = getsize(self.filename)
            try:
                np.savez(offsets_filename(self.filename),
                         offsets=offsets, size=size, ctime=ctime,
                         n_atoms=self._xdr.n_atoms)
            except Exception as e:
                warnings.warn("Couldn't save offsets because: {}".format(e))

    def rewind(self):
        """Read the first frame again"""
        self._read_frame(0)

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        return len(self._xdr)

    def _reopen(self):
        """reopen trajectory"""
        self.ts.frame = 0
        self._frame = -1
        self._xdr.close()
        self._xdr.open(self.filename.encode('utf-8'), 'r')

    def _read_frame(self, i):
        """read frame i"""
        self._frame = i - 1
        try:
            self._xdr.seek(i)
            timestep = self._read_next_timestep()
        except IOError:
            warnings.warn('seek failed, recalculating offsets and retrying')
            offsets = self._xdr.calc_offsets()
            self._xdr.set_offsets(offsets)
            self._read_offsets(store=True)
            self._xdr.seek(i)
            timestep = self._read_next_timestep()
        return timestep

    def _read_next_timestep(self, ts=None):
        """copy next frame into timestep"""
        if self._frame == self.n_frames - 1:
            raise IOError(errno.EIO, 'trying to go over trajectory limit')
        if ts is None:
            ts = self.ts
        frame = self._xdr.read()
        self._frame += 1
        self._frame_to_ts(frame, ts)
        return ts

    def Writer(self, filename, n_atoms=None, **kwargs):
        """Return writer for trajectory format"""
        if n_atoms is None:
            n_atoms = self.n_atoms
        return self._writer(filename, n_atoms=n_atoms, **kwargs)


class XDRBaseWriter(base.Writer):
    """Base class for libmdaxdr file formats xtc and trr"""

    def __init__(self, filename, n_atoms, convert_units=True, **kwargs):
        self.filename = filename
        self._convert_units = convert_units
        self.n_atoms = n_atoms
        self._xdr = self._file(self.filename, 'w')

    def close(self):
        """close trajectory"""
        self._xdr.close()

    def __del__(self):
        self.close()
