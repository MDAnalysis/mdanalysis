# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""XYZ trajectory reader --- :mod:`MDAnalysis.coordinates.XYZ`
==============================================================

Resources: the XYZ format was taken from
http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html and
is therefore compatible with VMD (you need a PDB or PSF file to define
the topology, just as here)

* comments are not allowed in the XYZ file
* the atom name (first column) is ignored
* the coordinates are assumed to be space-delimited rather than fixed width (this may cause issues - see below)
* all fields to the right of the z-coordinate are ignored
* it is assumed that the coordinates are in Angstroms
* the unitcell information is all zeros since this is not recorded in the XYZ format
* the length of a timestep can be set by passing the `delta` argument, it's assumed to be in ps (default: 1 ps)

There appears to be no rigid format definition so it is likely users
will need to tweak this Class.
"""

import os, errno
import numpy
import bz2
import gzip

import base
from base import Timestep


class XYZReader(base.Reader):
    """Reads from an XYZ file

    :Data:
        ts
          Timestep object containing coordinates of current frame

    :Methods:
        ``len(xyz)``
          return number of frames in xyz
        ``for ts in xyz:``
          iterate through trajectory

    .. Note: this can read both compressed (foo.xyz) and compressed
          (foo.xyz.bz2 or foo.xyz.gz) files; uncompression is handled
          on the fly

    File format: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html

    Validation: the geometric centre of 1284 atoms was calculated over
    500 frames using both MDAnalysis and a VMD Tcl script. There was
    exact agreement (measured to 3DP). bzipped and gzipped versions of
    the XYZ file were also tested

    """

    # this will be overidden when an instance is created and the file extension checked
    format = "XYZ"

    # these are assumed!
    units = {'time': 'ps', 'length': 'Angstrom'}

    def __init__(self, xyzfilename, **kwargs):
        self.filename = xyzfilename

        # the filename has been parsed to be either be foo.xyz or foo.xyz.bz2 by coordinates::core.py
        # so the last file extension will tell us if it is bzipped or not
        root, ext = os.path.splitext(self.filename)
        if ext[1:] == "bz2":
            self.compression = "bz2"
            self.xyzfile = bz2.BZ2File(self.filename, 'rb')
        elif ext[1:] == "gz":
            self.compression = "gz"
            self.xyzfile = gzip.open(self.filename, 'rb')
        elif ext[1:] == "xyz":
            self.compression = None
            self.xyzfile = open(self.filename, 'r')

        # note that, like for xtc and trr files, __numatoms and __numframes are used quasi-private variables
        # to prevent the properties being recalculated
        # this is because there is no indexing so the way it measures the number of frames is to read the whole file!
        self.__numatoms = None
        self.__numframes = None

        self.fixed = 0
        self.skip = 1
        self.periodic = False
        self.delta = kwargs.pop("delta", 1.0)   # can set delta manually, default is 1ps (taken from TRJReader)
        self.skip_timestep = 1

        self.ts = Timestep(self.numatoms)

        # Read in the first timestep
        self._read_next_timestep()

    @property
    def numatoms(self):
        """number of atoms in a frame"""
        if not self.__numatoms is None:   # return cached value
            return self.__numatoms
        try:
            self.__numatoms = self._read_xyz_natoms(self.filename)
        except IOError:
            return 0
        else:
            return self.__numatoms

    def _read_xyz_natoms(self,filename):

        # this assumes that this is only called once at startup and that the filestream is already open

        # read the first line
        n = self.xyzfile.readline()

        self.close()

        # need to check type of n
        return int(n)

    @property
    def numframes(self):
        if not self.__numframes is None:   # return cached value
            return self.__numframes
        try:
            self.__numframes = self._read_xyz_numframes(self.filename)
        except IOError:
            return 0
        else:
            return self.__numframes

    def _read_xyz_numframes(self, filename):
        self._reopen()

        # the number of lines in the XYZ file will be 2 greater than the number of atoms
        linesPerFrame = self.numatoms+2

        counter = 0
        # step through the file (assuming xyzfile has an iterator)
        for i in self.xyzfile:
            counter = counter + 1

        self.close()

        # need to check this is an integer!
        numframes = int(counter/linesPerFrame)

        return numframes


    def __iter__(self):
        self.ts.frame = 0  # start at 0 so that the first frame becomes 1
        self._reopen()
        while True:
            try:
                yield self._read_next_timestep()
            except EOFError:
                self.close()
                raise StopIteration

    def _read_next_timestep(self, ts=None):
        # check that the timestep object exists
        if ts is None:
            ts = self.ts
        # check that the xyzfile object exists; if not reopen the trajectory
        if self.xyzfile is None:
            self.open_trajectory()
        x = []
        y = []
        z = []

        # we assume that there are only two header lines per frame
        counter = -2

        for line in self.xyzfile:

            counter += 1

            if counter > 0:

                # assume the XYZ file is space delimited rather than being fixed format
                # (this could lead to problems where there is no gap e.g 9.768-23.4567)
                words = line.split()

                x.append(float(words[1]))
                y.append(float(words[2]))
                z.append(float(words[3]))

            # stop when the cursor has reached the end of that block
            if counter == self.numatoms:
                ts._unitcell = numpy.zeros((6), numpy.float32)
                ts._x[:] = x # more efficient to do it this way to avoid re-creating the numpy arrays
                ts._y[:] = y
                ts._z[:] = z
                ts.frame += 1
                return ts

        raise EOFError

    def rewind(self):
        """reposition on first frame"""
        self._reopen()
        # the next method is inherited from the Reader Class and calls _read_next_timestep
        self.next()

    def _reopen(self):
        self.close()
        self.open_trajectory()

    def open_trajectory(self):
        if not self.xyzfile is None:
            raise IOError(errno.EALREADY, 'XYZ file already opened', self.filename)
        if not os.path.exists(self.filename):
            # must check; otherwise might segmentation fault
            raise IOError(errno.ENOENT, 'XYZ file not found', self.filename)

        if self.compression == "bz2":
            self.xyzfile = bz2.BZ2File(self.filename, 'rb')
        elif self.compression == "gz":
            self.xyzfile = gzip.open(self.filename, 'rb')
        elif self.compression == None:
            self.xyzfile = open(self.filename, 'r')

        # reset ts
        ts = self.ts
        ts.status = 1
        ts.frame = 0
        ts.step = 0
        ts.time = 0
        return self.xyzfile

    def close(self):
        """Close xyz trajectory file if it was open."""
        if self.xyzfile is None:
            return
        self.xyzfile.close()
        self.xyzfile = None

    def __del__(self):
        if not self.xyzfile is None:
            self.close()

