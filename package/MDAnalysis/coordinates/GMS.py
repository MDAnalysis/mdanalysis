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

"""GMS trajectory reader --- :mod:`MDAnalysis.coordinates.GMS`
==============================================================

Resources: the GMS output format is a common output format for different
GAMESS distributions: US-GAMESS, Firefly (PC-GAMESS) and GAMESS-UK.

Current version was approbated with US-GAMESS & Firefly only.

There appears to be no rigid format definition so it is likely users
will need to tweak this Class.
"""

import MDAnalysis.core.util as util
import os, errno, re
import numpy

import base
from base import Timestep


class GMSReader(base.Reader):
    """Reads from an GAMESS output file

    :Data:
        ts
          Timestep object containing coordinates of current frame

    :Methods:
        ``len(out)``
          return number of frames in out
        ``for ts in out:``
          iterate through trajectory

    .. Note: this can read both compressed (foo.out) and compressed
          (foo.out.bz2 or foo.out.gz) files; uncompression is handled
          on the fly
    """

    # this will be overidden when an instance is created and the file extension checked
    format = "GMS"

    # these are assumed!
    units = {'time': 'ps', 'length': 'Angstrom'}

    def __init__(self, outfilename, **kwargs):
        self.filename = outfilename

        # the filename has been parsed to be either b(g)zipped or not
        self.outfile = util.anyopen(self.filename, 'r')

        # note that, like for xtc and trr files, _numatoms and _numframes are used quasi-private variables
        # to prevent the properties being recalculated
        # this is because there is no indexing so the way it measures the number of frames is to read the whole file!
        self._numatoms = None
        self._numframes = None
        self._runtyp = None

        self.ts = Timestep(0) # need for properties initial calculations
        self.fixed = 0
        self.skip = 1
        self.periodic = False
        self.delta = kwargs.pop("delta", 1.0)   # there is no time, so let delta be 1
        # update runtyp property
        self.runtyp 
        if not self.runtyp in ['optimize', 'surface']:
            raise AttributeError('Wrong RUNTYP= '+self.runtyp)
        self.skip_timestep = 1

        self.ts = Timestep(self.numatoms)
        # update numframes property
        self.numframes

        # Read in the first timestep

        self._read_next_timestep()

    @property
    def runtyp(self):
        """RUNTYP property of the GAMESS run"""
        if not self._runtyp is None:   # return cached value
            return self._runtyp
        try:
            self._runtyp = self._determine_runtyp()
        except IOError:
            return 0
        else:
            return self._runtyp


    def _determine_runtyp(self):
        self._reopen()
        counter = 0
        for line in self.outfile:
            m = re.match(r'^.*RUNTYP=([A-Z]+)\s+.*', line)
            if (m != None):
                self.close()
                return m.group(1).lower()

        self.close()
        raise EOFError



    @property
    def numatoms(self):
        """number of atoms in a frame"""
        if not self._numatoms is None:   # return cached value
            return self._numatoms
        try:
            self._numatoms = self._read_out_natoms()
        except IOError:
            return 0
        else:
            return self._numatoms

    def _read_out_natoms(self):
        self._reopen()
        # this assumes that this is only called once at startup and that the filestream is already open
        for line in self.outfile:           
            m = re.match(r'\s*TOTAL NUMBER OF ATOMS\s*=\s*([0-9]+)\s*',line)
            if m == None:
                continue
            self.close()
            return int(m.group(1))

        self.close()
        raise EOFError


    @property
    def numframes(self):
        if not self._numframes is None:   # return cached value
            return self._numframes
        try:
            self._numframes = self._read_out_numframes(self.filename)
        except IOError:
            return 0
        else:
            return self._numframes

    def _read_out_numframes(self, filename):
        self._reopen()
        counter = 0
        if self.runtyp == 'optimize':
            for line in self.outfile:           
                if re.match(r'^.NSERCH=.*', line):
                    counter += 1
        elif self.runtyp == 'surface':
            for line in self.outfile:           
                if re.match(r'^.COORD 1=.*', line):
                    counter += 1

        self.close()
        return int(counter)


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
        # check that the outfile object exists; if not reopen the trajectory
        if self.outfile is None:
            self.open_trajectory()
        x = []
        y = []
        z = []
        
        flag = 0
        counter = 0

        for line in self.outfile:
            if self.runtyp == 'optimize': 
                if (flag == 0) and (re.match(r'^.NSERCH=.*', line) != None):
                    flag = 1
                    continue
                if (flag == 1) and (re.match(r'^ COORDINATES OF ALL ATOMS ARE ',\
                    line) != None):
                    flag = 2
                    continue
                if (flag == 2) and (re.match(r'^\s*[-]+\s*', line) != None):
                    flag = 3
                    continue
                if flag == 3 and counter < self.numatoms:
                    words = line.split()
                    x.append(float(words[2]))
                    y.append(float(words[3]))
                    z.append(float(words[4]))
                    counter += 1

            elif self.runtyp == 'surface':
                if (flag == 0) and (re.match(\
                        r'^.COORD 1=\s*([-]?[0-9]+\.[0-9]+).*', line) != None):
                    flag = 1
                    continue
                if (flag == 1) and (re.match(\
                        r'^\s*HAS ENERGY VALUE\s*([-]?[0-9]+\.[0-9]+)\s*', line) != None):
                    flag = 3
                    continue
                if flag == 3 and counter < self.numatoms:
                    words = line.split()
                    x.append(float(words[1]))
                    y.append(float(words[2]))
                    z.append(float(words[3]))
                    counter += 1

            # stop when the cursor has reached the end of that block
            if counter == self._numatoms:
                ts._unitcell = numpy.zeros((6), numpy.float)
                ts._x[:] = x # more efficient to do it this way to avoid re-creating the numpy arrays
                ts._y[:] = y
                ts._z[:] = z
                #print ts.frame
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
        if not self.outfile is None:
            raise IOError(errno.EALREADY, 'GMS file already opened', self.filename)
        if not os.path.exists(self.filename):
            # must check; otherwise might segmentation fault
            raise IOError(errno.ENOENT, 'GMS file not found', self.filename)

        self.outfile = util.anyopen(self.filename, 'r')

        # reset ts
        ts = self.ts
        ts.status = 1
        ts.frame = 0
        ts.step = 0
        ts.time = 0
        return self.outfile

    def close(self):
        """Close out trajectory file if it was open."""
        if self.outfile is None:
            return
        self.outfile.close()
        self.outfile = None

    def __del__(self):
        self.close()

