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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""
TXYZ file format --- :mod:`MDAnalysis.coordinates.TXYZ`
=======================================================

Coordinate reader for Tinker_ xyz files .txyz and trajectory .arc files.
Differences between Tinker format_ and normal xyz files:

- the 1st header line contains both the number of atoms and a comment
- the 2nd header line (optional) has periodic boundary conditions information
- column 1 contains atom numbers (starting from 1)
- column 6 contains atoms types
- the following columns indicate connectivity (atoms to which that particular atom is
  bonded, according to numbering in column 1)

.. _format: http://chembytes.wikidot.com/tnk-tut00#toc2
.. _Tinker: https://dasher.wustl.edu/tinker/


Classes
-------

.. autoclass:: TXYZWriter
   :members:
   :inherited-members:

.. autoclass:: TXYZReader
   :members:
   :inherited-members:

"""
from __future__ import absolute_import, division
from six.moves import range

import numpy as np
import os
import errno

from . import base
from ..core import flags
from ..lib import util
from ..lib.util import openany, cached
from ..exceptions import NoDataError
from ..version import __version__


class TXYZWriter(base.WriterBase):
    """Writes a TXYZ file (single frame)
       or an ARC file (multiple frames)
    """

    format = 'TXYZ', 'ARC'
    multiframe = True
    # these are assumed!
    units = {'time': 'ps', 'length': 'Angstrom'}

    def __init__(self, filename, n_atoms=None, atoms=None, convert_units=None,
                 remark=None, **kwargs):
        """Initialize the TXYZ trajectory writer

        Parameters
        ----------
        filename: str
            filename of trajectory file. If it ends with "gz" then the file
            will be gzip-compressed; if it ends with "bz2" it will be bzip2
            compressed.
        n_atoms: int (optional)
            Number of atoms in trajectory. By default assume that this is None
            and that this file is used to store several different models
            instead of a single trajectory. If a number is provided each
            written TimeStep has to contain the same number of atoms.
        atoms: str | list (optional)
            Provide atom names: This can be a list of names or an
            :class:`AtomGroup`.  If none is provided, atoms will
            be called 'X' in the output. These atom names will be
            used when a trajectory is written from raw
            :class:`Timestep` objects which do not contain atom
            information. If you write a :class:`AtomGroup` with
            :meth:`XYZWriter.write` then atom information is taken
            at each step and *atoms* is ignored.
        remark: str (optional)
            text following n_atoms ("molecule name"). By default writes MDAnalysis
            version
        """
        self.filename = filename
        self.n_atoms = n_atoms
        if convert_units is not None:
            self.convert_units = convert_units
        else:
            self.convert_units = flags['convert_lengths']
        self.atomnames = self._get_atomnames(atoms)
        default_remark = "Written by {0} (release {1})".format(
            self.__class__.__name__, __version__)
        self.remark = default_remark if remark is None else remark
        # can also be gz, bz2
        self._txyz = util.anyopen(self.filename, 'wt')

    def _get_atomnames(self, atoms):
        """Return a list of atom names"""
        # Default case
        if atoms is None:
            return itertools.cycle(('X',))
        # Single atom name provided
        elif isinstance(atoms, six.string_types):
            return itertools.cycle((atoms,))
        # List of atom names providded
        elif isinstance(atoms, list):
            return atoms
        # AtomGroup or Universe, grab the names else default
        # (AtomGroup.atoms just returns AtomGroup)
        try:
            return atoms.atoms.names
        except (AttributeError, NoDataError):
            return itertools.cycle(('X',))

    def close(self):
        """Close the trajectory file and finalize the writing"""
        if self._txyz is not None:
            self._txyz.write("\n")
            self._txyz.close()
        self._txyz = None

    def write(self, obj):
        """Write object `obj` at current trajectory frame to file.

        Atom names in the output are taken from the `obj` or default
        to the value of the `atoms` keyword supplied to the
        :class:`TXYZWriter` constructor.

        Parameters
        ----------
        obj : Universe or AtomGroup
            The :class:`~MDAnalysis.core.groups.AtomGroup` or
            :class:`~MDAnalysis.core.universe.Universe` to write.
        """
        # prepare the Timestep and extract atom names if possible
        # (The way it is written it should be possible to write
        # trajectories with frames that differ in atom numbers
        # but this is not tested.)
        try:
            atoms = obj.atoms
        except AttributeError:
            if isinstance(obj, base.Timestep):
                ts = obj
            else:
                raise TypeError("No Timestep found in obj argument")
        else:
            if hasattr(obj, 'universe'):
                # For AtomGroup and children (Residue, ResidueGroup, Segment)
                ts_full = obj.universe.trajectory.ts
                if ts_full.n_atoms == atoms.n_atoms:
                    ts = ts_full
                else:
                    # Only populate a time step with the selected atoms.
                    ts = ts_full.copy_slice(atoms.indices)
            elif hasattr(obj, 'trajectory'):
                # For Universe only --- get everything
                ts = obj.trajectory.ts
            # update atom names
            self.atomnames = self._get_atomnames(atoms)

        self.write_next_timestep(ts)

    def write_next_timestep(self, ts=None):
        """Write coordinate information in *ts* to the trajectory"""
        if ts is None:
            if not hasattr(self, 'ts'):
                raise NoDataError('TXYZWriter: no coordinate data to write to '
                                  'trajectory file')
            else:
                ts = self.ts

        if self.n_atoms is not None:
            if self.n_atoms != ts.n_atoms:
                raise ValueError('n_atoms keyword was specified indicating '
                                 'that this should be a trajectory of the '
                                 'same model. But the provided TimeStep has a '
                                 'different number ({}) then expected ({})'
                                 ''.format(ts.n_atoms, self.n_atoms))
        else:
            if (not isinstance(self.atomnames, itertools.cycle) and
                len(self.atomnames) != ts.n_atoms):
                logger.info('Trying to write a TimeStep with unkown atoms. '
                            'Expected {}, got {}. Try using "write" if you are '
                            'using "write_next_timestep" directly'.format(
                                len(self.atomnames), ts.n_atoms))
                self.atomnames = np.array([self.atomnames[0]] * ts.n_atoms)

        if self.convert_units:
            coordinates = self.convert_pos_to_native(
                ts.positions, inplace=False)
        else:
            coordinates = ts.positions

        self._xyz.write("{0:d} frame {1}\n".format(ts.n_atoms, ts.frame))
        self._xyz.write("{0:10.5f} {1:10.5f} {2:10.5f} {3:10.5f}
                         {4:10.5f} {5:10.5f} {6:10.5f}\n".format(ts.dimensions))
        for atom, (x, y, z), typ  in zip(self.atomnames, coordinates, self.atomtypes):
            self._xyz.write("{0!s:>8}  {1:10.5f} {2:10.5f} {3:10.5f} {4}"
                            "".format(atom, x, y, z, typ))

class TXYZReader(base.ReaderBase):
    """Reads from a TXYZ file"""


    format = ['TXYZ', 'ARC']
    # these are assumed!
    units = {'time': 'ps', 'length': 'Angstrom'}
    _Timestep = base.Timestep

    def __init__(self, filename, **kwargs):
        super(TXYZReader, self).__init__(filename, **kwargs)

        # the filename has been parsed to be either be foo.xyz or foo.xyz.bz2 by
        # coordinates::core.py so the last file extension will tell us if it is
        # bzipped or not
        root, ext = os.path.splitext(self.filename)
        self.txyzfile = util.anyopen(self.filename)
        self._cache = dict()
        with util.openany(self.filename) as inp:
           inp.readline()
           line=inp.readline()
           try:
               float(line.split()[1])
           except ValueError:
               self.periodic=False
           else:
               self.periodic=True
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        # Haven't quite figured out where to start with all the self._reopen()
        # etc.
        # (Also cannot just use seek() or reset() because that would break
        # with urllib2.urlopen() streams)
        self._read_next_timestep()

    @property
    @cached('n_atoms')
    def n_atoms(self):
        """number of atoms in a frame"""
        with util.anyopen(self.filename) as f:
            n = f.readline().split()[0]
        # need to check type of n
        return int(n)

    @property
    @cached('n_frames')
    def n_frames(self):
        try:
            return self._read_txyz_n_frames()
        except IOError:
            return 0

    def _read_txyz_n_frames(self):
        # the number of lines in the XYZ file will be 1 greater than the
        # number of atoms
        linesPerFrame = self.n_atoms + 1
        if self.periodic:
            linesPerFrame +=1
        counter = 0
        offsets = []

        with util.anyopen(self.filename) as f:
            line = True
            while line:
                if not counter % linesPerFrame:
                    offsets.append(f.tell())
                line = f.readline()
                counter += 1

        # need to check this is an integer!
        n_frames = int(counter / linesPerFrame)
        self._offsets = offsets
        return n_frames

    def _read_frame(self, frame):
        self.txyzfile.seek(self._offsets[frame])
        self.ts.frame = frame - 1  # gets +1'd in next
        return self._read_next_timestep()

    def _read_next_timestep(self, ts=None):
        # check that the timestep object exists
        if ts is None:
            ts = self.ts

        f = self.txyzfile

        try:
            # we assume that there is only one header line per frame
            f.readline()
            if self.periodic:
                ts.dimensions=f.readline().split()
            # convert all entries at the end once for optimal speed
            tmp_buf = []
            for i in range(self.n_atoms):
                tmp_buf.append(f.readline().split()[2:5])
            ts.positions = tmp_buf
            ts.frame += 1
            return ts
        except (ValueError, IndexError) as err:
            raise EOFError(err)

    def _reopen(self):
        self.close()
        self.open_trajectory()

    def open_trajectory(self):
        if self.txyzfile is not None:
            raise IOError(
                errno.EALREADY, 'TXYZ file already opened', self.filename)

        self.txyzfile = util.anyopen(self.filename)

        # reset ts
        ts = self.ts
        ts.frame = -1

        return self.txyzfile

    def close(self):
        """Close arc trajectory file if it was open."""
        if self.txyzfile is None:
            return
        self.txyzfile.close()
        self.txyzfile = None
