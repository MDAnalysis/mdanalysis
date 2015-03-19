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


"""XYZ trajectory reader --- :mod:`MDAnalysis.coordinates.XYZ`
==============================================================

The :ref:`XYZ format <xyz-format>` is a loosely defined, simple
coordinate trajectory format. The implemented format definition was
taken from the `VMD xyzplugin`_ and is therefore compatible with VMD.

Note the following:
* comments are not allowed in the XYZ file (we neither read nor write
  them to remain compatible with VMD)
* the atom name (first column) is ignored during reading
* the coordinates are assumed to be space-delimited rather than fixed
  width (this may cause issues - see below)
* all fields to the right of the z-coordinate are ignored
* the unitcell information is all zeros since this is not recorded in
  the XYZ format

**Units**
* Coordinates are in Angstroms.
* The length of a timestep can be set by passing the *delta* argument,
  it's assumed to be in ps (default: 1 ps).

There appears to be no rigid format definition so it is likely users
will need to tweak this Class.

.. _xyz-format:

XYZ File format
---------------

Definiton used by the :class:`XYZReader` and :class:`XYZWriter` (and
the `VMD xyzplugin`_ from whence the definition was taken)::

    [ comment line            ] !! NOT IMPLEMENTED !! DO NOT INCLUDE
    [ N                       ] # of atoms, required by this xyz reader plugin  line 1
    [ molecule name           ] name of molecule (can be blank)                 line 2
    atom1 x y z [optional data] atom name followed by xyz coords                line 3
    atom2 x y z [ ...         ] and (optionally) other data.
    ...
    atomN x y z [ ...         ]                                                 line N+2

* comment lines not implemented
* optional data ignored
* molecule name



.. Links
.. _`VMD xyzplugin`: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html

"""

import os
import errno
import numpy
import itertools

import base
import MDAnalysis
import MDAnalysis.core
import MDAnalysis.core.util as util
from MDAnalysis import NoDataError


class XYZWriter(base.Writer):
    """Writes an XYZ file

    The XYZ file format is not formally defined. This writer follows
    the implement
    http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html .
    """

    format = 'XYZ'
    # these are assumed!
    units = {'time': 'ps', 'length': 'Angstrom'}

    def __init__(self, *args, **kwargs):
        """Initialize the XYZ trajectory writer

        :Arguments:
            *filename*
                file name of trajectory file. If it ends with "gz" then the file
                will be gzip-compressed; if it ends with "bz2" it will be bzip2
                compressed.
         :Keywords:
             *atoms*
                Provide atom names: This can be a list of names or an :class:`AtomGroup`.
                If none is provided, atoms will be called 'X' in the output. These atom
                names will be used when a trajectory is written from raw :class:`Timestep`
                objects which do not contain atom information.

                If you write a :class:`AtomGroup` with :meth:`XYZWriter.write` then atom
                information is taken at each step and *atoms* is ignored.
             *remark*
                single line of text ("molecule name")
        """
        # numatoms is ignored ...
        self.filename = args[0]
        # convert length and time to base units on the fly?
        convert_units = kwargs.pop('convert_units', None)
        self.convert_units = convert_units if convert_units is not None else MDAnalysis.core.flags['convert_lengths']
        self.atomnames = self._get_atomnames(kwargs.pop('atoms', "X"))
        self.remark = kwargs.pop('remark',
                                 "Written by {0} (release {1})".format(self.__class__.__name__, MDAnalysis.__version__))

        self.xyz = util.anyopen(self.filename, 'w')  # can also be gz, bz2

    def _get_atomnames(self, atoms):
        """Return a list of atom names"""
        # AtomGroup
        try:
            return atoms.names()
        except AttributeError:
            pass
        # universe?
        try:
            return atoms.atoms.names()
        except AttributeError:
            pass
        # list or string (can be a single atom name... deal with this in write_next_timestep() once we know numatoms)
        return numpy.asarray(util.asiterable(atoms))

    def close(self):
        """Close the trajectory file and finalize the writing"""
        if self.xyz is not None:
            self.xyz.write("\n")
            self.xyz.close()
        self.xyz = None

    def write(self, obj):
        """Write object *obj* at current trajectory frame to file.

        *obj* can be a :class:`~MDAnalysis.core.AtomGroup.AtomGroup`)
        or a whole :class:`~MDAnalysis.core.AtomGroup.Universe`.

        Atom names in the output are taken from the *obj* or default
        to the value of the *atoms* keyword supplied to the
        :class:`XYZWriter` constructor.

        :Arguments:
          *obj*
            :class:`~MDAnalysis.core.AtomGroup.AtomGroup` or
            :class:`~MDAnalysis.core.AtomGroup.Universe`
        """
        # prepare the Timestep and extract atom names if possible
        # (The way it is written it should be possible to write
        # trajectories with frames that differ in atom numbers
        # but this is not tested.)
        try:
            atoms = obj.atoms
        except AttributeError:
            atoms = None
        if atoms:  # have a AtomGroup
            if hasattr(obj, 'universe'):
                # For AtomGroup and children (Residue, ResidueGroup, Segment)
                ts_full = obj.universe.trajectory.ts
                if ts_full.numatoms == atoms.numberOfAtoms():
                    ts = ts_full
                else:
                    # Only populate a time step with the selected atoms.
                    ts = ts_full.copy_slice(atoms.indices())
            elif hasattr(obj, 'trajectory'):
                # For Universe only --- get everything
                ts = obj.trajectory.ts
            # update atom names
            self.atomnames = atoms.names()
        else:
            ts = obj

        self.write_next_timestep(ts)

    def write_next_timestep(self, ts=None):
        """Write coordinate information in *ts* to the trajectory"""
        if ts is None:
            if not hasattr(self, "ts"):
                raise NoDataError("XYZWriter: no coordinate data to write to trajectory file")
            else:
                ts = self.ts

        if len(self.atomnames) != ts.numatoms:
            self.atomnames = numpy.array([self.atomnames[0]] * ts.numatoms)

        if self.convert_units:
            coordinates = self.convert_pos_to_native(ts._pos, inplace=False)
        else:
            coordinates = ts._pos

        self.xyz.write("{0:d}\n".format(ts.numatoms))
        self.xyz.write("frame {0}\n".format(ts.frame))
        for atom, (x, y, z) in itertools.izip(self.atomnames, coordinates):
            self.xyz.write("%8s  %10.5f %10.5f %10.5f\n" % (atom, x, y, z))


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
          on the fly and also reads streams via
          :class:`~MDAnalysis.core.util.NamedStream`.

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
    _Timestep = base.Timestep

    def __init__(self, filename, **kwargs):
        self.filename = filename

        # the filename has been parsed to be either be foo.xyz or foo.xyz.bz2 by coordinates::core.py
        # so the last file extension will tell us if it is bzipped or not
        root, ext = os.path.splitext(self.filename)
        self.xyzfile = util.anyopen(self.filename, "r")
        self.compression = ext[1:] if ext[1:] != "xyz" else None

        # note that, like for xtc and trr files, __numatoms and __numframes are used quasi-private variables
        # to prevent the properties being recalculated
        # this is because there is no indexing so the way it measures the number of frames is to read the whole file!
        self.__numatoms = None
        self.__numframes = None

        self.fixed = 0
        self.skip = 1
        self.periodic = False
        self.delta = kwargs.pop("delta", 1.0)  # can set delta manually, default is 1ps (taken from TRJReader)
        self.skip_timestep = 1

        self.ts = self._Timestep(self.numatoms)  # numatoms has sideeffects: read trajectory... (FRAGILE)

        # Read in the first timestep (FRAGILE);
        # FIXME: Positions on frame 0 (whatever that means) instead of 1 (as all other readers do).
        #        Haven't quite figured out where to start with all the self._reopen() etc.
        #        (Also cannot just use seek() or reset() because that would break with urllib2.urlopen() streams)
        self._read_next_timestep()

    @property
    def numatoms(self):
        """number of atoms in a frame"""
        if not self.__numatoms is None:  # return cached value
            return self.__numatoms
        try:
            self.__numatoms = self._read_xyz_natoms(self.filename)
        except IOError:
            return 0
        else:
            return self.__numatoms

    def _read_xyz_natoms(self, filename):
        # this assumes that this is only called once at startup and that the filestream is already open
        # (FRAGILE)
        n = self.xyzfile.readline()
        self.close()
        # need to check type of n
        return int(n)

    @property
    def numframes(self):
        if not self.__numframes is None:  # return cached value
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
        linesPerFrame = self.numatoms + 2
        counter = 0
        # step through the file (assuming xyzfile has an iterator)
        for i in self.xyzfile:
            counter = counter + 1
        self.close()

        # need to check this is an integer!
        numframes = int(counter / linesPerFrame)
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
                ts._x[:] = x  # more efficient to do it this way to avoid re-creating the numpy arrays
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
        if self.xyzfile is not None:
            raise IOError(errno.EALREADY, 'XYZ file already opened', self.filename)

        self.xyzfile = util.anyopen(self.filename, "r")

        # reset ts
        ts = self.ts
        ts.status = 1
        ts.frame = 0
        ts.step = 0
        ts.time = 0
        return self.xyzfile

    def Writer(self, filename, **kwargs):
        """Returns a XYZWriter for *filename* with the same parameters as this XYZ.

        All values can be changed through keyword arguments.

        :Arguments:
          *filename*
              filename of the output DCD trajectory
        :Keywords:
          *atoms*
              names of the atoms (if not taken from atom groups)

        :Returns: :class:`XYZWriter` (see there for more details)
        """
        return XYZWriter(filename, **kwargs)

    def close(self):
        """Close xyz trajectory file if it was open."""
        if self.xyzfile is None:
            return
        self.xyzfile.close()
        self.xyzfile = None
