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


"""XYZ trajectory reader --- :mod:`MDAnalysis.coordinates.XYZ`
==============================================================

The :ref:`XYZ format <xyz-format>` is a loosely defined, simple
coordinate trajectory format. The implemented format definition was
taken from the `VMD xyzplugin`_ and is therefore compatible with VMD.

Note the following:

* Comments are not allowed in the XYZ file (we neither read nor write
  them to remain compatible with VMD).
* The atom name (first column) is ignored during reading.
* The coordinates are assumed to be space-delimited rather than fixed
  width (this may cause issues - see below).
* All fields to the right of the z-coordinate are ignored.
* The unitcell information is all zeros since this is not recorded in
  the XYZ format.

.. rubric:: Units

* Coordinates are in Angstroms.
* The length of a timestep can be set by passing the *dt* argument,
  it's assumed to be in ps (default: 1 ps).

There appears to be no rigid format definition so it is likely users
will need to tweak this class.

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

.. Note::
   * comment lines not implemented (do not include them)
   * molecule name: the line is required but the content is ignored
     at the moment
   * optional data (after the coordinates) are presently ignored


.. Links
.. _`VMD xyzplugin`: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html

"""

from six.moves import range, zip
import os
import errno
import numpy as np
import logging
logger = logging.getLogger('MDAnalysis.analysis.psa')

from . import base
from ..core import flags
from ..lib import util
from ..lib.util import cached
from ..exceptions import NoDataError
from ..version import __version__


class XYZWriter(base.Writer):
    """Writes an XYZ file

    The XYZ file format is not formally defined. This writer follows
    the VMD implementation for the molfile `xyzplugin`_.

    .. _xyzplugin:
       http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html
    """

    format = 'XYZ'
    multiframe = True
    # these are assumed!
    units = {'time': 'ps', 'length': 'Angstrom'}

    def __init__(self, filename, n_atoms=None, atoms='X', convert_units=None,
                 remark='default', **kwargs):
        """Initialize the XYZ trajectory writer

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
            single line of text ("molecule name"). By default write MDAnalysis
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
        self.remark = default_remark if remark == 'default' else remark
        # can also be gz, bz2
        self._xyz = util.anyopen(self.filename, 'wt')

    def _get_atomnames(self, atoms):
        """Return a list of atom names"""
        # AtomGroup
        try:
            return atoms.names
        except AttributeError:
            pass
        # universe?
        try:
            return atoms.atoms.names
        except AttributeError:
            pass
        # list or string (can be a single atom name... deal with this in
        # write_next_timestep() once we know n_atoms)
        if self.n_atoms is None:
            return np.asarray(util.asiterable(atoms))
        if isinstance(atoms, list):
            return atoms
        else:
            return np.asarray([atoms for _ in range(self.n_atoms)])

    def close(self):
        """Close the trajectory file and finalize the writing"""
        if self._xyz is not None:
            self._xyz.write("\n")
            self._xyz.close()
        self._xyz = None

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
                if ts_full.n_atoms == atoms.n_atoms:
                    ts = ts_full
                else:
                    # Only populate a time step with the selected atoms.
                    ts = ts_full.copy_slice(atoms.indices)
            elif hasattr(obj, 'trajectory'):
                # For Universe only --- get everything
                ts = obj.trajectory.ts
            # update atom names
            self.atomnames = atoms.names
        else:
            if isinstance(obj, base.Timestep):
                ts = obj
            else:
                raise TypeError("No Timestep found in obj argument")

        self.write_next_timestep(ts)

    def write_next_timestep(self, ts=None):
        """Write coordinate information in *ts* to the trajectory"""
        if ts is None:
            if not hasattr(self, 'ts'):
                raise NoDataError('XYZWriter: no coordinate data to write to '
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
            if len(self.atomnames) != ts.n_atoms:
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

        self._xyz.write("{0:d}\n".format(ts.n_atoms))
        self._xyz.write("frame {0}\n".format(ts.frame))
        for atom, (x, y, z) in zip(self.atomnames, coordinates):
            self._xyz.write("{0!s:>8}  {1:10.5f} {2:10.5f} {3:10.5f}\n"
                            "".format(atom, x, y, z))


class XYZReader(base.Reader):
    """Reads from an XYZ file

    :Data:
        :attr:`ts`
          Timestep object containing coordinates of current frame

    :Methods:
        ``len(xyz)``
          return number of frames in xyz
        ``for ts in xyz:``
          iterate through trajectory

    .. Note: this can read both compressed (foo.xyz) and compressed
          (foo.xyz.bz2 or foo.xyz.gz) files; uncompression is handled
          on the fly and also reads streams via
          :class:`~MDAnalysis.lib.util.NamedStream`.

    The XYZ file format follows VMD's xyzplugin_ and is also described
    under :ref:`XYZ format <xyz-format>`.

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based. Added *dt* and
       *time_offset* keywords (passed to :class:`Timestep`)
    """

    # Phil Fowler:
    # Validation: the geometric centre of 1284 atoms was calculated over
    # 500 frames using both MDAnalysis and a VMD Tcl script. There was
    # exact agreement (measured to three decimal places). bzipped and
    # gzipped versions of the XYZ file were also tested

    format = "XYZ"
    # these are assumed!
    units = {'time': 'ps', 'length': 'Angstrom'}
    _Timestep = base.Timestep

    def __init__(self, filename, **kwargs):
        super(XYZReader, self).__init__(filename, **kwargs)

        # the filename has been parsed to be either be foo.xyz or foo.xyz.bz2 by
        # coordinates::core.py so the last file extension will tell us if it is
        # bzipped or not
        root, ext = os.path.splitext(self.filename)
        self.xyzfile = util.anyopen(self.filename, "r")
        self.compression = ext[1:] if ext[1:] != "xyz" else None
        self._cache = dict()

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
        with util.anyopen(self.filename, 'r') as f:
            n = f.readline()
        # need to check type of n
        return int(n)

    @property
    @cached('n_frames')
    def n_frames(self):
        try:
            return self._read_xyz_n_frames()
        except IOError:
            return 0

    def _read_xyz_n_frames(self):
        # the number of lines in the XYZ file will be 2 greater than the
        # number of atoms
        linesPerFrame = self.n_atoms + 2
        counter = 0
        offsets = []

        with util.anyopen(self.filename, 'r') as f:
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
        self.xyzfile.seek(self._offsets[frame])
        self.ts.frame = frame - 1  # gets +1'd in next
        return self._read_next_timestep()

    def _read_next_timestep(self, ts=None):
        # check that the timestep object exists
        if ts is None:
            ts = self.ts

        f = self.xyzfile

        try:
            # we assume that there are only two header lines per frame
            f.readline()
            f.readline()
            for i in range(self.n_atoms):
                self.ts._pos[i] = list(map(float, f.readline().split()[1:4]))
            ts.frame += 1
            return ts
        except (ValueError, IndexError) as err:
            raise EOFError(err)

    def rewind(self):
        """reposition on first frame"""
        self._reopen()
        # the next method calls _read_next_timestep
        self.next()

    def _reopen(self):
        self.close()
        self.open_trajectory()

    def open_trajectory(self):
        if self.xyzfile is not None:
            raise IOError(
                errno.EALREADY, 'XYZ file already opened', self.filename)

        self.xyzfile = util.anyopen(self.filename, "r")

        # reset ts
        ts = self.ts
        ts.frame = -1

        return self.xyzfile

    def Writer(self, filename, n_atoms=None, **kwargs):
        """Returns a XYZWriter for *filename* with the same parameters as this
        XYZ.

        Parameters
        ----------
        filename: str
            filename of the output trajectory
        n_atoms: int (optional)
            number of atoms. If none is given use the same number of atoms from
            the reader instance is used
        **kwargs:
            See :class:`XYZWriter` for additional kwargs

        Returns
        -------
        :class:`XYZWriter` (see there for more details)

        See Also
        --------
        :class: `XYZWriter`
        """
        if n_atoms is None:
            n_atoms = self.n_atoms
        return XYZWriter(filename, n_atoms=n_atoms, **kwargs)

    def close(self):
        """Close xyz trajectory file if it was open."""
        if self.xyzfile is None:
            return
        self.xyzfile.close()
        self.xyzfile = None
