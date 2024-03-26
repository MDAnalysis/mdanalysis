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

Definition used by the :class:`XYZReader` and :class:`XYZWriter` (and
the `VMD xyzplugin`_ from whence the definition was taken)::

    [ comment line            ] !! NOT IMPLEMENTED !! DO NOT INCLUDE
    [ N                       ] # of atoms, required by this xyz reader plugin  line 1
    [ molecule name           ] name of molecule (can be blank)                 line 2
    atom1 x y z [optional data] atom name followed by xyz coords                line 3
    atom2 x y z [ ...         ] and (optionally) other data.
    ...
    atomN x y z [ ...         ]                                                 line N+2


Note
----
* comment lines not implemented (do not include them)
* molecule name: the line is required but the content is ignored
  at the moment
* optional data (after the coordinates) are presently ignored


.. Links
.. _`VMD xyzplugin`: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html

Classes
-------

"""
import itertools
import os
import errno
import numpy as np
import warnings
import logging
logger = logging.getLogger('MDAnalysis.coordinates.XYZ')

from . import base
from .timestep import Timestep
from ..lib import util
from ..lib.util import cached, store_init_arguments
from ..exceptions import NoDataError
from ..version import __version__


class XYZWriter(base.WriterBase):
    """Writes an XYZ file

    The XYZ file format is not formally defined. This writer follows
    the VMD implementation for the molfile `xyzplugin`_.


    Notes
    -----
    By default, the XYZ writer will attempt to use the input
    :class:`~MDAnalysis.core.groups.AtomGroup` or
    :class:`~MDAnalysis.core.universe.Universe` ``elements`` record to assign
    atom names in the XYZ file. If the ``elements`` record is missing, then
    the ``name`` record will be used. In the event that neither of these are
    available, the atoms will all be named ``X``. Please see, the
    `User Guide`_ for more information on how to add topology attributes if
    you wish to add your own elements / atom names to a
    :class:`~MDAnalysis.core.universe.Universe`.


    .. Links

    .. _xyzplugin:
           http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html
    .. _User Guide:
           https://userguide.mdanalysis.org/examples/constructing_universe.html#Adding-topology-attributes


    .. versionchanged:: 1.0.0
       Use elements attribute instead of names attribute, if present.
    .. versionchanged:: 2.0.0
       Support for passing timestep to the writer was deprecated in 1.0 and
       has now been removed. As a consequence, custom names can no longer be
       passed to the writer, these should be added to the
       :class:`~MDAnalysis.core.universe.Universe`, or
       :class:`~MDAnalysis.core.groups.AtomGroup` before invoking the writer.
    """

    format = 'XYZ'
    multiframe = True
    # these are assumed!
    units = {'time': 'ps', 'length': 'Angstrom'}

    def __init__(self, filename, n_atoms=None, convert_units=True,
                 remark=None, **kwargs):
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
        convert_units : bool (optional)
            convert quantities to default MDAnalysis units of Angstrom upon
            writing  [``True``]
        remark: str (optional)
            single line of text ("molecule name"). By default writes MDAnalysis
            version and frame


        .. versionchanged:: 1.0.0
           Removed :code:`default_remark` variable (Issue #2692).
        .. versionchanged:: 2.0.0
           Due to the removal of timestep as an input for writing, the atoms
           parameter is no longer relevant and has been removed. If passing
           an empty universe, please use ``add_TopologyAttr`` to add in the
           required elements or names.
        """
        self.filename = filename
        self.remark = remark
        self.n_atoms = n_atoms
        self.convert_units = convert_units

        # can also be gz, bz2
        self._xyz = util.anyopen(self.filename, 'wt')

    def _get_atoms_elements_or_names(self, atoms):
        """Return a list of atom elements (if present) or fallback to atom names"""
        try:
            return atoms.atoms.elements
        except (AttributeError, NoDataError):
            try:
                return atoms.atoms.names
            except (AttributeError, NoDataError):
                wmsg = ("Input AtomGroup or Universe does not have atom "
                        "elements or names attributes, writer will default "
                        "atom names to 'X'")
                warnings.warn(wmsg)
                return itertools.cycle(('X',))

    def close(self):
        """Close the trajectory file and finalize the writing"""
        if self._xyz is not None:
            self._xyz.write("\n")
            self._xyz.close()
        self._xyz = None

    def write(self, obj):
        """Write object `obj` at current trajectory frame to file.

        Atom elements (or names) in the output are taken from the `obj` or
        default to the value of the `atoms` keyword supplied to the
        :class:`XYZWriter` constructor.

        Parameters
        ----------
        obj : Universe or AtomGroup
            The :class:`~MDAnalysis.core.groups.AtomGroup` or
            :class:`~MDAnalysis.core.universe.Universe` to write.


        .. versionchanged:: 2.0.0
           Deprecated support for Timestep argument has now been removed.
           Use AtomGroup or Universe as an input instead.
        """
        # prepare the Timestep and extract atom names if possible
        # (The way it is written it should be possible to write
        # trajectories with frames that differ in atom numbers
        # but this is not tested.)
        try:
            atoms = obj.atoms
        except AttributeError:
            errmsg = "Input obj is neither an AtomGroup or Universe"
            raise TypeError(errmsg) from None
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
            self.atomnames = self._get_atoms_elements_or_names(atoms)

        self._write_next_frame(ts)

    def _write_next_frame(self, ts=None):
        """
        Write coordinate information in *ts* to the trajectory

        .. versionchanged:: 1.0.0
           Print out :code:`remark` if present, otherwise use generic one 
           (Issue #2692).
           Renamed from `write_next_timestep` to `_write_next_frame`.
        """
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
            if (not isinstance(self.atomnames, itertools.cycle) and
                len(self.atomnames) != ts.n_atoms):
                logger.info('Trying to write a TimeStep with unknown atoms. '
                            'Expected {} atoms, got {}. Try using "write" if you are '
                            'using "_write_next_frame" directly'.format(
                                len(self.atomnames), ts.n_atoms))
                self.atomnames = np.array([self.atomnames[0]] * ts.n_atoms)

        if self.convert_units:
            coordinates = self.convert_pos_to_native(
                ts.positions, inplace=False)
        else:
            coordinates = ts.positions

        # Write number of atoms
        self._xyz.write("{0:d}\n".format(ts.n_atoms))

        # Write remark
        if self.remark is None:
            remark = "frame {} | Written by MDAnalysis {} (release {})\n".format(
                ts.frame, self.__class__.__name__, __version__)

            self._xyz.write(remark)
        else:
            self._xyz.write(self.remark.strip() + "\n")

        # Write content
        for atom, (x, y, z) in zip(self.atomnames, coordinates):
            self._xyz.write("{0!s:>8}  {1:10.5f} {2:10.5f} {3:10.5f}\n"
                            "".format(atom, x, y, z))


class XYZReader(base.ReaderBase):
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
    _Timestep = Timestep

    @store_init_arguments
    def __init__(self, filename, **kwargs):
        super(XYZReader, self).__init__(filename, **kwargs)

        # the filename has been parsed to be either be foo.xyz or foo.xyz.bz2 by
        # coordinates::core.py so the last file extension will tell us if it is
        # bzipped or not
        root, ext = os.path.splitext(self.filename)
        self.xyzfile = util.anyopen(self.filename)
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
        with util.anyopen(self.filename) as f:
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
        self.xyzfile.seek(self._offsets[frame])
        self.ts.frame = frame - 1  # gets +1'd in next
        return self._read_next_timestep()

    def _read_next_timestep(self, ts=None):
        # check that the timestep object exists
        if ts:
            warnings.warn("ts argument to _read_next_timestep is deprecated as of 2.7.0 and will be removed in 3.0.0, see #3928")

        if ts is None:
            ts = self.ts

        f = self.xyzfile

        try:
            # we assume that there are only two header lines per frame
            f.readline()
            f.readline()
            # convert all entries at the end once for optimal speed
            tmp_buf = []
            for i in range(self.n_atoms):
                tmp_buf.append(f.readline().split()[1:4])
            ts.positions = tmp_buf
            ts.frame += 1
            return ts
        except (ValueError, IndexError) as err:
            raise EOFError(err) from None

    def _reopen(self):
        self.close()
        self.open_trajectory()

    def open_trajectory(self):
        if self.xyzfile is not None:
            raise IOError(
                errno.EALREADY, 'XYZ file already opened', self.filename)

        self.xyzfile = util.anyopen(self.filename)

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
        :class:`XYZWriter`
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
