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
"""GSD trajectory reader  --- :mod:`MDAnalysis.coordinates.GSD`
============================================================

Class to read the GSD trajectory, output of `HOOMD-blue`_. The GSD format
specifies both the topology and the trajectory of the particles in the
simulation. The topology is read by the
:class:`~MDAnalysis.topology.GSDParser.GSDParser` class.

The GSD format was developed having in mind the possibility of changing number
of particles, particle types, particle identities and changing topology.
Currently this class has limited functionality, due to the fact that the number
of particles and the topology are kept fixed in most MD simulations. The user
will get an error only if at any time step the number of particles is detected
to be different to the one that was set at the first time step. No check on
changes in particle identity or topology is currently implemented.

.. _`HOOMD-blue`: https://glotzerlab.engin.umich.edu/hoomd-blue/

Classes
-------

.. autoclass:: GSDReader
   :inherited-members:

.. autoclass:: GSDPicklable
   :members:

.. autofunction:: gsd_pickle_open

"""
import numpy as np

try:
    import gsd
    import gsd.fl
except ImportError:
    HAS_GSD = False
    import types

    class MockHOOMDTrajectory:
        pass

    gsd = types.ModuleType("gsd")
    gsd.hoomd = types.ModuleType("hoomd")
    gsd.hoomd.HOOMDTrajectory = MockHOOMDTrajectory
else:
    HAS_GSD = True

from . import base
from MDAnalysis.lib.util import store_init_arguments


class GSDReader(base.ReaderBase):
    """Reader for the GSD format."""

    format = "GSD"
    units = {"time": None, "length": None}

    @store_init_arguments
    def __init__(self, filename, **kwargs):
        """
        Parameters
        ----------
        filename : str
            trajectory filename
        **kwargs : dict
            General reader arguments.


        .. versionadded:: 0.17.0
        .. versionchanged:: 2.0.0
            Now use a picklable :class:`gsd.hoomd.HOOMDTrajectory`--
            :class:`GSDPicklable`
        .. versionchanged:: 2.6.0
           Support for GSD versions below 3.0.1 have been dropped. This
           includes support for schema 1.3.
        """
        if not HAS_GSD:
            errmsg = "GSDReader: To read GSD files, please install gsd"
            raise ImportError(errmsg)

        super(GSDReader, self).__init__(filename, **kwargs)
        self.filename = filename
        self.open_trajectory()
        self.n_atoms = self._file[0].particles.N
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self._read_next_timestep()

    def open_trajectory(self):
        """opens the trajectory file using gsd.hoomd module"""
        self._frame = -1
        self._file = gsd_pickle_open(self.filename, mode="r")

    def close(self):
        """close reader"""
        self._file.file.close()

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        return len(self._file)

    def _reopen(self):
        """reopen trajectory"""
        self.close()
        self.open_trajectory()

    def _read_frame(self, frame):
        try:
            myframe = self._file[frame]
        except IndexError:
            raise IOError from None

        # set frame number
        self._frame = frame

        # sets the Timestep object
        self.ts.frame = frame
        self.ts.data["step"] = myframe.configuration.step

        # set frame box dimensions
        self.ts.dimensions = myframe.configuration.box
        self.ts.dimensions[3:] = np.rad2deg(np.arccos(self.ts.dimensions[3:]))

        # set particle positions
        frame_positions = myframe.particles.position
        n_atoms_now = frame_positions.shape[0]
        if n_atoms_now != self.n_atoms:
            raise ValueError(
                "Frame %d has %d atoms but the initial frame has %d"
                " atoms. MDAnalysis in unable to deal with variable"
                " topology!" % (frame, n_atoms_now, self.n_atoms)
            )
        else:
            self.ts.positions = frame_positions
        return self.ts

    def _read_next_timestep(self):
        """read next frame in trajectory"""
        return self._read_frame(self._frame + 1)


class GSDPicklable(gsd.hoomd.HOOMDTrajectory):
    """Hoomd GSD file object (read-only) that can be pickled.

    This class provides a file-like object (as by :func:`gsd.hoomd.open`,
    namely :class:`gsd.hoodm.HOOMDTrajectory`) that, unlike file objects,
    can be pickled. Only read mode is supported.

    When the file is pickled, filename and mode of :class:`gsd.fl.GSDFile` in
    the file are saved. On unpickling, the file is opened by filename.
    This means that for a successful unpickle, the original file still has to
    be accessible with its filename.

    Note
    ----
    Open hoomd GSD files with `gsd_pickle_open`.
    After pickling, the current frame is reset. `universe.trajectory[i]` has
    to be used to return to its original frame.

    Parameters
    ----------
    file: :class:`gsd.fl.GSDFile`
        File to access.

    Example
    -------
    ::

        gsdfileobj = gsd.fl.open(name=filename,
                                     mode='r',
                                     application='gsd.hoomd '+ gsd.version.version,
                                     schema='hoomd',
                                     schema_version=[1, 3])
        file = GSDPicklable(gsdfileobj)
        file_pickled = pickle.loads(pickle.dumps(file))

    See Also
    ---------
    :func:`MDAnalysis.lib.picklable_file_io.FileIOPicklable`
    :func:`MDAnalysis.lib.picklable_file_io.BufferIOPicklable`
    :func:`MDAnalysis.lib.picklable_file_io.TextIOPicklable`
    :func:`MDAnalysis.lib.picklable_file_io.GzipPicklable`
    :func:`MDAnalysis.lib.picklable_file_io.BZ2Picklable`


    .. versionadded:: 2.0.0
    """

    def __getstate__(self):
        return self.file.name, self.file.mode

    def __setstate__(self, args):
        gsd_version = gsd.version.version
        schema_version = [1, 4]
        gsdfileobj = gsd.fl.open(
            name=args[0],
            mode=args[1],
            application="gsd.hoomd " + gsd_version,
            schema="hoomd",
            schema_version=schema_version,
        )
        self.__init__(gsdfileobj)


def gsd_pickle_open(name: str, mode: str = "r"):
    """Open hoomd schema GSD file with pickle function implemented.

    This function returns a GSDPicklable object. It can be used as a
    context manager, and replace the built-in :func:`gsd.hoomd.open` function
    in read mode that only returns an unpicklable file object.

    Schema version will depend on the version of gsd module.

    Note
    ----
    Can be only used with read mode.

    Parameters
    ----------
    name : str
        a filename given a text or byte string.
    mode: str, optional
        'r':  open for reading.

    Returns
    -------
    stream-like object: GSDPicklable

    Raises
    ------
    ValueError
        if `mode` is not one of the allowed read modes

    Examples
    -------
    open as context manager::

        with gsd_pickle_open('filename') as f:
            line = f.readline()

    open as function::

        f = gsd_pickle_open('filename')
        line = f.readline()
        f.close()

    See Also
    --------
    :func:`MDAnalysis.lib.util.anyopen`
    :func:`MDAnalysis.lib.picklable_file_io.pickle_open`
    :func:`MDAnalysis.lib.picklable_file_io.bz2_pickle_open`
    :func:`MDAnalysis.lib.picklable_file_io.gzip_pickle_open`
    :func:`gsd.hoomd.open`


    .. versionadded:: 2.0.0
    .. versionchanged:: 2.6.0
       Only GSD versions 3.0.1+ are supported. 'rb' mode support
       has been replaced with 'r' mode.
    """
    gsd_version = gsd.version.version
    schema_version = [1, 4]
    if mode != "r":
        raise ValueError("Only read mode 'r' " "files can be pickled.")
    gsdfileobj = gsd.fl.open(
        name=name,
        mode=mode,
        application="gsd.hoomd " + gsd_version,
        schema="hoomd",
        schema_version=schema_version,
    )
    return GSDPicklable(gsdfileobj)
