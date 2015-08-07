# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
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


"""
Gromacs XTC file IO --- :mod:`MDAnalysis.coordinates.XTC`
=========================================================

The Gromacs `XTC trajectory format`_ is a format with lossy
compression. Coordinates are only stored with a fixed precision (by
default, 1/1000 of a nm). The XTC format can only store
*coordinates*. Its main advantage is that it requires less disk space
than e.g. TRR or DCD trajectories and the loss of precision is usually
not a problem.

If one wants to store Gromacs trajectories without loss of precision
or with velocities and/or forces then one should use the TRR format
(see module :mod:`~MDAnalysis.coordinates.TRR`).

The XTC I/O interface uses
:mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile2` to implement random
access to frames. This works by initially building an internal index
of all frames and then using this index for direct seeks. Building the
index is triggered by
:func:`~MDAnalysis.coordinates.xdrfile.libxdrfile2.read_xtc_n_frames`,
which typically happens when one accesses the
:attr:`XTCReader.n_frames` attribute for the first time. Building the
index may take many minutes for large trajectories but afterwards
access is faster than with native Gromacs tools.

.. _XTC trajectory format:
   http://www.gromacs.org/Documentation/File_Formats/.xtc_File

.. versionchanged:: 0.8.0
   The XTC I/O interface now uses
   :mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile2`, which has
   seeking and indexing capabilities. Note that unlike
   :mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile` before it,
   :mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile2` is distributed
   under the GNU GENERAL PUBLIC LICENSE, version 2 (or higher).
   :class:`~MDAnalysis.coordinates.XTC.Timestep` now correctly
   deals with presence/absence of coordinate/velocity/force
   information on a per-frame basis.

Module reference
----------------

.. autoclass:: Timestep
   :members:
   :inherited-members:

.. autoclass:: XTCReader
   :members:
   :inherited-members:

.. autoclass:: XTCWriter
   :members:
   :inherited-members:

"""

from .xdrfile.XTC import XTCReader, XTCWriter, Timestep
