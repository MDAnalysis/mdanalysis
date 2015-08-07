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
Gromacs TRR file IO --- :mod:`MDAnalysis.coordinates.TRR`
=========================================================

The Gromacs `TRR trajectory format`_ is a lossless format like
e.g. the DCD format (see :mod:`~MDAnalysis.coordinates.DCD`) and
unlike the :mod:`~MDAnalysis.coordinates.XTC` format, which stores
reduced precision coordinates. Therefore, if one wants to convert
*to* Gromacs trajectories without loss of precision then one should
use the TRR format.

The TRR format can store *velocities* and *forces* in addition to
coordinates. It is also used by other Gromacs tools to store and
process other data such as modes from a principal component analysis.

The TRR I/O interface uses :mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile2`
to implement random access to frames. This works by initially building an
internal index of all frames and then using this index for direct
seeks. Building the index is triggered by
:func:`~MDAnalysis.coordinates.xdrfile.libxdrfile2.read_trr_n_frames`, which
typically happens when one accesses the :attr:`TRRReader.n_frames` attribute
for the first time. Building the index may take many minutes for large
trajectories but afterwards access is faster than with native Gromacs tools.


.. _TRR trajectory format:
   http://www.gromacs.org/Documentation/File_Formats/.trr_File

.. versionchanged:: 0.8.0
   The TRR I/O interface now uses
   :mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile2`, which has seeking and
   indexing capabilities. Note that unlike
   :mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile` before it,
   :mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile2` is distributed under the
   GNU GENERAL PUBLIC LICENSE, version 2 (or higher).
   :class:`~MDAnalysis.coordinates.TRR.Timestep` now correctly deals
   with presence/absence of coordinate/velocity/force information on a
   per-frame basis.


Tips and Tricks
---------------

Filling a TRR with PCA modes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following `recipe by Ramon Crehuet`_ shows how to convert modes
stored in a NumPy-like array (e.g. from a PCA analysis with MMTK_) to
a TRR usable by Gromacs. The idea is to manually fill a
:class:`~MDAnalysis.coordinates.xdrfile.TRR.Timestep` with the desired
values and then write it to a file with the appropriate
:class:`~MDAnalysis.coordinates.xdrfile.TRR.TRRWriter`. In order to
respect the Gromacs format for modes in a TRR file, one must write the
average coordinates in the first frame of the TRR and the modes into
subsequent ones. The mode number is stored in the
:attr:`~MDAnalysis.coordinates.xdrfile.TRR.Timestep.step` attribute
and the mode coordinates are filling the
:attr:`~MDAnalysis.coordinates.xdrfile.TRR.Timestep.positions` attribute of
:class:`~MDAnalysis.coordinates.xdrfile.TRR.Timestep`::

   # 'modes' is a mode object with M PCs, similar to a MxNx3 array
   # 'xav' the average coordinates, a Nx3 array for N atoms

   N = len(xav)   # number of atoms, i.e. number of coordinates

   W = Writer('pca.trr', n_atoms=N)            # TRR writer
   ts = MDAnalysis.coordinates.TRR.Timestep(N)  # TRR time step
                                                # N of atoms is passed.
   for frame,mode in enumerate(modes[4:16]):
       ts.lmbda = -1

       ts.frame = frame         # manually change the frame number
       ts._frame = frame - 1

       if frame<=1:
          ts.positions = xav
       else:
          ts.positions = mode.scaledToNorm(1.).array*10   # nm to angstroms
       if frame <= 1:
          ts.time = frame-1
       else:
          ts.time = mode.frequency
       W.write(ts)             # converts angstrom to nm for gmx

    W.close()

.. _MMTK: http://dirac.cnrs-orleans.fr/Manuals/MMTK/index.html

.. _`recipe by Ramon Crehuet`: https://github.com/MDAnalysis/mdanalysis/issues/79

Module reference
----------------

.. autoclass:: Timestep
   :members:
   :inherited-members:

.. autoclass:: TRRReader
   :members:
   :inherited-members:

.. autoclass:: TRRWriter
   :members:
   :inherited-members:
"""

from .xdrfile.TRR import TRRReader, TRRWriter, Timestep
