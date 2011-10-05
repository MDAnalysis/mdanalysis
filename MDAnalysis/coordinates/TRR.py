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

.. _TRR trajectory format:
   http://www.gromacs.org/Documentation/File_Formats/.trr_File


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
:attr:`~MDAnalysis.coordinates.xdrfile.TRR.Timestep._pos` attribute of
:class:`~MDAnalysis.coordinates.xdrfile.TRR.Timestep`::

   # 'modes' is a mode object with M PCs, similar to a MxNx3 array
   # 'xav' the average coordinates, a Nx3 array for N atoms

   N = len(xav)   # number of atoms, i.e. number of coordinates

   W = Writer('pca.trr', numatoms=N)            # TRR writer
   ts = MDAnalysis.coordinates.TRR.Timestep(N)  # TRR time step

   for frame,mode in enumerate(modes[4:16]):
       ts.lmbda = -1
       if frame<=1:
          ts._pos[:] = xav
       else:
          ts._pos[:] = mode.scaledToNorm(1.).array*10   # nm to angstroms
       ts.frame = frame         # manually change the frame number
       ts.step = frame - 1
       if frame <= 1:
          ts.time = frame-1
       else:
          ts.time = mode.frequency
       W.write(ts)             # converts angstrom to nm for gmx

    W.close()

.. _MMTK: http://dirac.cnrs-orleans.fr/Manuals/MMTK/index.html

.. _`recipe by Ramon Crehuet`: http://code.google.com/p/mdanalysis/issues/detail?id=79

Module reference
----------------

.. automodule:: MDAnalysis.coordinates.xdrfile.TRR
   :members:
"""

from xdrfile.TRR import TRRReader, TRRWriter, Timestep
