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

"""Reading of `Gromacs XTC trajectories`_.

.. _Gromacs XTC trajectories: http://www.gromacs.org/Documentation/File_Formats/.xtc_File
.. _Gromacs: http://www.gromacs.org


.. SeeAlso:: :mod:`MDAnalysis.coordinates.xdrfile.libxdrfile` for low-level
   bindings to the Gromacs trajectory file formats

Classes
-------

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

import core

class Timestep(core.Timestep):
    """Timestep for a Gromacs XTC trajectory."""

class XTCWriter(core.TrjWriter):
    """Write a Gromacs_ XTC trajectory."""
    format = "XTC"

class XTCReader(core.TrjReader):
    """Read Gromacs_ XTC trajectory."""
    format = "XTC"
    _Timestep = Timestep
    _Writer = XTCWriter

