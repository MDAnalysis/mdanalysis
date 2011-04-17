# MDAnalysis.selections
# Copyright (c) 2010 Oliver Beckstein
# Published under the GNU Publice License version 2 (or higher)
"""
Selections
==========

Functions to write a :class:`MDAnalysis.AtomGroup.AtomGroup` selection
to a file so that it can be used in another programme.

:mod:`MDAnalysis.selections.vmd`
    VMD_ selections
:mod:`MDAnalysis.selections.pymol`
    PyMol_ selections
:mod:`MDAnalysis.selections.gromacs`
    Gromacs_ selections
:mod:`MDAnalysis.selections.charmm`
    CHARMM_ selections

The :class:`MDAnalysis.selections.base.SelectionWriter` base class and
helper functions are in :mod:`MDAnalysis.selections.base`.
"""

import vmd, pymol, gromacs, charmm
from base import get_writer

# Signature:
#   W = SelectionWriter(filename, **kwargs)
#   W.write(AtomGroup)
_selection_writers = {'vmd': vmd.SelectionWriter,
                      'charmm': charmm.SelectionWriter, 'str': charmm.SelectionWriter,
                      'pymol': pymol.SelectionWriter, 'pml': pymol.SelectionWriter,
                      'gromacs': gromacs.SelectionWriter, 'ndx': gromacs.SelectionWriter,
                      }
