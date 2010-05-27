# MDAnalysis.selections
# Copyright (c) 2010 Oliver Beckstein
# Published under the GNU Publice License version 2 (or higher)
"""
Selections
==========

Functions to write a :class:`MDAnalysis.AtomGroup.AtomGroup` selection
to a file so that it can be used in another programme such as VMD,
pymol, gromacs, or CHARMM.
"""

import vmd

# frame writers: export to single frame formats such as PDB, gro, crd
# Signature:
#   W = SelectionWriter(filename)
#   W.write(AtomGroup)
_selection_writers = {'vmd': vmd.SelectionWriter,
                      #'str': charmm.SelectionWriter,
                      #'pml': pymol.SelectionWriter,
                      #'ndx': gromacs.SelectionWriter,
                      }
