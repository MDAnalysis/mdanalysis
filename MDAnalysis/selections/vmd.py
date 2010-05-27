# MDAnalysis.selections.vmd
# Copyright (c) 2010 Oliver Beckstein
# Published under the GNU Publice License version 2 (or higher)
"""
VMD selections
==============

Write :class:`MDAnalysis.AtomGroup.AtomGroup` selection to a tcl file
that defines VMD atomselect macros. To be used in VMD like this::

  source macros.vmd
  set sel [atomselect top mdanalysis001]  # use macro 001

"""
from __future__ import with_statement

import MDAnalysis.core.util as util

class SelectionWriter(object):
    format = "tcl"
    ext = "vmd"
    def __init__(self, filename, mode="w", numcols=8):
        """Set up SelectionWriter to write to *filename*.

        :Arguments:
           *filename*
               tcl file
           *mode*
               overwrite ("w") or append ("a") ["w"]
           *numcols*
               number of individual index numbers per line [8]
        """
        self.filename = util.filename(filename,ext=self.ext)
        if not mode in ('a', 'w'):
            raise ValueError("mode must be either 'w' or 'a', not %r" % mode)
        self.mode = mode        
        self.numcols = numcols
        self.number = 0

    def write(self,selection,number=None,name=None,frame=None):
        """Write selection to a VMD tcl file.

        :Arguments:
           *selection*
               a :class:`MDAnalysis.AtomGroup.AtomGroup`
           *number*
               selection will be named "mdanalysis<number>"
               (``None`` auto increments between writes; useful
               when appending) [``None``]
           *name*
               selection will be named *name* (instead of numbered)
               [``None``]        
           *frame*
               write selection of this frame (or the current one if
               ``None`` [``None``]
        """
        u = selection.universe
        if frame is not None:            
            u.trajectory[frame]  # advance to frame
        else:
            try:
                frame = u.trajectory.ts.frame
            except AttributeError:
                frame = 1   # should catch cases when we are analyzing a single PDB (?)
        if name is None:
            if number is None:
                self.number += 1
                number = self.number
            name = "mdanalysis%(number)03d" % vars()
        with file(self.filename, self.mode) as tcl:
            tcl.write("# MDAnalysis selection\n")
            tcl.write("atomselect macro %(name)s {index " % vars())
            for iatom, atom in enumerate(selection.atoms):
                index = atom.number  # VMD index is 0-based (as is MDAnalysis?)
                tcl.write(" %(index)d" % vars())
                if (iatom+1) % self.numcols == 0:
                    tcl.write("\\\n\t")
            tcl.write(" }\n")

        
