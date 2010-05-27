# MDAnalysis.selections.vmd
# Part of MDAnalysis http://mdanalysis.googlecode.com
# Copyright (c) 2006-2010 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# Released under the GNU Public Licence, v2
"""
VMD selections
==============

Write :class:`MDAnalysis.AtomGroup.AtomGroup` selection to a tcl file
that defines VMD atomselect macros. To be used in VMD like this::

  source macros.vmd
  set sel [atomselect top mdanalysis001]  # use macro 001

"""
import base

class SelectionWriter(base.SelectionWriter):
    format = "VMD"
    ext = "vmd"
    continuation = '\\'
    commentfmt = "# %s"

    def _write_head(self, out, **kwargs):
        out.write(self.comment("MDAnalysis VMD selection"))
        out.write("atomselect macro %(name)s {index " % kwargs)

    def _translate(self, atoms, **kwargs):
        # VMD index is 0-based (as is MDAnalysis)
        return [str(atom.number) for atom in atoms]

    def _write_tail(self, out, **kwargs):
        out.write("}")

