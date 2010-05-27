# MDAnalysis.selections.charmm
# Part of MDAnalysis http://mdanalysis.googlecode.com
# Copyright (c) 2006-2010 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# Released under the GNU Public Licence, v2
"""
CHARMM selections
=================

Write :class:`MDAnalysis.AtomGroup.AtomGroup` selection to a str file
that defines a CHARMM selection. To be used in CHARMM like this::

  stream macro.str

The selection is name *mdanalysis001*.
"""
import base

class SelectionWriter(base.SelectionWriter):
    format = "CHARMM"
    ext = "str"
    continuation = '-'
    commentfmt = "! %s"
    default_numterms = 4  # be conservative because CHARMM only reads 72 columns

    def _translate(self, atoms, **kwargs):
        # CHARMM index is 1-based
        def _index(atom):
            return "BYNUM %d" % (atom.number + 1)
        return base.join(atoms, ' .or.', _index)

    def _write_head(self, out, **kwargs):
        out.write(self.comment("MDAnalysis CHARMM selection"))
        out.write("DEFINE %(name)s SELECT " % kwargs + self.continuation + '\n')

    def _write_tail(self, out, **kwargs):
        out.write("END")

        
