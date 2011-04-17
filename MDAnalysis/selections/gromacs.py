# MDAnalysis.selections.gromacs
# Part of MDAnalysis http://mdanalysis.googlecode.com
# Copyright (c) 2006-2010 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# Released under the GNU Public Licence, v2
"""
Gromacs selections
==================

Write :class:`MDAnalysis.AtomGroup.AtomGroup` selection to a ndx file
that defines a Gromacs index group. To be used in Gromacs like this::

  <GROMACS_COMMAND> -n macro.ndx

The index groups are named *mdanalysis001*, *mdanalysis002*, etc.
"""
import base

class SelectionWriter(base.SelectionWriter):
    format = "Gromacs"
    ext = "ndx"
    default_numterms = 12

    def _translate(self, atoms, **kwargs):
        # Gromacs index is 1-based; MDAnalysis is 0-based
        return [str(atom.number + 1) for atom in atoms]

    def _write_head(self, out, **kwargs):
        out.write("[ %(name)s ]\n" % kwargs)


