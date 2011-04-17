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
#     Molecular Dynamics Simulations. J. Comput. Chem. (2011),
#     in press.
#

"""
PyMol selections
=================

Write :class:`MDAnalysis.AtomGroup.AtomGroup` selection to a pml file
that defines PyMol atomselect macros. To be used in PyMol like this::

  @macros.pml

The selection appears.
"""
import base

class SelectionWriter(base.SelectionWriter):
    format = "PyMol"
    ext = "pml"
    continuation = '\\' # quoted backslash!
    commentfmt = "# %s"
    default_numterms = 6

    def _translate(self, atoms, **kwargs):
        # PyMol index is 1-based
        def _index(atom):
            return "index %d" % (atom.number + 1)
        return base.join(atoms, ' |', _index)

    def _write_head(self, out, **kwargs):
        out.write(self.comment("MDAnalysis PyMol selection"))
        out.write("select %(name)s, " % kwargs + self.continuation + '\n')


