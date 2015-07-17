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
Selection exporters
===================

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
helper functions are in :mod:`MDAnalysis.selections.base`, with the
exception of `:func:get_writer`:

.. autofunction:: get_writer
"""
from __future__ import absolute_import

import os.path

from . import vmd
from . import pymol
from . import gromacs
from . import charmm
from . import jmol

# Signature:
#   W = SelectionWriter(filename, **kwargs)
#   W.write(AtomGroup)
_selection_writers = {
    'vmd': vmd.SelectionWriter,
    'charmm': charmm.SelectionWriter, 'str': charmm.SelectionWriter,
    'pymol': pymol.SelectionWriter, 'pml': pymol.SelectionWriter,
    'gromacs': gromacs.SelectionWriter, 'ndx': gromacs.SelectionWriter,
    'jmol': jmol.SelectionWriter, 'spt': jmol.SelectionWriter
}

def get_writer(filename, defaultformat):
    """Return a :class:`SelectionWriter` for *filename* or a *defaultformat*."""

    if filename:
        format = os.path.splitext(filename)[1][1:]  # strip initial dot!
    format = format or defaultformat  # use default if no fmt from fn
    format = format.strip().lower()  # canonical for lookup
    try:
        return _selection_writers[format]
    except KeyError:
        raise NotImplementedError("Writing as %r is not implemented; only %r will work."
                                  % (format, _selection_writers.keys()))
