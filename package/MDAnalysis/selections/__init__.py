# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
Selection exporters
===================

Functions to write a :class:`MDAnalysis.core.groups.AtomGroup` selection
to a file so that it can be used in another programme.

:mod:`MDAnalysis.selections.vmd`
    VMD_ selections
:mod:`MDAnalysis.selections.pymol`
    PyMol_ selections
:mod:`MDAnalysis.selections.gromacs`
    Gromacs_ selections
:mod:`MDAnalysis.selections.charmm`
    CHARMM_ selections

The :class:`MDAnalysis.selections.base.SelectionWriterBase` base class and
helper functions are in :mod:`MDAnalysis.selections.base`, with the
exception of `:func:get_writer`:

.. autofunction:: get_writer
"""
import os.path

from .. import _SELECTION_WRITERS

from . import base
from . import vmd
from . import pymol
from . import gromacs
from . import charmm
from . import jmol


def get_writer(filename: str, defaultformat: str) -> base.SelectionWriterBase:
    """Return a SelectionWriter for `filename` or a `defaultformat`.

    Parameters
    ----------
    filename : str
       name of the output file; the extension is used to guess the file format
    defaultformat : str
       if `filename` does not have an extension, use `defaultformat` instead

    Returns
    -------
    SelectionWriterBase : `type`
        the writer *class* for the detected format

    Raises
    ------
    :exc:`NotImplementedError`
        for any format that is not defined
    """
    format = None
    if filename:
        format = os.path.splitext(filename)[1][1:]  # strip initial dot!
    format = format or defaultformat  # use default if no fmt from fn
    format = format.strip().upper()  # canonical for lookup
    try:
        return _SELECTION_WRITERS[format]
    except KeyError:
        errmsg = (f"Writing as {format} is not implemented; only "
                  f"{ _SELECTION_WRITERS.keys()} will work.")
        raise NotImplementedError(errmsg) from None
