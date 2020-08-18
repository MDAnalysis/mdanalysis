# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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

import warnings

try:
    shell = get_ipython()
except NameError:
    shell = None
    warnings.warn("You should be in an interactive python shell (IPython or "
                  "Notebook) to use this properly")

from .. import _FORMATTERS


class _Formattermeta(type):
    """Automatic Formatter registration metaclass

    .. versionadded:: 2.0.0
    """
    # Auto register upon class creation
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            fmt = classdict['format']
        except KeyError:
            pass
        else:
            _FORMATTERS[fmt] = {}


class FormatterBase(metaclass=_Formattermeta):
    """Base class for formatters

    .. versionadded:: 2.0.0
    """

    def add_repr(self, obj, mime, func):
        """Add a new formatter to an object

        Parameters
        ----------
        obj : object
            An MDAnalysis object
        mime : str
            The MIME type to add, e.g. "image/svg+xml"
        func : callable
            A function that returns the data that will be displayed by the
            formatter, e.g. raw SVG text for an SVG representation
        """
        # add formatter
        if shell:
            shell.display_formatter.formatters[mime].for_type(obj, func)
        # register it
        _FORMATTERS[self.format][(obj, mime)] = func

    def reset_repr(self, obj, mime):
        """Drop the IPython formatter of an object

        If no other formatter is available to IPython for the object, it will
        use the object's :meth:`__repr__` method

        Parameters
        ----------
        obj : object
            An MDAnalysis object
        mime : str
            The MIME type to drop, e.g. "image/png"
        """
        if shell:
            shell.display_formatter.formatters[mime].pop(obj)
        _FORMATTERS[self.format].pop((obj, mime))

    def reset_all_repr(self):
        """Reset all registered formatters of the formatter class"""
        for obj, mime in list(_FORMATTERS[self.format].keys()):
            self.reset_repr(obj, mime)
