# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""
Base topology reader classes --- :mod:`MDAnalysis.topology.base`
================================================================

Derive topology reader classes from the base class in this module. All
topology readers raise :exc:`IOError` upon failing to read a topology
file and :exc:`ValueError` upon failing to make sense of the read data.

Classes
-------

.. autoclass:: TopologyReader
   :members:
   :inherited-members:

"""
import six

from . import _PARSERS
from ..coordinates.base import IObase
from ..lib import util


class _Topologymeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            fmt = util.asiterable(classdict['format'])
        except KeyError:
            pass
        else:
            for f in fmt:
                _PARSERS[f] = cls


class TopologyReader(six.with_metaclass(_Topologymeta, IObase)):
    """Base class for topology readers

    Parameters
    ----------
    filename : str
        name of the topology file
    universe : Universe, optional
        Supply a Universe to the Parser.  This then passes it to the
        atom instances that are created within parsers.
    kwargs : optional
        Other keyword arguments that can vary with the specific format.
        These are stored as self.kwargs

    All topology readers must define a `parse` method which
    returns a topology dict.

    Raises
    ------
    * :exc:`IOError` upon failing to read a topology file
    * :exc:`ValueError` upon failing to make sense of the read data

    .. versionadded:: 0.9.0
    .. versionchanged:: 0.9.2
       Added keyword 'universe' to pass to Atom creation.
    """
    def __init__(self, filename, universe=None, **kwargs):
        self.filename = filename
        self._u = universe
        self.kwargs = kwargs

    def parse(self):
        raise NotImplementedError("Override this in each subclass")
