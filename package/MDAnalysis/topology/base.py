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
Base topology reader classes --- :mod:`MDAnalysis.topology.base`
================================================================

Derive topology reader classes from the base class in this module. All
topology readers raise :exc:`IOError` upon failing to read a topology
file and :exc:`ValueError` upon failing to make sense of the read data.

Classes
-------

.. autoclass:: TopologyReaderBase
   :members:
   :inherited-members:

"""
from functools import reduce

import itertools
import numpy as np
import warnings

from .. import _PARSERS, _PARSER_HINTS
from ..coordinates.base import IOBase
from ..lib import util



class _Topologymeta(type):
    """Internal: Topology Parser registration voodoo

    When classes which inherit from TopologyReaderBase are *defined*
    this metaclass makes it known to MDAnalysis.  The optional `format`
    attribute and `_format_hint` staticmethod are read:
     - `format` defines the file extension this Parser targets.
     - `_format_hint` defines a function which returns a boolean if the
       Parser can process a particular object

    Eg::

      class ThingParser(TopologyReaderBase):
          format = ['foo', 'bar']

          @staticmethod
          _format_hint(thing):
              try:
                  import WeirdPackage
              except ImportError:
                  return False
              return isinstance(thing, WeirdPackage.Thing)

    This way there is no strict dependency on "WeirdPackage", but if
    a user supplies a WeirdPackage.Thing the "ThingParser' will be able
    to step up and read it.

    .. versionchanged:: 1.0.0
       Added format_hint functionality
    """
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            fmt = util.asiterable(classdict['format'])
        except KeyError:
            pass
        else:
            for fmt_name in fmt:
                fmt_name = fmt_name.upper()
                _PARSERS[fmt_name] = cls

                if '_format_hint' in classdict:
                    _PARSER_HINTS[fmt_name] = classdict['_format_hint'].__func__

class TopologyReaderBase(IOBase, metaclass=_Topologymeta):
    """Base class for topology readers

    Parameters
    ----------
    filename : str
        name of the topology file
    universe : Universe, optional
        Supply a Universe to the Parser.  This then passes it to the
        atom instances that are created within parsers.

    All topology readers must define a `parse` method which
    returns a Topology object

    Raises
    ------
    * :exc:`IOError` upon failing to read a topology file
    * :exc:`ValueError` upon failing to make sense of the read data

    .. versionadded:: 0.9.0
    .. versionchanged:: 0.9.2
       Added keyword 'universe' to pass to Atom creation.
    """
    def __init__(self, filename):
        
        if isinstance(filename, util.NamedStream):
            self.filename = filename
        else:
            self.filename = str(filename)


    def parse(self, **kwargs):  # pragma: no cover
        raise NotImplementedError("Override this in each subclass")


def squash_by(child_parent_ids, *attributes):
    """Squash a child-parent relationship

    Arguments
    ---------
    child_parent_ids - array of ids (unique values that identify the parent)
    *attributes - other arrays that need to follow the sorting of ids

    Returns
    -------
    child_parents_idx - an array of len(child) which points to the index of
                        parent
    parent_ids - len(parent) of the ids
    *parent_attrs - len(parent) of the other attributes
    """
    unique_resids, sort_mask, atom_idx = np.unique(
        child_parent_ids, return_index=True, return_inverse=True)

    return atom_idx, unique_resids, [attr[sort_mask] for attr in attributes]


def change_squash(criteria, to_squash):
    """Squash per atom data to per residue according to changes in resid

    Parameters
    ----------
    criteria : list of numpy ndarray
      Arrays which when changing indicate a new residue
    to_squash : list of numpy arrays
      Arrays which get squashed according to the criteria arrays

    Returns
    -------
    residx : numpy array
      The Residue *index* that each Atom gets assigned to. [len(resids)]
    squashed : numpy array
      The to_squash arrays reduced down to per Residue values


    Example
    -------
    resids = np.array([2, 2, 3, 3, 2, 2])
    resnames = np.array(['RsA', 'RsA', 'RsB', 'RsB', 'RsC', 'RsC'])
    segids = np.array(['A', 'A', 'A', 'A', 'B', 'B'])

    residx, (new_resids, new_resnames, new_segids) = resid_change_squash(
                                                      (resids,), (resids, resnames, segids))

    # Per atom res index
    residx: [0, 0, 1, 1, 2, 2]
    # Per residue record of each attribute
    new_resids: [2, 3, 2]
    new_resnames: ['RsA', 'RsB', 'RsC']
    new_segids: ['A', 'A', 'B']
    """
    def get_borders(*arrays):
        """Generator of indices to slice arrays when they change"""
        borders = np.nonzero(reduce(np.logical_or,
                                    (a[:-1] != a[1:] for a in arrays)))
        # Add Nones so we can slice from start to end
        return [None] + list(borders[0] + 1) + [None]

    l0 = len(criteria[0])
    if not all(len(other) == l0
               for other in itertools.chain(criteria[1:], to_squash)):
        raise ValueError("All arrays must be equally sized")

    # 1) Detect where resids change
    borders = get_borders(*criteria)
    # Number of groups = number of changes + 1
    # 2 `None`s have been added, so -1
    nres = len(borders) - 1

    # 2) Allocate new arrays
    # Per atom record of what residue they belong to
    residx = np.zeros_like(criteria[0], dtype=int)
    # Per residue record of various attributes
    new_others = [np.zeros(nres, dtype=o.dtype) for o in to_squash]

    # 3) Slice through resids and others to find values
    for i, (x, y) in enumerate(zip(borders[:-1], borders[1:])):
        residx[x:y] = i  # atoms between x & y are in the i'th residue
        for old, new in zip(to_squash, new_others):
            new[i] = old[x:y][0]  # TODO: Check that x:y is the same
            # Should be the same for self consistency...

    return residx, new_others


def reduce_singular(values):
    """Returns the value in an array of length 1, or
    the tuple of an array with a longer lengh.

    Parameters
    ----------
    values: array-like
        Array to squash

    Returns
    -------
    values: tuple or single value
    """
    if len(values) == 1:
        return values[0]
    else:
        return tuple(values)
