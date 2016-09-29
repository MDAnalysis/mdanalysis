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
import numpy as np
import warnings
import six
from six.moves import zip

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
    returns a Topology object

    Raises
    ------
    * :exc:`IOError` upon failing to read a topology file
    * :exc:`ValueError` upon failing to make sense of the read data

    .. versionadded:: 0.9.0
    .. versionchanged:: 0.9.2
       Added keyword 'universe' to pass to Atom creation.
    """
    def __init__(self, filename, **kwargs):
        self.filename = filename
        self.kwargs = kwargs

    def parse(self):  # pragma: no cover
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


def _string_diff(array):
    """Drop in replacement for np.diff that works on strings"""
    return np.array([0 if x == y else 1
                     for x, y in zip(array[:-1], array[1:])])


def resid_change_squash(resids, *other_attrs):
    """Squash per atom data to per residue according to changes in resid

    Parameters
    ----------
    resids : numpy ndarray
      The array of values which define the identity of the Residue. When
      this changes a new residue has been started
    other_attrs : multiple numpy ndarrays
      Other attributes which get squashed 

    Returns
    -------
    residx
      The Residue *index* that each Atom gets assigned to. [len(resids)]
    new_resids
      The resid of each new 
    other_attrs
      

    Example
    -------
    resids = np.array([2, 2, 3, 3, 2, 2])
    resnames = np.array(['RsA', 'RsA', 'RsB', 'RsB', 'RsC', 'RsC'])
    segids = np.array(['A', 'A', 'A', 'A', 'B', 'B'])

    residx, new_resids, (new_resnames, new_segids) = resid_change_squash(resids, resnames, segids)

    # Per atom res index
    residx: [0, 0, 1, 1, 2, 2]
    # Per residue record of attribute
    resids: [2, 3, 2]
    resnames: ['RsA', 'RsB', 'RsC']
    # Segids now per residue, can squash again in same function to find per-segment
    segids: ['A', 'A', 'B']
    """
    # TODO: Add error checking to inputs ^_^"

    # 1) Detect where resids change
    try:
        diff = np.diff(resids)
    except TypeError:
        # np.diff uses '-', use slower string version
        diff = _string_diff(resids)
    finally:
        diff = np.nonzero(diff)[0]
    # Number of unique residues
    nres = len(diff) + 1

    # 2) Allocate new arrays
    # Per atom record of what residue they belong to
    residx = np.zeros_like(resids, dtype=np.int)
    # Per residue record of various attributes
    new_resids = np.zeros(nres, dtype=resids.dtype)
    new_others = [np.zeros(nres, dtype=o.dtype) for o in other_attrs]

    # 3) Slice through resids and others to find values
    chops = [None] + list(diff + 1) + [None]
    for i, (x, y) in enumerate(zip(chops[:-1], chops[1:])):
        residx[x:y] = i  # atoms between x & y are in the i'th residue
        new_resids[i] = resids[x:y][0]
        for attr, new in zip(other_attrs, new_others):
            new[i] = attr[x:y][0]  # TODO: Check that x:y is the same
            # Should be the same for self consistency...

    return residx, new_resids, new_others
            
