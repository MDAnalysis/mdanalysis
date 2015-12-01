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
from itertools import izip
import numpy as np
import warnings

from ..coordinates.base import IObase


class TopologyReader(IObase):
    """Base class for topology readers

    All topology readers must:
      * Be initialised with a filename
      * Return a struct dict by calling their parse function

    :Raises:
       * :exc:`IOError` upon failing to read a topology file
       * :exc:`ValueError` upon failing to make sense of the read data

    .. versionadded:: 0.9.0
    .. versionchanged:: 0.9.2
       Added keyword 'universe' to pass to Atom creation.
    """
    def __init__(self, filename, universe=None, **kwargs):
        """Standard arguments for a TopologyReader:

        :Arguments:

           *filename*
               name of the topology file

        :Keywords:
           *universe*
               Supply a Universe to the Parser.  This then passes it to the
               atom instances that are created within parsers.
           *kwargs*
               Other keyword arguments that can vary with the specific format.
               These are stored as self.kwargs

        """
        self.filename = filename
        self._u = universe
        self.kwargs = kwargs

    def parse(self):
        raise NotImplementedError("Override this in each subclass")


def squash_by(resids, *attributes):
    """Produce len(n_residues) arrays from len(n_atoms) arrays

    Groups elements in resids according to unique values and
    then compresses all *attributes* arrays

    Arguments
    ---------
    resids - array of resids
    *attributes - other arrays that need to follow the sorting of resids

    Returns
    -------
    new_resids, new_attributes - arrays that are now per-residue
    """
    unique_resids = set(resids)
    n_residues = len(unique_resids)

    # Array of my new resids
    new_resids = np.asarray(np.sort(list(unique_resids)), dtype=np.int32)

    new_order = np.argsort(resids)
    sorted_resids = resids[new_order]
    # Find borders in the sorted resids
    changes = np.where(np.diff(sorted_resids) != 0)[0]
    borders = np.concatenate([[0], changes + 1, [len(resids) + 1]])
    # Sort original attributes according to the new order
    sorted_attributes = [att[new_order] for att in attributes]
    
    # new array for each attribute, keep original dtype
    new_atts = [np.zeros(n_residues, dtype=att.dtype) for att in attributes]

    for att, new_att in izip(attributes, new_atts):
        for i, (x, y) in enumerate(izip(borders[:-1], borders[1:])):
            view = att[x:y]
            if not np.unique(view).size == 1:
                warnings.warn("Nonconsistent resid attribute")
            new_att[i] = view[0]
        
    return new_resids, new_atts
