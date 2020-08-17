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
Calculating RDKit descriptors and fingerprints --- :mod:`MDAnalysis.analysis.RDKit`
==========================================================================

This module contains a wrapper class to calculate descriptors and fingerprints
from RDKit molecules. Technically, it is not limited to descriptors provided
by RDKit since any function that can take an RDKit molecule as input will work.


.. versionadded:: 2.0.0

"""

import warnings

import numpy as np
try:
    from rdkit.Chem import (AllChem, DataStructs, rdMolDescriptors, 
                            Descriptors, Descriptors3D, EState, 
                            GraphDescriptors, Lipinski, MolSurf)
except ImportError:
    raise ImportError("RDKit is needed to use the RDKit descriptors code, but "
                      "it doesn't appear to be installed")

from .base import AnalysisBase


_RDKIT_FP = {
    "AtomPair": AllChem.GetAtomPairFingerprint,
    "hashed_AtomPair": AllChem.GetHashedAtomPairFingerprint,
    "Morgan": AllChem.GetMorganFingerprint,
    "hashed_Morgan": AllChem.GetHashedMorganFingerprint,
    "TopologicalTorsion": AllChem.GetTopologicalTorsionFingerprint,
    "hashed_TopologicalTorsion": AllChem.GetHashedTopologicalTorsionFingerprint,
    "RDKit": AllChem.UnfoldedRDKFingerprintCountBased,
    "hashed_RDKit": AllChem.RDKFingerprint,
    "MACCSKeys": AllChem.GetMACCSKeysFingerprint,
}


_RDKIT_DESCRIPTORS = {}
for module in [rdMolDescriptors, Descriptors, Descriptors3D, EState, 
    EState.EState_VSA, GraphDescriptors, Lipinski, MolSurf]:
    for arg in [x for x in dir(module)
                if not (x.startswith("_") or x.endswith("_"))]:
        thing = getattr(module, arg)
        if callable(thing):
            _RDKIT_DESCRIPTORS[arg] = thing


def get_fingerprint(ag, kind, hashed=True, as_array=True, **kwargs):
    """Generate a fingerprint for an AtomGroup using RDKit

    Parameters
    ----------
    ag : MDAnalysis.core.groups.AtomGroup
        The AtomGroup used to generate the fingerprint
    kind : str
        The kind of fingerprint to generate. One of:
        * AtomPair
        * Morgan
        * TopologicalTorsion
        * RDKit
        * MACCSKeys
    hashed : bool
        Return a hashed version of the fingerprint
    as_array : bool
        Return the fingerprint as a numpy array
    \**kwargs : object
        Arguments passed to the fingerprint function

    Returns
    -------
    fp : numpy.array or rdkit.DataStructs.cDataStructs object
        The fingerprint in the desired shape and format

    Notes
    -----
    To generate a Morgan fingerprint, don't forget to specify the radius::

        get_fingerprint(ag, 'Morgan', radius=2)
    
    """
    kind = f"hashed_{kind}" if hashed else kind
    try:
        fp_function = _RDKIT_FP[kind]
    except KeyError:
        if kind == "hashed_MACCSKeys":
            raise KeyError(f"MACCSKeys is not available in a hashed version. "
                           "Please use `hashed=False`") from None
        raise KeyError(f"Could not find {kind!r} in the available "
                       "fingerprints") from None
    mol = ag.convert_to("RDKIT")
    fp = fp_function(mol, **kwargs)
    if not as_array:
        return fp
    array = np.zeros((0, ), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, array)
    return array


class RDKitDescriptors(AnalysisBase):
    """TODO"""
    def __init__(self, atomgroup, *args, **kwargs):
        super().__init__(atomgroup.universe.trajectory,
                         **kwargs)
        self._ag = atomgroup
        self._functions = {}
        for thing in args:
            if callable(thing):
                self._functions[thing.__name__] = thing
            else:
                func = _RDKIT_DESCRIPTORS.get(thing, 
                    _RDKIT_DESCRIPTORS.get(f"Calc{thing}"))
                if func is None:
                    raise KeyError(f"Could not find {thing!r} in RDKit. "
                        "Please try passing the required function directly "
                        "instead of a string")
                self._functions[thing] = func

    def _prepare(self):
        self.results = {name: None for name in self._functions.keys()}
        self._mol = self._ag.convert_to("RDKIT")

    def _single_frame(self):
        for name, func in self._functions.items():
            self.results[name] = func(self._mol)
