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
===================================================================================

This module contains a wrapper class to calculate descriptors and fingerprints
from RDKit molecules. Technically, it is not limited to descriptors provided
by RDKit since any function that can take an RDKit molecule as input will work.


.. versionadded:: 2.0.0

"""

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
            _RDKIT_DESCRIPTORS[arg] = (module, thing)


def get_fingerprint(ag, kind, hashed=False, dtype=None, **kwargs):
    r"""Generate a fingerprint for an AtomGroup using RDKit

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
    dtype : None or str
        Return type for the fingerprint: One of:

        * None for the original RDKit object
        * "dict" for a dictionary containing the bits that are set
        * "array" for the sequence of on and off bits

    **kwargs : object
        Arguments passed to the fingerprint function

    Returns
    -------
    fp : rdkit.DataStructs.cDataStructs object or dict or numpy.ndarray
        The fingerprint in the desired shape and format

    Notes
    -----
    To generate a Morgan fingerprint, don't forget to specify the radius::

        get_fingerprint(ag, 'Morgan', radius=2)

    ``dtype="array"`` is not recommended for non-hashed fingerprints
    (except MACCS) as it is likely to overflow the available memory and crash.

    """
    key = f"hashed_{kind}" if hashed else kind
    try:
        fp_function = _RDKIT_FP[key]
    except KeyError:
        if key == "hashed_MACCSKeys":
            raise ValueError(f"MACCSKeys is not available in a hashed version."
                             " Please use `hashed=False`") from None
        raise ValueError(f"Could not find {kind!r} in the available "
                         "fingerprints") from None
    mol = ag.convert_to("RDKIT")
    fp = fp_function(mol, **kwargs)
    if dtype is None:
        return fp
    elif dtype == "array":
        array = np.zeros((0, ), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, array)
        return array
    elif dtype == "dict":
        classname = fp.__class__.__name__
        if classname.endswith("BitVect"):
            return {bit: 1 for bit in fp.GetOnBits()}
        return fp.GetNonzeroElements()
    else:
        raise ValueError(f"{dtype!r} is not a supported output type")


class RDKitDescriptors(AnalysisBase):
    r"""Class to compute molecular descriptors through RDKit

    Parameters
    ----------
    atomgroup : MDAnalysis.core.groups.AtomGroup
        The AtomGroup used to calculate descriptors
    *descriptors : str or function
        Either a descriptor name in RDKit or a function that takes an RDKit
        molecule as argument

    Attributes
    ----------
    n_descriptors : int
        Number of descriptors requested
    results : numpy.ndarray
        2D array of shape (n_frames, n_descriptors)

    Example
    -------
    Here's an example with a custom function for a trajectory with 3 frames::

        >>> def num_atoms(mol):
        ...    return mol.GetNumAtoms()
        >>> desc = RDKitDescriptors(u.atoms,
        ...                         "MolWt", "RadiusOfGyration", num_atoms)
        >>> desc.run()
        >>> desc.results
        array([[46.06900000000002, 1.161278342193013, 'C2H6O', 9],
               [46.06900000000002, 1.175492972121405, 'C2H6O', 9],
               [46.06900000000002, 1.173230936577319, 'C2H6O', 9]],
              dtype=object)

    """

    def __init__(self, atomgroup, *descriptors, **kwargs):
        super().__init__(atomgroup.universe.trajectory,
                         **kwargs)
        self._ag = atomgroup
        self._functions = []
        for thing in descriptors:
            if callable(thing):
                self._functions.append(thing)
            else:
                obj = _RDKIT_DESCRIPTORS.get(
                    thing, _RDKIT_DESCRIPTORS.get(f"Calc{thing}"))
                if obj is None:
                    raise ValueError(
                        f"Could not find {thing!r} in RDKit. "
                        "Please try passing the required function directly "
                        "instead of a string")
                module, func = obj
                self._functions.append(func)
        self.n_descriptors = len(self._functions)

    def _prepare(self):
        self.results = np.empty(shape=(self.n_frames, self.n_descriptors),
                                dtype=object)

    def _single_frame(self):
        mol = self._ag.convert_to("RDKIT")
        for i, func in enumerate(self._functions):
            self.results[self._frame_index][i] = func(mol)

    @staticmethod
    def list_available(flat=False):
        """List the available RDKit functions to calculate descriptors

        Parameters
        ----------
        flat : bool
            Return a flat array instead of a dictionary

        Returns
        -------
        descriptors : dict or numpy.ndarray
            A dictionary of list of descriptors indexed by the RDKit module
            where the descriptor can be found, or a flat array of descriptors

        Notes
        -----
        This is NOT a curated list of descriptors, some of the names listed
        here might correspond to helper functions available in RDKit that do
        not represent any kind of molecular property.

        """
        if flat:
            return np.array(list(_RDKIT_DESCRIPTORS.keys()), dtype=object)
        t = [(module.__name__, desc)
             for desc, (module, func) in _RDKIT_DESCRIPTORS.items()]
        descriptors = {}
        for module, desc in sorted(t):
            try:
                descriptors[module].append(desc)
            except KeyError:
                descriptors[module] = [desc]
        return descriptors
