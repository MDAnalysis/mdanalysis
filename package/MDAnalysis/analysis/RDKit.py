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

:Author: Cédric Bouysset
:Year: 2020 
:Copyright: GNU Public License, v2 or any higher version 

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
from ..core.groups import AtomGroup


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
for module in [rdMolDescriptors, Descriptors, Descriptors3D, EState.EState,
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
        Arguments passed to the fingerprint function, i.e. ``nBits=2048`` for
        hashed fingerprints, ``radius=2`` for Morgan fingerprints...etc.

    Returns
    -------
    fp : rdkit.DataStructs.cDataStructs object or dict or numpy.ndarray
        The fingerprint in the desired shape and format

    Examples
    --------
    Hashed fingerprint::

        >>> fp = get_fingerprint(u.atoms, 'AtomPair', hashed=True)
        >>> fp.GetNonzeroElements()
        {619: 1, 745: 5, 938: 6, 1060: 1, 1130: 1, 1321: 1, 1593: 3, 1648: 8,
         1658: 2, 1723: 1, 1843: 4, 1969: 3}

    Non-hashed fingerprint::

        >>> fp = get_fingerprint(u.atoms, 'MACCSKeys', dtype="array")
        >>> np.where(fp == 1)[0]
        array([ 82, 109, 112, 114, 126, 138, 139, 153, 155, 157, 160, 164])

    Notes
    -----
    Links to the documentation of the available fingerprints:

        * AtomPair: :func:`~rdkit.Chem.rdMolDescriptors.GetAtomPairFingerprint`
        * Hashed AtomPair: :func:`~rdkit.Chem.rdMolDescriptors.GetHashedAtomPairFingerprint`
        * Morgan: :func:`~rdkit.Chem.rdMolDescriptors.GetMorganFingerprint`
        * Hashed Morgan: :func:`~rdkit.Chem.rdMolDescriptors.GetHashedMorganFingerprint`
        * TopologicalTorsion: :func:`~rdkit.Chem.rdMolDescriptors.GetTopologicalTorsionFingerprint`
        * Hashed TopologicalTorsion: :func:`~rdkit.Chem.rdMolDescriptors.GetHashedTopologicalTorsionFingerprint`
        * RDKit: :func:`~rdkit.Chem.rdmolops.UnfoldedRDKFingerprintCountBased`
        * Hashed RDKit: :func:`~rdkit.Chem.rdmolops.RDKFingerprint`
        * MACCSKeys: :func:`~rdkit.Chem.rdMolDescriptors.GetMACCSKeysFingerprint`

    To generate a Morgan fingerprint, don't forget to specify the radius::

        get_fingerprint(ag, 'Morgan', radius=2)

    ``dtype="array"`` is not recommended for non-hashed fingerprints
    (except MACCS) as it is likely to overflow the available memory and crash.

    """
    if not isinstance(ag, AtomGroup):
        raise ValueError("The first argument must be an AtomGroup")
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
        >>> desc = RDKitDescriptors(u.atoms, "MolWt", "RadiusOfGyration",
        ...                         "MolFormula", num_atoms)
        >>> desc.run()
        >>> desc.results
        array([[46.06900000000002, 1.161278342193013, 'C2H6O', 9],
               [46.06900000000002, 1.175492972121405, 'C2H6O', 9],
               [46.06900000000002, 1.173230936577319, 'C2H6O', 9]],
              dtype=object)

    Notes
    -----
    Links to the modules from which descriptors are taken:

        * :mod:`rdkit.Chem.Descriptors`: Molecular descriptors
        * :mod:`rdkit.Chem.Descriptors3D`: Descriptors derived from a 
          molecule's 3D structure
        * :mod:`rdkit.Chem.EState.EState`: Basic EState definitions
        * :mod:`rdkit.Chem.EState.EState_VSA`: Hybrid EState-VSA descriptors
        * :mod:`rdkit.Chem.GraphDescriptors`: Topological/topochemical
          descriptors
        * :mod:`rdkit.Chem.Lipinski`: Lipinski parameters for molecules
        * :mod:`rdkit.Chem.MolSurf`: Approximate molecular surface area
          descriptors
        * :mod:`rdkit.Chem.rdMolDescriptors`: Molecular descriptors 
          (redundancies with :mod:`rdkit.Chem.Descriptors`)
    
    To get a list of all available descriptors, see :meth:`list_available`

    """

    def __init__(self, atomgroup, *descriptors, **kwargs):
        super().__init__(atomgroup.universe.trajectory,
                         **kwargs)
        if not isinstance(atomgroup, AtomGroup):
            raise ValueError("The first argument must be an AtomGroup")
        self.ag = atomgroup
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
        mol = self.ag.convert_to("RDKIT")
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
