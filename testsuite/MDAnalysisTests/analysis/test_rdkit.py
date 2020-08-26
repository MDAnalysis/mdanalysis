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

import numpy as np
import MDAnalysis as mda
from MDAnalysisTests.util import import_not_available
try:
    from MDAnalysis.analysis.RDKit import RDKitDescriptors, get_fingerprint
    from rdkit.Chem import AllChem
except ImportError:
    pass

from numpy.testing import assert_equal

import pytest


requires_rdkit = pytest.mark.skipif(import_not_available("rdkit"),
                                    reason="requires RDKit")


@pytest.mark.skipif(not import_not_available("rdkit"),
                    reason="only for min dependencies build")
class TestRequiresRDKit(object):
    def test_requires_rdkit(self):
        with pytest.raises(ImportError, match="RDKit is needed"):
            from MDAnalysis.analysis.RDKit import RDKitDescriptors


@requires_rdkit
class TestDescriptorsRDKit:
    @pytest.fixture
    def u(self):
        return mda.Universe.from_smiles("CCO")

    def test_not_ag(self, u):
        with pytest.raises(ValueError,
                           match="must be an AtomGroup"):
            RDKitDescriptors(u, "MolWt").run()

    def test_arg_function(self, u):
        def func(mol):
            return mol.GetNumAtoms()
        desc = RDKitDescriptors(u.atoms, func).run()
        assert desc.results[0, 0] == 9

    def test_arg_string(self, u):
        desc = RDKitDescriptors(u.atoms, "MolWt", "MolFormula").run()
        assert round(desc.results[0, 0], 10) == 46.069
        # MolFormula is known as "CalcMolFormula" in RDKit but should still
        # match even with the missing Calc
        assert desc.results[0, 1] == "C2H6O"

    def test_multi_frame(self):
        u = mda.Universe.from_smiles("CCO", numConfs=3)
        desc = RDKitDescriptors(u.atoms, "RadiusOfGyration").run()
        assert desc.results.shape == (3, 1)

    def test_unknown_desc(self, u):
        with pytest.raises(ValueError, match="Could not find 'foo' in RDKit"):
            desc = RDKitDescriptors(u.atoms, "foo")

    def test_list_available(self):
        avail = RDKitDescriptors.list_available()
        assert "rdkit.Chem.Descriptors" in avail.keys()
        assert "MolWt" in avail["rdkit.Chem.Descriptors"]
        flat = RDKitDescriptors.list_available(flat=True)
        assert "CalcAUTOCORR2D" in flat


@requires_rdkit
class FingerprintFactory:
    def GetAtomPairFingerprint(*args, **kwargs):
        return AllChem.GetAtomPairFingerprint(*args, **kwargs)

    def GetHashedAtomPairFingerprint(*args, **kwargs):
        return AllChem.GetHashedAtomPairFingerprint(*args, **kwargs)

    def GetMorganFingerprint(*args, **kwargs):
        return AllChem.GetMorganFingerprint(*args, **kwargs)

    def GetHashedMorganFingerprint(*args, **kwargs):
        return AllChem.GetHashedMorganFingerprint(*args, **kwargs)

    def GetTopologicalTorsionFingerprint(*args, **kwargs):
        return AllChem.GetTopologicalTorsionFingerprint(*args, **kwargs)

    def GetHashedTopologicalTorsionFingerprint(*args, **kwargs):
        return AllChem.GetHashedTopologicalTorsionFingerprint(*args, **kwargs)

    def UnfoldedRDKFingerprintCountBased(*args, **kwargs):
        return AllChem.UnfoldedRDKFingerprintCountBased(*args, **kwargs)

    def RDKFingerprint(*args, **kwargs):
        return AllChem.RDKFingerprint(*args, **kwargs)

    def GetMACCSKeysFingerprint(*args, **kwargs):
        return AllChem.GetMACCSKeysFingerprint(*args, **kwargs)

    def null(*args, **kwargs):
        pass

@pytest.fixture(scope="function")
def fp_func(request):
    return getattr(FingerprintFactory, request.param)


@requires_rdkit
class TestFingerprintsRDKit:
    @pytest.fixture
    def u(self):
        return mda.Universe.from_smiles("CCO")

    def test_not_ag(self, u):
        with pytest.raises(ValueError,
                           match="must be an AtomGroup"):
            get_fingerprint(u, "MACCSKeys")

    def test_hashed_maccs(self, u):
        with pytest.raises(ValueError,
                           match="MACCSKeys is not available in a hashed version"):
            get_fingerprint(u.atoms, "MACCSKeys", hashed=True)
        fp = get_fingerprint(u.atoms, "MACCSKeys", hashed=False, dtype="dict")
        assert list(fp.keys()) == [82, 109, 112, 114, 126, 138, 139, 153, 155, 157, 160, 164]

    def test_unknown_fp(self, u):
        with pytest.raises(ValueError,
                           match="Could not find 'foo' in the available fingerprints"):
            get_fingerprint(u.atoms, "foo")

    @pytest.mark.parametrize("kind, hashed, dtype, out_type", [
        ("MACCSKeys", False, None, "ExplicitBitVect"),
        ("AtomPair", True, None, "IntSparseIntVect"),
        ("AtomPair", True, "array", "ndarray"),
        ("AtomPair", False, "dict", "dict"),
    ])
    def test_return_types(self, u, kind, hashed, dtype, out_type):
        fp = get_fingerprint(u.atoms, kind, hashed, dtype)
        assert fp.__class__.__name__ == out_type

    def test_unknown_dtype(self, u):
        with pytest.raises(ValueError,
                           match="'foo' is not a supported output type"):
            get_fingerprint(u.atoms, "MACCSKeys", dtype="foo")

    def test_kwargs(self, u):
        fp = get_fingerprint(u.atoms, "AtomPair", hashed=True, dtype="array",
                             nBits=128)
        assert len(fp) == 128

    @pytest.mark.parametrize("kind, kwargs, hashed, dtype, n_on_bits, fp_func", [
        ("MACCSKeys", {}, False, "array", 12, "null"),
        ("MACCSKeys", {}, False, "dict", 12, "null"),
        ("MACCSKeys", {}, False, None, 12, "GetMACCSKeysFingerprint"),
        ("AtomPair", {}, True, "array", 5, "null"),
        ("AtomPair", {}, False, "dict", 12, "null"),
        ("AtomPair", {}, False, None, 12, "GetAtomPairFingerprint"),
        ("AtomPair", {}, True, None, 5, "GetHashedAtomPairFingerprint"),
        ("Morgan", dict(radius=2), True, "array", 8, "null"),
        ("Morgan", dict(radius=2), False, "dict", 11, "null"),
        ("Morgan", dict(radius=2), False, None, 11, "GetMorganFingerprint"),
        ("Morgan", dict(radius=2), True, None, 8, "GetHashedMorganFingerprint"),
        ("RDKit", {}, True, "array", 84, "null"),
        ("RDKit", {}, False, "dict", 42, "null"),
        ("RDKit", {}, False, None, 42, "UnfoldedRDKFingerprintCountBased"),
        ("RDKit", {}, True, None, 84, "RDKFingerprint"),
        ("TopologicalTorsion", {}, True, "array", 0, "null"),
        ("TopologicalTorsion", {}, False, "dict", 0, "null"),
        ("TopologicalTorsion", {}, False, None, 0, "GetTopologicalTorsionFingerprint"),
        ("TopologicalTorsion", {}, True, None, 0, "GetHashedTopologicalTorsionFingerprint"),
    ], indirect=["fp_func"])
    def test_fp(self, u, kind, kwargs, hashed, dtype, n_on_bits, fp_func):
        fp = get_fingerprint(u.atoms, kind,
                             hashed=hashed, dtype=dtype, **kwargs)
        if fp_func != "null":
            mol = u.atoms.convert_to("RDKIT")
            rdkit_fp = fp_func(mol, **kwargs)
        classname = fp.__class__.__name__
        if classname.endswith("BitVect"):
            assert list(fp.GetOnBits()) == list(rdkit_fp.GetOnBits())
        elif classname.endswith("IntVect"):
            assert fp.GetNonzeroElements() == rdkit_fp.GetNonzeroElements()
        elif classname == "dict":
            assert len(fp.keys()) == n_on_bits
        else:
            assert len(np.where(fp == 1)[0]) == n_on_bits
