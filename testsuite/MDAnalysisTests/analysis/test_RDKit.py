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
class TestFingerprintsRDKit:
    @pytest.fixture
    def u(self):
        return mda.Universe.from_smiles("CCO")

    def test_hashed_maccs(self, u):
        with pytest.raises(ValueError,
                           match="MACCSKeys is not available in a hashed version"):
            get_fingerprint(u.atoms, "MACCSKeys", hashed=True)
        fp = get_fingerprint(u.atoms, "MACCSKeys", hashed=False, 
                             as_array=False)
        assert list(fp.GetOnBits()) == [82, 109, 112, 114, 126, 138, 139, 153, 155, 157, 160, 164]

    def test_unknown_fp(self, u):
        with pytest.raises(ValueError,
                           match="Could not find 'foo' in the available fingerprints"):
            get_fingerprint(u.atoms, "foo", hashed=False)

    @pytest.mark.parametrize("kind, hashed, as_array, dtype", [
        ("MACCSKeys", False, False, "ExplicitBitVect"),
        ("AtomPair", True, False, "IntSparseIntVect"),
        ("AtomPair", True, True, "ndarray"),
    ])
    def test_return_types(self, u, kind, hashed, as_array, dtype):
        fp = get_fingerprint(u.atoms, kind, hashed, as_array)
        assert fp.__class__.__name__ == dtype

    def test_kwargs(self, u):
        fp = get_fingerprint(u.atoms, "AtomPair", hashed=True, as_array=True,
                             nBits=128)
        assert len(fp) == 128

    @pytest.mark.parametrize("kind, kwargs, hashed, as_array, n_on_bits", [
        ("MACCSKeys", {}, False, False, 12),
        ("AtomPair", {}, True, True, 5),
        ("AtomPair", {}, False, False, 12),
        ("Morgan", dict(radius=2), True, True, 8),
        ("Morgan", dict(radius=2), False, False, 11),
        ("RDKit", {}, True, True, 84),
        ("RDKit", {}, False, False, 42),
        ("TopologicalTorsion", {}, True, True, 0),
        ("TopologicalTorsion", {}, False, False, 0),
    ])
    def test_fp(self, u, kind, kwargs, hashed, as_array, n_on_bits):
        fp = get_fingerprint(u.atoms, kind,
                             hashed=hashed, as_array=as_array, **kwargs)
        classname = fp.__class__.__name__
        if classname.endswith("BitVect"):
            assert fp.GetNumOnBits() == n_on_bits
        elif classname.endswith("IntVect"):
            assert len(fp.GetNonzeroElements()) == n_on_bits
        else:
            assert len(np.where(fp == 1)[0]) == n_on_bits
