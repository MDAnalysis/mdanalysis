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
import pytest
import MDAnalysis as mda
from MDAnalysisTests.util import import_not_available


requires_rdkit = pytest.mark.skipif(import_not_available("rdkit"),
                                    reason="requires RDKit")


@requires_rdkit
class TestConvertTo:
    @pytest.fixture(scope="class")
    def u(self):
        return mda.Universe.from_smiles("CCO")

    def test_convert_to_case_insensitive(self, u):
        mol = u.atoms.convert_to("rdkit")

    def test_convert_to_lib_as_method(self, u):
        mol = u.atoms.convert_to.rdkit()

    def test_convert_to_kwargs(self, u):
        mol = u.atoms.convert_to("RDKIT", NoImplicit=False)
        assert mol.GetAtomWithIdx(0).GetNoImplicit() is False

    def test_convert_to_lib_method_kwargs(self, u):
        mol = u.atoms.convert_to.rdkit(NoImplicit=False)
        assert mol.GetAtomWithIdx(0).GetNoImplicit() is False


class TestAccessor:
    def test_access_from_class(self):
        assert (mda.core.AtomGroup.convert_to is
                mda.core.accessors.ConverterWrapper)


class TestConverterWrapper:
    def test_raises_valueerror(self):
        u = mda.Universe.empty(1)
        with pytest.raises(ValueError,
                           match="No 'mdanalysis' converter found"):
            u.atoms.convert_to("mdanalysis")

    @requires_rdkit
    def test_single_instance(self):
        u1 = mda.Universe.from_smiles("C")
        u2 = mda.Universe.from_smiles("CC")
        assert (u1.atoms.convert_to.rdkit.__wrapped__ is
                u2.atoms.convert_to.rdkit.__wrapped__)
