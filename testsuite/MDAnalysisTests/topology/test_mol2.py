# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import
from numpy.testing import (
    assert_,
    assert_equal,
)

import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    mol2_molecule,
    mol2_molecules,
)


class MOL2Base(ParserBase):
    parser = mda.topology.MOL2Parser.MOL2Parser
    expected_attrs = ['ids', 'names', 'types', 'charges',
                      'resids', 'resnames', 'bonds']
    guessed_attrs = ['elements', 'masses']
    expected_n_atoms = 49
    expected_n_residues = 1
    expected_n_segments = 1

    def test_attr_size(self):
        assert_(len(self.top.ids) == self.top.n_atoms)
        assert_(len(self.top.names) == self.top.n_atoms)
        assert_(len(self.top.types) == self.top.n_atoms)
        assert_(len(self.top.charges) == self.top.n_atoms)
        assert_(len(self.top.resids) == self.top.n_residues)
        assert_(len(self.top.resnames) == self.top.n_residues)

    def test_bonds(self):
        assert_(len(self.top.bonds) == 49)  # bonds for 49 atoms
        assert_(len(self.top.bonds.values) == 51)  # this many bonds

    
class TestMOL2Parser(MOL2Base):
    filename = mol2_molecule


class TestMOL2Parser2(MOL2Base):
    filename = mol2_molecules


def test_bond_orders():
    ref_orders = ('am 1 1 2 1 2 1 1 am 1 1 am 2 2 '
                  '1 1 1 1 1 1 1 1 1 1 1 1 1 1 '
                  'ar ar ar 1 ar 1 ar 1 ar 1 1 1 '
                  '2 1 1 1 1 2 1 1 2 1 1').split()
    u = mda.Universe(mol2_molecule)
    orders = [bond.order for bond in u.atoms.bonds]
    assert_equal(orders, ref_orders)
