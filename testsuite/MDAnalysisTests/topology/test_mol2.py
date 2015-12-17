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
from numpy.testing import (
    assert_,
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
