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
)

import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    PQR,
    PQR_icodes,
)


class TestPQRParser(ParserBase):
    parser = mda.topology.PQRParser.PQRParser
    filename = PQR
    expected_attrs = ['ids', 'names', 'charges', 'radii',
                      'resids', 'resnames', 'icodes',
                      'segids']
    guessed_attrs = ['masses', 'types']
    expected_n_atoms = 3341
    expected_n_residues = 214
    expected_n_segments = 1

    def test_attr_size(self):
        assert_(len(self.top.ids) == self.top.n_atoms)
        assert_(len(self.top.names) == self.top.n_atoms)
        assert_(len(self.top.charges) == self.top.n_atoms)
        assert_(len(self.top.radii) == self.top.n_atoms)
        assert_(len(self.top.resids) == self.top.n_residues)
        assert_(len(self.top.resnames) == self.top.n_residues)
        assert_(len(self.top.segids) == self.top.n_segments)

class TestPQRParser2(TestPQRParser):
    filename = PQR_icodes

    expected_n_atoms = 5313
    expected_n_residues = 474
