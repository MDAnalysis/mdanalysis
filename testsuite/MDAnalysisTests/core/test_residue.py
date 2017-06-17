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
    dec,
    assert_,
    assert_equal,
)

import MDAnalysis as mda

from MDAnalysisTests import parser_not_found
from MDAnalysisTests.datafiles import PSF, DCD


class TestResidue(object):
    # Legacy tests from before 363
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        self.res = self.universe.residues[100]

    def test_type(self):
        assert_(isinstance(self.res, mda.core.groups.Residue))
        assert_equal(self.res.resname, "ILE")
        assert_equal(self.res.resid, 101)

    def test_index(self):
        atom = self.res.atoms[2]
        assert_(isinstance(atom, mda.core.groups.Atom))
        assert_equal(atom.name, "CA")
        assert_equal(atom.index, 1522)
        assert_equal(atom.resid, 101)

    def test_atom_order(self):
        assert_equal(self.res.atoms.indices,
                     sorted(self.res.atoms.indices))


