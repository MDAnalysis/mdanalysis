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

"""32 bit compat tests

Tests for making sure that integer arrays used for indexing use `np.intp`.
This dtype is important for platform independent indexing of other arrays.

"""
from __future__ import absolute_import

import numpy as np
from numpy.testing import (
    assert_,
)
from MDAnalysisTests import make_Universe


class TestIndexDtype(object):
    def setUp(self):
        self.u = make_Universe()

    def tearDown(self):
        del self.u

    def test_ag_ix(self):
        assert_(self.u.atoms.ix.dtype == np.intp)

    def test_rg_ix(self):
        assert_(self.u.residues.ix.dtype == np.intp)

    def test_sg_ix(self):
        assert_(self.u.segments.ix.dtype == np.intp)

    def test_atom_ix_array(self):
        assert_(self.u.atoms[0].ix_array.dtype == np.intp)

    def test_residue_ix_array(self):
        assert_(self.u.residues[0].ix_array.dtype == np.intp)

    def test_segment_ix_array(self):
        assert_(self.u.segments[0].ix_array.dtype == np.intp)

    def test_atomgroup_indices(self):
        assert_(self.u.atoms.indices.dtype == np.intp)

    def test_atomgroup_residue_upshift(self):
        assert_(self.u.atoms.resindices.dtype == np.intp)

    def test_atomgroup_segment_upshift(self):
        assert_(self.u.atoms.segindices.dtype == np.intp)

    def test_residuegroup_atom_downshift(self):
        # downshift arrays are a list (one for each residue)
        assert_(all((arr.dtype == np.intp)
                    for arr in self.u.residues.indices))

    def test_residuegroup_resindices(self):
        assert_(self.u.residues.resindices.dtype == np.intp)

    def test_residuegroup_segment_upshift(self):
        assert_(self.u.residues.segindices.dtype == np.intp)

    def test_segmentgroup_atom_downshift(self):
        assert_(all((arr.dtype == np.intp)
                    for arr in self.u.segments.indices))

    def test_segmentgroup_residue_downshift(self):
        assert_(all((arr.dtype == np.intp)
                    for arr in self.u.segments.resindices))

    def test_segmentgroup_segindices(self):
        assert_(self.u.segments.segindices.dtype == np.intp)
