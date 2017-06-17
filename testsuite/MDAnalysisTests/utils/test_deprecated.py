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

# test deprecated code
# - in particular stubs introduced in 0.11.0 (and which
#   will be removed in 1.0)
from __future__ import absolute_import

from numpy.testing import assert_raises

class TestImports(object):
    def test_core_units(self):
        try:
            import MDAnalysis.core.units
        except ImportError:
            raise AssertionError("MDAnalysis.core.units not available")

    def test_core_util(self):
        try:
            import MDAnalysis.core.util
        except ImportError:
            raise AssertionError("MDAnalysis.core.util not available")

    def test_core_log(self):
        try:
            import MDAnalysis.core.log
        except ImportError:
            raise AssertionError("MDAnalysis.core.log not available")

    def test_core_distances(self):
        try:
            import MDAnalysis.core.distances
        except ImportError:
            raise AssertionError("MDAnalysis.core.distances not available")

    def test_core_transformations(self):
        try:
            import MDAnalysis.core.transformations
        except ImportError:
            raise AssertionError("MDAnalysis.core.transformations not available")

    def test_core_qcprot(self):
        try:
            import MDAnalysis.core.qcprot
        except ImportError:
            raise AssertionError("MDAnalysis.core.qcprot not available")

    def test_KDTree(self):
        try:
            import MDAnalysis.KDTree
        except ImportError:
            raise AssertionError("MDAnalysis.KDTree not available")

    def test_analysis_x3dna(self):
        try:
            import MDAnalysis.analysis.x3dna
            from MDAnalysis.analysis.x3dna import X3DNA
        except ImportError:
            raise AssertionError("MDAnalysis.analysis.x3dna not available")

def test_collections_NotImplementedError():
    import MDAnalysis
    with assert_raises(NotImplementedError):
        MDAnalysis.collection.clear()





