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
import sys
from numpy.testing import assert_, assert_raises

import MDAnalysis as mda
from MDAnalysis.lib import lazy
from MDAnalysisTests import block_import


def _check_all_present(modnames):
    for modname in modnames: 
        assert_(modname in sys.modules)

class TestLazyExisting(object):
    modnames = ('MDAnalysis', 'MDAnalysis.analysis',
                'MDAnalysis.analysis.distances')

    # We attempt to run module functions (without arguments, which triggers
    #  TypeError exceptions) to see whether we can reach them.
    def test_load_base(self):
        MDAnalysis = lazy.import_module("MDAnalysis.analysis.distances",
                                        level='base')
        _check_all_present(self.modnames)
        assert_raises(TypeError, MDAnalysis.analysis.distances.dist)

    def test_load_leaf(self):
        distances = lazy.import_module("MDAnalysis.analysis.distances")
        _check_all_present(self.modnames)
        assert_raises(TypeError, distances.dist)

    def test_load_function(self):
        dist = lazy.import_function("MDAnalysis.analysis.distances.dist")
        _check_all_present(self.modnames)
        assert_raises(TypeError, dist)

    def test_load_functions(self):
        dist, dist_nonexistent = lazy.import_function("MDAnalysis.analysis.distances",
                                    "dist", "dist_nonexistent")
        _check_all_present(self.modnames)
        assert_raises(TypeError, dist)
        assert_raises(AttributeError, dist_nonexistent)


class TestLazyMissing(object):
    modnames = ('scipy', 'scipy.stats')

    # In this case failure occurs on accession, so we must test for that,
    #  rather than function behavior.
    @block_import('scipy')
    def test_load_base(self):
        scipy = lazy.import_module("scipy.stats", level='base')
        _check_all_present(self.modnames)
        assert_raises(ImportError, getattr, scipy, 'stats')

    @block_import('scipy')
    def test_load_leaf(self):
        stats = lazy.import_module("scipy.stats")
        _check_all_present(self.modnames)
        assert_raises(ImportError, getattr, stats, 'anderson')

    @block_import('scipy')
    def test_load_function(self):
        func1 = lazy.import_function("scipy.stats.anderson")
        _check_all_present(self.modnames)
        assert_raises(ImportError, func1)

    @block_import('scipy')
    def test_load_functions(self):
        func1, func2 = lazy.import_function("scipy.stats",
                                            "anderson", "whatever_")
        _check_all_present(self.modnames)
        assert_raises(ImportError, func1)
        assert_raises(ImportError, func2)

