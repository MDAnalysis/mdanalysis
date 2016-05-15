# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import print_function

import MDAnalysis
import MDAnalysis.analysis.psa

from numpy.testing import (TestCase, dec, assert_array_less,
                           assert_array_almost_equal, assert_)
import numpy as np

from MDAnalysisTests.datafiles import PSF, DCD, DCD2
from MDAnalysisTests import parser_not_found, tempdir


class TestPSAnalysis(TestCase):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.iu1 = np.triu_indices(3, k=1)
        self.universe1 = MDAnalysis.Universe(PSF, DCD)
        self.universe2 = MDAnalysis.Universe(PSF, DCD2)
        self.universe_rev = MDAnalysis.Universe(PSF, DCD)
        self.universes = [self.universe1, self.universe2, self.universe_rev]
        self.psa = MDAnalysis.analysis.psa.PSAnalysis(self.universes,           \
                                               path_select='name CA',           \
                                                      targetdir=self.tmpdir.name)
        self.psa.generate_paths(align=True)
        self.psa.paths[-1] = self.psa.paths[-1][::-1,:,:] # reverse third path
        self._run()

    def _run(self):
        self.psa.run(metric='hausdorff')
        self.hausd_matrix = self.psa.get_pairwise_distances()
        self.psa.run(metric='discrete_frechet')
        self.frech_matrix = self.psa.get_pairwise_distances()
        self.hausd_dists = self.hausd_matrix[self.iu1]
        self.frech_dists = self.frech_matrix[self.iu1]

    def tearDown(self):
        del self.universe1
        del self.universe2
        del self.universe_rev
        del self.psa
        del self.tmpdir

    def test_hausdorff_bound(self):
        err_msg = "Some Frechet distances are smaller than corresponding "      \
                + "Hausdorff distances"
        assert_array_less(self.hausd_dists, self.frech_dists, err_msg)

    def test_reversal_hausdorff(self):
        err_msg = "Hausdorff distances changed after path reversal"
        assert_array_almost_equal(self.hausd_matrix[1,2],                       \
                                  self.hausd_matrix[0,1],                       \
                                  decimal=3, err_msg=err_msg)

    def test_reversal_frechet(self):
        err_msg = "Frechet distances did not increase after path reversal"
        assert_(self.frech_matrix[1,2] >= self.frech_matrix[0,1], err_msg)
