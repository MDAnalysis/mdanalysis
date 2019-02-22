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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import
import MDAnalysis as mda
import pytest

from MDAnalysisTests.datafiles import PSF_TRICLINIC, DCD_TRICLINIC
from MDAnalysis.analysis.dielectric import DielectricConstant
from numpy.testing import assert_almost_equal
    
class TestDielectric(object):
    @pytest.fixture()
    def ag(self):
        return mda.Universe(PSF_TRICLINIC, DCD_TRICLINIC).atoms
        
    def test_water(self, ag):
        eps = DielectricConstant(ag).run()
        assert_almost_equal(eps.results['eps_mean'], 3.874, decimal=1)
        
    def test_make_whole(self, ag):
        # cut molecules apart
        ag.universe.transfer_to_memory()
        for ts in ag.universe.trajectory:
            ag.wrap()
    
        eps_broken = DielectricConstant(ag, make_whole=False).run()
        assert_almost_equal(eps_broken.results['eps_mean'], 721.711, decimal=1)
        
        eps_whole = DielectricConstant(ag, make_whole=True).run()
        assert_almost_equal(eps_whole.results['eps_mean'], 3.874, decimal=1)
        
    def test_temperature(self, ag):        
        eps = DielectricConstant(ag, temperature=100).run()
        assert_almost_equal(eps.results['eps_mean'], 9.621, decimal=1)
        
    def test_non_neutral(self, ag):        
        with pytest.raises(NotImplementedError):
            DielectricConstant(ag[:-1]).run()
            
    def test_free_charges(self, ag):
        ag.fragments[0].charges += 1
        ag.fragments[1].charges -= 1
        
        with pytest.raises(NotImplementedError):
            DielectricConstant(ag).run()
