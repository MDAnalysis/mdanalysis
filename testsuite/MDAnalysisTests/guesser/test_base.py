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
import numpy as np
import MDAnalysis as mda
from MDAnalysis.guesser.base import GuesserBase, get_guesser
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import Masses, Atomnames
import MDAnalysis.tests.datafiles as datafiles
from numpy.testing import assert_allclose

class TesttBaseGuesser():
    
    def test_get_guesser(self):
        class TestGuesser1(GuesserBase):
            context = 'test1'

        class TestGuesser2(TestGuesser1):
            context = 'test2'

        assert get_guesser(TestGuesser1()).context == 'test1'
        assert get_guesser(TestGuesser2()).context == 'test2'

    def test_guess_invalid_attribute(self):
        with pytest.raises(ValueError, match='default guesser can not guess one or more of the provided attributes'):
            mda.Universe(datafiles.PDB, to_guess=['foo'])

    def test_guess_attribute_with_missing_parent_attr(self):
        names = Atomnames(np.array(['C', 'HB', 'HA', 'O'], dtype=object))
        masses = Masses(np.array([np.nan, np.nan, np.nan, np.nan], dtype=np.float64))
        top = Topology(4, 1, 1, attrs=[names,masses,])
        u = mda.Universe(top, to_guess=['masses'])
        assert hasattr(u.atoms, 'elements')
        assert_allclose(u.atoms.masses, np.array([12.01100, 1.00800, 1.00800, 15.99900]), rtol=1e-3, atol=0)
        
    def test_failed_to_find_parent_attr(self):
        masses = Masses(np.array([np.nan, np.nan, np.nan, np.nan], dtype=np.float64))
        top = Topology(4, 1, 1, attrs=[masses,])
        with pytest.raises(ValueError, match='Failed to find/guess a parent attribute to guess masses from'):
            mda.Universe(top, to_guess=['masses'])
        
          
        
        
        