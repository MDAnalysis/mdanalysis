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

import pytest
from numpy.testing import (
    assert_,
    assert_array_equal,
    assert_equal,
    assert_warns,
)
import numpy as np

from MDAnalysis.topology import guessers
from MDAnalysis.core.topologyattrs import Angles

from MDAnalysisTests import make_Universe
from MDAnalysisTests.core.test_fragments import make_starshape

class TestGuessMasses(object):
    def test_guess_masses(self):
        out = guessers.guess_masses(['C', 'C', 'H'])
        
        assert_(isinstance(out, np.ndarray))
        assert_array_equal(out,
                           np.array([12.011, 12.011, 1.008]))

    def test_guess_masses_warn(self):
        assert_warns(UserWarning, guessers.guess_masses, ['X'])

    def test_guess_masses_miss(self):
        out = guessers.guess_masses(['X', 'Z'])
        assert_array_equal(out, np.array([0.0, 0.0]))

    def test_get_atom_mass(self):
        assert_(guessers.get_atom_mass('H') == 1.008)

    def test_get_atom_mass_miss(self):
        assert_(guessers.get_atom_mass('XYZ') == 0.0)

    def test_guess_atom_mass(self):
        assert_(guessers.guess_atom_mass('1H') == 1.008)


class TestGuessTypes(object):
    # guess_types
    # guess_atom_type
    # guess_atom_element
    def test_guess_types(self):
        out = guessers.guess_types(['MG2+', 'C12'])

        assert_(isinstance(out, np.ndarray))
        assert_array_equal(out, np.array(['MG', 'C'], dtype=object))

    def test_guess_atom_element(self):
        assert_(guessers.guess_atom_element('MG2+') == 'MG')

    def test_guess_atom_element_empty(self):
        assert_(guessers.guess_atom_element('') == '')

    def test_guess_atom_element_singledigit(self):
        assert_(guessers.guess_atom_element('1') == '1')

    def test_guess_atom_element_1H(self):
        assert_(guessers.guess_atom_element('1H') == 'H')
        assert_(guessers.guess_atom_element('2H') == 'H')

def test_guess_charge():
    # this always returns 0.0
    assert_(guessers.guess_atom_charge('this') == 0.0)

def test_guess_bonds_Error():
    u = make_Universe(trajectory=True)
    with pytest.raises(ValueError):
        guessers.guess_bonds( u.atoms[:4],
                u.atoms.positions[:5])
    
def test_guess_impropers():
    u = make_starshape()

    ag = u.atoms[:5]

    u.add_TopologyAttr(Angles(guessers.guess_angles(ag.bonds)))

    vals = guessers.guess_improper_dihedrals(ag.angles)
    assert_equal(len(vals), 12)
