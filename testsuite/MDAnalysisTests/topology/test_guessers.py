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
    assert_array_equal,
    assert_warns,
)
import numpy as np

from MDAnalysis.topology import guessers

class TestGuessMasses(object):
    @staticmethod
    def test_guess_masses():
        out = guessers.guess_masses(['C', 'C', 'H'])
        
        assert_(isinstance(out, np.ndarray))
        assert_array_equal(out,
                           np.array([12.011, 12.011, 1.008]))

    @staticmethod
    def test_guess_masses_warn():
        assert_warns(UserWarning, guessers.guess_masses, ['X'])

    @staticmethod
    def test_guess_masses_miss():
        out = guessers.guess_masses(['X', 'Z'])
        assert_array_equal(out, np.array([0.0, 0.0]))

    @staticmethod
    def test_get_atom_mass():
        assert_(guessers.get_atom_mass('H') == 1.008)

    @staticmethod
    def test_get_atom_mass_miss():
        assert_(guessers.get_atom_mass('XYZ') == 0.0)

    @staticmethod
    def test_guess_atom_mass():
        assert_(guessers.guess_atom_mass('1H') == 1.008)


class TestGuessTypes(object):
    # guess_types
    # guess_atom_type
    # guess_atom_element
    @staticmethod
    def test_guess_types():
        out = guessers.guess_types(['MG2+', 'C12'])

        assert_(isinstance(out, np.ndarray))
        assert_array_equal(out, np.array(['MG', 'C'], dtype=object))
    
    @staticmethod
    def test_guess_atom_element():
        assert_(guessers.guess_atom_element('MG2+') == 'MG')
    
    @staticmethod
    def test_guess_atom_element_empty():
        assert_(guessers.guess_atom_element('') == '')

    @staticmethod
    def test_guess_atom_element_singledigit():
        assert_(guessers.guess_atom_element('1') == '1')

    @staticmethod
    def test_guess_atom_element_1H():
        assert_(guessers.guess_atom_element('1H') == 'H')
        assert_(guessers.guess_atom_element('2H') == 'H')

def test_guess_charge():
    # this always returns 0.0
    assert_(guessers.guess_atom_charge('this') == 0.0)

