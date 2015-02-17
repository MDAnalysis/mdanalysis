# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

import numpy as np
from numpy.testing import *

import MDAnalysis.core.units as units
from MDAnalysis.core import flags


class TestDefaultUnits(TestCase):
    def testLength(self):
        assert_equal(flags['length_unit'], 'Angstrom',
                     "The default length unit should be Angstrom (in core.flags)")

    def testTime(self):
        assert_equal(flags['time_unit'], 'ps',
                     "The default length unit should be pico seconds (in core.flags)")

    def testConvertGromacsTrajectories(self):
        assert_equal(flags['convert_lengths'], True,
                     "The default behaviour should be to auto-convert Gromacs trajectories")

class TestConstants(object):
    # CODATA 2010 (NIST): http://physics.nist.gov/cuu/Constants/
    # (accessed 2015-02-15)
    # Add a reference value to this dict for every entry in
    # units.constants
    constants_reference = {
        'N_Avogadro': 6.02214129e+23,         # mol**-1
        'elementary_charge': 1.602176565e-19, # As
        }

    def test_constant(self):
        for name,value in self.constants_reference.iteritems():
            yield self.check_physical_constant, name, value

    @staticmethod
    def check_physical_constant(name, reference):
        assert_almost_equal(units.constants[name], reference)


class TestConversion(TestCase):
    def testLength(self):
        nm = 12.34567
        A = nm * 10.
        assert_almost_equal(units.convert(nm, 'nm', 'A'), A,
                            err_msg="Conversion nm -> A failed")
        assert_almost_equal(units.convert(A, 'Angstrom', 'nm'), nm,
                            err_msg="Conversion A -> nm failed")

    def testTime(self):
        assert_almost_equal(units.convert(1, 'ps', 'AKMA'), 20.45482949774598,
                            err_msg="Conversion ps -> AKMA failed")
        assert_almost_equal(units.convert(1, 'AKMA', 'ps'), 0.04888821,
                            err_msg="Conversion AKMA -> ps failed")

    def testEnergy(self):
        assert_almost_equal(units.convert(1, 'kcal/mol', 'kJ/mol'), 4.184,
                            err_msg="Conversion kcal/mol -> kJ/mol failed")
        assert_almost_equal(units.convert(1, 'kcal/mol', 'eV'), 0.0433641,
                            err_msg="Conversion kcal/mol -> eV failed")
