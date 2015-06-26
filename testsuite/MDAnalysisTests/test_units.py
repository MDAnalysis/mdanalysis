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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import unicode_literals

import numpy as np
from numpy.testing import *

from MDAnalysis import units
from MDAnalysis.core import flags


class TestDefaultUnits(TestCase):
    def testLength(self):
        assert_equal(flags['length_unit'], 'Angstrom',
                     b"The default length unit should be Angstrom (in core.flags)")

    def testTime(self):
        assert_equal(flags['time_unit'], 'ps',
                     b"The default length unit should be pico seconds (in core.flags)")

    def testConvertGromacsTrajectories(self):
        assert_equal(flags['convert_lengths'], True,
                     b"The default behaviour should be to auto-convert Gromacs trajectories")


class TestUnitEncoding(TestCase):
    def testUnicode(self):
        try:
            assert_equal(units.lengthUnit_factor["\u212b"], 1.0)
        except KeyError:
            raise AssertionError("Unicode symbol for Angtrom not supported")

    def testUTF8Encoding(self):
        try:
            assert_equal(units.lengthUnit_factor[b"\xe2\x84\xab"], 1.0)
        except KeyError:
            raise AssertionError("UTF-8-encoded symbol for Angtrom not supported")


class TestConstants(object):
    # CODATA 2010 (NIST): http://physics.nist.gov/cuu/Constants/
    # (accessed 2015-02-15)
    # Add a reference value to this dict for every entry in
    # units.constants
    constants_reference = {
        'N_Avogadro': 6.02214129e+23,          # mol**-1
        'elementary_charge': 1.602176565e-19,  # As
        'calorie': 4.184,                      # J
        }

    def test_constant(self):
        for name, value in self.constants_reference.iteritems():
            yield self.check_physical_constant, name, value

    @staticmethod
    def check_physical_constant(name, reference):
        assert_almost_equal(units.constants[name], reference)


class TestConversion(object):
    @staticmethod
    def _assert_almost_equal_convert(value, u1, u2, ref):
        assert_almost_equal(units.convert(value, u1, u2), ref,
                            err_msg="Conversion {0} --> {1} failed".format(u1, u2))

    # generate individual test cases using nose's test generator mechanism
    def testLength(self):
        nm = 12.34567
        A = nm * 10.
        yield self._assert_almost_equal_convert, nm, 'nm', 'A', A
        yield self._assert_almost_equal_convert, A, 'Angstrom', 'nm', nm

    def testTime(self):
        yield self._assert_almost_equal_convert, 1, 'ps', 'AKMA', 20.45482949774598
        yield self._assert_almost_equal_convert, 1, 'AKMA', 'ps', 0.04888821

    def testEnergy(self):
        yield self._assert_almost_equal_convert, 1, 'kcal/mol', 'kJ/mol', 4.184
        yield self._assert_almost_equal_convert, 1, 'kcal/mol', 'eV', 0.0433641

    def testForce(self):
        yield self._assert_almost_equal_convert, 1, 'kJ/(mol*A)', 'J/m', 1.66053892103219e-11
        yield self._assert_almost_equal_convert, 2.5, 'kJ/(mol*nm)', 'kJ/(mol*A)', 0.25
        yield self._assert_almost_equal_convert, 1, 'kcal/(mol*Angstrom)', 'kJ/(mol*Angstrom)', 4.184


