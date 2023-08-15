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
from numpy.testing import assert_equal, assert_almost_equal

from MDAnalysis import units


class TestUnitEncoding(object):
    def test_unicode(self):
        try:
            assert_equal(units.lengthUnit_factor[u"\u212b"], 1.0)
        except KeyError:
            raise AssertionError("Unicode symbol for Angtrom not supported")

    def test_unicode_encoding_with_symbol(self):
        try:
            assert_equal(units.lengthUnit_factor[u"Å"], 1.0)
        except KeyError:
            raise AssertionError("UTF-8-encoded symbol for Angtrom not supported")


class TestConstants(object):
    # CODATA 2010 (NIST): http://physics.nist.gov/cuu/Constants/
    # (accessed 2015-02-15)
    # Add a reference value to this dict for every entry in
    # units.constants
    constants_reference = (
        ('N_Avogadro', 6.02214129e+23),  # mol**-1
        ('elementary_charge', 1.602176565e-19),  # As
        ('calorie', 4.184),  # J
        ('Boltzmann_constant', 8.314462159e-3),  # KJ (mol K)**-1
        ('Boltzman_constant', 8.314462159e-3),  # remove in 2.8.0
        ('electric_constant', 5.526350e-3),  # As (Angstroms Volts)**-1
    )

    @pytest.mark.parametrize('name, value', constants_reference)
    def test_constant(self, name, value):
        assert_almost_equal(units.constants[name], value)

    def test_boltzmann_typo_deprecation(self):
        wmsg = ("Please use 'Boltzmann_constant' henceforth. The key "
                "'Boltzman_constant' was a typo and will be removed "
                "in MDAnalysis 2.8.0.")
        with pytest.warns(DeprecationWarning, match=wmsg):
            units.constants['Boltzman_constant']


class TestConversion(object):
    @staticmethod
    def _assert_almost_equal_convert(value, u1, u2, ref):
        val = units.convert(value, u1, u2)
        assert_almost_equal(val, ref,
                            err_msg="Conversion {0} --> {1} failed".format(u1, u2))

    nm = 12.34567
    A = nm * 10.
    @pytest.mark.parametrize('quantity, unit1, unit2, ref', (
        (nm, 'nm', 'A', A),
        (A, 'Angstrom', 'nm', nm),
    ))
    def test_length(self, quantity, unit1, unit2, ref):
        self._assert_almost_equal_convert(quantity, unit1, unit2, ref)

    @pytest.mark.parametrize('quantity, unit1, unit2, ref', (
        (1, 'ps', 'AKMA', 20.45482949774598),
        (1, 'AKMA', 'ps', 0.04888821),
        (1, 'ps', 'ms', 1e-9),
        (1, 'ms', 'ps', 1e9),
        (1, 'ps', 'us', 1e-6),
        (1, 'us', 'ps', 1e6),
    ))
    def test_time(self, quantity, unit1, unit2, ref):
        self._assert_almost_equal_convert(quantity, unit1, unit2, ref)

    @pytest.mark.parametrize('quantity, unit1, unit2, ref', (
        (1, 'kcal/mol', 'kJ/mol', 4.184),
        (1, 'kcal/mol', 'eV', 0.0433641),
    ))
    def test_energy(self, quantity, unit1, unit2, ref):
        self._assert_almost_equal_convert(quantity, unit1, unit2, ref)

    @pytest.mark.parametrize('quantity, unit1, unit2, ref', (
        (1, 'kJ/(mol*A)', 'J/m', 1.66053892103219e-11),
        (2.5, 'kJ/(mol*nm)', 'kJ/(mol*A)', 0.25),
        (1, 'kcal/(mol*Angstrom)', 'kJ/(mol*Angstrom)', 4.184),
    ))
    def test_force(self, quantity, unit1, unit2, ref):
        self._assert_almost_equal_convert(quantity, unit1, unit2, ref)

    @pytest.mark.parametrize('quantity, unit1, unit2, ref', (
        (1, 'A/ps', 'm/s', 1e-10/1e-12),
        (1, 'A/ps', 'nm/ps', 0.1),
        (1, 'A/ps', 'pm/ps', 1e2),
        (1, 'A/ms', 'A/ps', 1e9),
        (1, 'A/us', 'A/ps', 1e6),
        (1, 'A/fs', 'A/ps', 1e-3),
        (1, 'A/AKMA', 'A/ps', 1/4.888821e-2),
    ))
    def test_speed(self, quantity, unit1, unit2, ref):
        self._assert_almost_equal_convert(quantity, unit1, unit2, ref)

    @pytest.mark.parametrize('quantity, unit1, unit2', ((nm, 'Stone', 'nm'),
                                                        (nm, 'nm', 'Stone')))
    def test_unit_unknown(self, quantity, unit1, unit2):
        with pytest.raises(ValueError):
            units.convert(quantity, unit1, unit2)

    def test_unit_unconvertable(self):
        nm = 12.34567
        A = nm * 10.
        with pytest.raises(ValueError):
            units.convert(A, 'A', 'ps')


class TestBaseUnits:
    @staticmethod
    @pytest.fixture
    def ref():
        # This is a copy of the dictionary we expect.
        # We want to know if base units are added or altered.
        ref = {"length": "A",
               "time": "ps",
               "energy": "kJ/mol",
               "charge": "e",
               "force": "kJ/(mol*A)",
               "speed": "A/ps"}
        return ref

    def test_MDANALYSIS_BASE_UNITS_correct(self, ref):
        assert ref == units.MDANALYSIS_BASE_UNITS
