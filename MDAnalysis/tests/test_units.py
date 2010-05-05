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
        assert_equal(flags['convert_gromacs_lengths'], True,
                     "The default behaviour should be to auto-convert Gromacs trajectories")

class TestConversion(TestCase):
    def testLength(TestCase):
        nm = 12.34567
        A = nm * 10.
        assert_almost_equal(units.convert(nm, 'nm', 'A'), A,
                            err_msg="Conversion nm -> A failed")
        assert_almost_equal(units.convert(A, 'Angstrom', 'nm'), nm,
                            err_msg="Conversion A -> nm failed")

    def testTime(TestCase):
        assert_almost_equal(units.convert(1, 'ps', 'AKMA'), 20.45482949774598,
                            err_msg="Conversion ps -> AKMA failed")
        assert_almost_equal(units.convert(1, 'AKMA', 'ps'), 0.04888821,
                            err_msg="Conversion AKMA -> ps failed")
