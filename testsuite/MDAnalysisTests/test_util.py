# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

import numpy
import numpy.random
from numpy.testing import *
from numpy import pi, sin, cos

import MDAnalysis.core.util as util

def check_parse_residue(rstring, residue):
    assert_equal(util.parse_residue(rstring), residue)

def check_convert_aa_3to1(resname3, resname1):
    assert_equal(util.convert_aa_code(resname3), resname1)

def check_convert_aa_1to3(resname1, resname3_canonical):
    assert_equal(util.convert_aa_code(resname1), resname3_canonical)

class TestStringFunctions(object):
    # (1-letter, (canonical 3 letter, other 3/4 letter, ....))
    aa = [('H', ('HIS', 'HISA', 'HISB', 'HSE', 'HSD', 'HIS1', 'HIS2', 'HIE', 'HID')),
          ('K', ('LYS', 'LYSH', 'LYN')),
          ('A', ('ALA',)),
          ('D', ('ASP', 'ASPH', 'ASH')),
          ('E', ('GLU', 'GLUH', 'GLH')),
          ('N', ('ASN',)),
          ('Q', ('GLN',)),
          ('C', ('CYS', 'CYSH', 'CYS1', 'CYS2')),
          ]

    residues = [("LYS300:HZ1", ("LYS", 300, "HZ1")),
                ("K300:HZ1", ("LYS", 300, "HZ1")),
                ("K300", ("LYS", 300, None)),
                ("LYS 300:HZ1", ("LYS", 300, "HZ1")),
                ("M1:CA", ("MET", 1, "CA")),
                ]

    def testParse_residue(self):
        for rstring, residue in self.residues:
            yield check_parse_residue, rstring, residue

    def test_convert_aa_code_long(self):
        for resname1, strings in self.aa:
            for resname3 in strings:
                yield check_convert_aa_3to1, resname3, resname1

    def test_convert_aa_code_short(self):
        for resname1, strings in self.aa:
            yield check_convert_aa_1to3, resname1, strings[0]


class TestIterable(TestCase):
    def test_lists(self):
        assert_equal(util.iterable([1, 2, 3]), True)
        assert_equal(util.iterable([]), True)
    def test_tuples(self):
        assert_equal(util.iterable((1, 2, 3)), True)
        assert_equal(util.iterable(()), True)
    def test_iterator(self):
        assert_equal(util.iterable(xrange(3)), True)
    def test_arrays(self):
        assert_equal(util.iterable(numpy.array([1,2,3])), True)
    def test_scalars(self):
        assert_equal(util.iterable(123), False)
    def test_strings(self):
        """Test that iterable() works on any string (Fixed bug)"""
        assert_equal(util.iterable("byte string"), False)
        assert_equal(util.iterable(u"unicode string"), False)

class TestGeometryFunctions(TestCase):
    def setUp(self):
        self.e1 = numpy.array([1.,0,0])
        self.e2 = numpy.array([0,1.,0])
        self.e3 = numpy.array([0,0,1.])
        self.a  = numpy.array([cos(pi/3),sin(pi/3),0])
        self.null = numpy.zeros(3)

    def testAngleUnitvectors(self):
        assert_equal(util.angle(self.e1, self.e2), pi/2)
        assert_equal(util.angle(self.e1, self.a), pi/3)

    def testAngleVectors(self):
        assert_equal(util.angle(2*self.e1, self.e2), pi/2)
        assert_equal(util.angle(-2*self.e1, self.e2), pi - pi/2)
        assert_equal(util.angle(23.3*self.e1, self.a), pi/3)

    def testAngleNullVector(self):
        assert_equal(util.angle(self.e1, self.null), numpy.nan)

    def testAngleColinear(self):
        assert_equal(util.angle(self.a, self.a), 0.0)

    def testAnglePi(self):
        assert_almost_equal(util.angle(-2.3456e7*self.e1, 3.4567e-6*self.e1), pi)
        assert_almost_equal(util.angle(2.3456e7*self.e1, 3.4567e-6*self.e1), 0.0)

    def testAngleRandom(self):
        for x in numpy.random.uniform(0, pi, 20):
            r = numpy.random.uniform(0,1000)
            v = r * numpy.array([cos(x), sin(x), 0])
            assert_almost_equal(util.angle(self.e1, v), x, 6)

    def testNorm(self):
        assert_equal(util.norm(self.e3), 1)
        assert_equal(util.norm(self.a), numpy.linalg.norm(self.a))

    def testNormNullVector(self):
        assert_equal(util.norm(self.null), 0.0)

    def testNormRandom(self):
        for x in numpy.random.uniform(0, pi, 20):
            r = numpy.random.uniform(0,1000)
            v = r * numpy.array([cos(x), sin(x), 0])
            assert_almost_equal(util.norm(v), r, 6)

    def testNormal(self):
        assert_equal(util.normal(self.e1, self.e2), self.e3)
        # add more non-trivial tests

    def testNormalNullVector(self):
        assert_equal(util.normal(self.e1, self.null), 0.0)

    def testStp(self):
        assert_equal(util.stp(self.e1, self.e2, self.e3), 1.0)
        # add more non-trivial tests

    def testDihedral(self):
        ab = self.e1
        bc = ab + self.e2
        cd = bc + self.e3
        assert_almost_equal(util.dihedral(ab, bc, cd), -pi/2)
