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
import numpy
import numpy.random
from numpy.testing import *
from numpy import pi, sin, cos

import MDAnalysis.core.util as util
from MDAnalysis.core.util import cached
from MDAnalysis.tests.datafiles import PSF

import cStringIO
import os.path
import tempfile


def check_parse_residue(rstring, residue):
    assert_equal(util.parse_residue(rstring), residue)


def check_convert_aa_3to1(resname3, resname1):
    assert_equal(util.convert_aa_code(resname3), resname1)


def check_convert_aa_1to3(resname1, resname3_canonical):
    assert_equal(util.convert_aa_code(resname1), resname3_canonical)


class TestStringFunctions(object):
    # (1-letter, (canonical 3 letter, other 3/4 letter, ....))
    aa = [
        ('H', ('HIS', 'HISA', 'HISB', 'HSE', 'HSD', 'HIS1', 'HIS2', 'HIE', 'HID')),
        ('K', ('LYS', 'LYSH', 'LYN')),
        ('A', ('ALA',)),
        ('D', ('ASP', 'ASPH', 'ASH')),
        ('E', ('GLU', 'GLUH', 'GLH')),
        ('N', ('ASN',)),
        ('Q', ('GLN',)),
        ('C', ('CYS', 'CYSH', 'CYS1', 'CYS2')),
    ]

    residues = [
        ("LYS300:HZ1", ("LYS", 300, "HZ1")),
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
        assert_equal(util.iterable(numpy.array([1, 2, 3])), True)

    def test_scalars(self):
        assert_equal(util.iterable(123), False)

    def test_strings(self):
        """Test that iterable() works on any string (Fixed bug)"""
        assert_equal(util.iterable("byte string"), False)
        assert_equal(util.iterable(u"unicode string"), False)


class TestFilename(TestCase):
    def setUp(self):
        self.root = "foo"
        self.filename = "foo.psf"
        self.ext = "pdb"
        self.filename2 = "foo.pdb"

    def testStringNoExt(self):
        fn = util.filename(self.filename)
        assert_equal(fn, self.filename)

    def testStringExt(self):
        fn = util.filename(self.filename, ext=self.ext)
        assert_equal(fn, self.filename2)

    def testStringKeep(self):
        fn = util.filename(self.filename, ext=self.ext, keep=True)
        assert_equal(fn, self.filename)

    def testStringRootExt(self):
        fn = util.filename(self.root, ext=self.ext)
        assert_equal(fn, self.filename2)

    def testStringRootExtKeep(self):
        fn = util.filename(self.root, ext=self.ext, keep=True)
        assert_equal(fn, self.filename2)

    def testNamedStream(self):
        ns = util.NamedStream(cStringIO.StringIO(), self.filename)
        fn = util.filename(ns, ext=self.ext)
        # assert_equal replace by this if loop to avoid segfault on some systems
        if fn != ns:
            raise AssertionError("fn and ns are different")
        assert_equal(str(fn), self.filename2)
        assert_equal(ns.name, self.filename2)


class TestGeometryFunctions(TestCase):
    def setUp(self):
        self.e1 = numpy.array([1., 0, 0])
        self.e2 = numpy.array([0, 1., 0])
        self.e3 = numpy.array([0, 0, 1.])
        self.a = numpy.array([cos(pi / 3), sin(pi / 3), 0])
        self.null = numpy.zeros(3)

    def testAngleUnitvectors(self):
        assert_equal(util.angle(self.e1, self.e2), pi / 2)
        assert_equal(util.angle(self.e1, self.a), pi / 3)

    def testAngleVectors(self):
        assert_equal(util.angle(2 * self.e1, self.e2), pi / 2)
        assert_equal(util.angle(-2 * self.e1, self.e2), pi - pi / 2)
        assert_equal(util.angle(23.3 * self.e1, self.a), pi / 3)

    def testAngleNullVector(self):
        assert_equal(util.angle(self.e1, self.null), numpy.nan)

    def testAngleColinear(self):
        assert_equal(util.angle(self.a, self.a), 0.0)

    def testAnglePi(self):
        assert_almost_equal(util.angle(-2.3456e7 * self.e1, 3.4567e-6 * self.e1), pi)
        assert_almost_equal(util.angle(2.3456e7 * self.e1, 3.4567e-6 * self.e1), 0.0)

    def testAngleRandom(self):
        for x in numpy.random.uniform(0, pi, 20):
            r = numpy.random.uniform(0, 1000)
            v = r * numpy.array([cos(x), sin(x), 0])
            assert_almost_equal(util.angle(self.e1, v), x, 6)

    def testNorm(self):
        assert_equal(util.norm(self.e3), 1)
        assert_equal(util.norm(self.a), numpy.linalg.norm(self.a))

    def testNormNullVector(self):
        assert_equal(util.norm(self.null), 0.0)

    def testNormRandom(self):
        for x in numpy.random.uniform(0, pi, 20):
            r = numpy.random.uniform(0, 1000)
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
        assert_almost_equal(util.dihedral(ab, bc, cd), -pi / 2)


class Class_with_Caches(object):
    def __init__(self):
        self._cache = dict()
        self.ref1 = 1.0
        self.ref2 = 2.0
        self.ref3 = 3.0
        self.ref4 = 4.0
        self.ref5 = 5.0

    @cached('val1')
    def val1(self):
        return self.ref1

    # Do one with property decorator as these are used together often
    @property
    @cached('val2')
    def val2(self):
        return self.ref2

    # Check use of property setters
    @property
    @cached('val3')
    def val3(self):
        return self.ref3

    @val3.setter
    def val3(self, new):
        self._clear_caches('val3')
        self._fill_cache('val3', new)

    @val3.deleter
    def val3(self):
        self._clear_caches('val3')

    # Check that args are passed through to underlying functions
    @cached('val4')
    def val4(self, n1, n2):
        return self._init_val_4(n1, n2)

    def _init_val_4(self, m1, m2):
        return self.ref4 + m1 + m2

    # Args and Kwargs
    @cached('val5')
    def val5(self, n, s=None):
        return self._init_val_5(n, s=s)

    def _init_val_5(self, n, s=None):
        return n * s

    # These are designed to mimic the AG and Universe cache methods
    def _clear_caches(self, *args):
        if len(args) == 0:
            self._cache = dict()
        else:
            for name in args:
                try:
                    del self._cache[name]
                except KeyError:
                    pass

    def _fill_cache(self, name, value):
        self._cache[name] = value


class TestCachedDecorator(TestCase):
    def setUp(self):
        self.obj = Class_with_Caches()

    def tearDown(self):
        del self.obj

    def test_val1_lookup(self):
        self.obj._clear_caches()
        assert_equal('val1' in self.obj._cache, False)
        assert_equal(self.obj.val1(), self.obj.ref1)
        ret = self.obj.val1()
        assert_equal('val1' in self.obj._cache, True)
        assert_equal(self.obj._cache['val1'], ret)
        assert_equal(self.obj.val1() is self.obj._cache['val1'], True)

    def test_val1_inject(self):
        # Put something else into the cache and check it gets returned
        # this tests that the cache is blindly being used
        self.obj._clear_caches()
        ret = self.obj.val1()
        assert_equal('val1' in self.obj._cache, True)
        assert_equal(ret, self.obj.ref1)
        new = 77.0
        self.obj._fill_cache('val1', new)
        assert_equal(self.obj.val1(), new)

    # Managed property
    def test_val2_lookup(self):
        self.obj._clear_caches()
        assert_equal('val2' in self.obj._cache, False)
        assert_equal(self.obj.val2, self.obj.ref2)
        ret = self.obj.val2
        assert_equal('val2' in self.obj._cache, True)
        assert_equal(self.obj._cache['val2'], ret)

    def test_val2_inject(self):
        self.obj._clear_caches()
        ret = self.obj.val2
        assert_equal('val2' in self.obj._cache, True)
        assert_equal(ret, self.obj.ref2)
        new = 77.0
        self.obj._fill_cache('val2', new)
        assert_equal(self.obj.val2, new)

        # Setter on cached attribute

    def test_val3_set(self):
        self.obj._clear_caches()
        assert_equal(self.obj.val3, self.obj.ref3)
        new = 99.0
        self.obj.val3 = new
        assert_equal(self.obj.val3, new)
        assert_equal(self.obj._cache['val3'], new)

    def test_val3_del(self):
        # Check that deleting the property removes it from cache,
        self.obj._clear_caches()
        assert_equal(self.obj.val3, self.obj.ref3)
        assert_equal('val3' in self.obj._cache, True)
        del self.obj.val3
        assert_equal('val3' in self.obj._cache, False)
        # But allows it to work as usual afterwards
        assert_equal(self.obj.val3, self.obj.ref3)
        assert_equal('val3' in self.obj._cache, True)

    # Pass args
    def test_val4_args(self):
        self.obj._clear_caches()
        assert_equal(self.obj.val4(1, 2), 1 + 2 + self.obj.ref4)
        # Further calls should yield the old result
        # this arguably shouldn't be cached...
        assert_equal(self.obj.val4(3, 4), 1 + 2 + self.obj.ref4)

    # Pass args and kwargs
    def test_val5_kwargs(self):
        self.obj._clear_caches()
        assert_equal(self.obj.val5(5, s='abc'), 5 * 'abc')

        assert_equal(self.obj.val5(5, s='!!!'), 5 * 'abc')
