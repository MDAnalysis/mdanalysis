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
from __future__ import absolute_import, division

from six.moves import range, StringIO
import six

import numpy as np
from numpy.testing import (assert_raises, assert_equal, assert_almost_equal,
                           assert_array_almost_equal, assert_, assert_array_equal)

import MDAnalysis as mda
import MDAnalysis.lib.util as util
import MDAnalysis.lib.mdamath as mdamath
from MDAnalysis.lib.util import cached
from MDAnalysis.core.topologyattrs import Bonds
from MDAnalysis.exceptions import NoDataError


from MDAnalysisTests.datafiles import Make_Whole


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

    def test_parse_residue(self):
        for rstring, residue in self.residues:
            yield check_parse_residue, rstring, residue

    def test_parse_residue_VE(self):
        assert_raises(ValueError, util.parse_residue, 'ZZZ')

    def test_convert_aa_code_long(self):
        for resname1, strings in self.aa:
            for resname3 in strings:
                yield check_convert_aa_3to1, resname3, resname1

    def test_convert_aa_code_short(self):
        for resname1, strings in self.aa:
            yield check_convert_aa_1to3, resname1, strings[0]

    def test_VE_1(self):
        assert_raises(ValueError, util.convert_aa_code, 'XYZXYZ')

    def test_VE_2(self):
        assert_raises(ValueError, util.convert_aa_code, '£')

def test_greedy_splitext(inp="foo/bar/boing.2.pdb.bz2",
                         ref=("foo/bar/boing", ".2.pdb.bz2")):
    root, ext = util.greedy_splitext(inp)
    assert_equal(root, ref[0], err_msg="root incorrect")
    assert_equal(ext, ref[1], err_msg="extension incorrect")

class TestIterable(object):
    def test_lists(self):
        assert_equal(util.iterable([1, 2, 3]), True)
        assert_equal(util.iterable([]), True)

    def test_tuples(self):
        assert_equal(util.iterable((1, 2, 3)), True)
        assert_equal(util.iterable(()), True)

    def test_iterator(self):
        assert_equal(util.iterable(range(3)), True)

    def test_arrays(self):
        assert_equal(util.iterable(np.array([1, 2, 3])), True)

    def test_scalars(self):
        assert_equal(util.iterable(123), False)

    def test_strings(self):
        """Test that iterable() works on any string (Fixed bug)"""
        assert_equal(util.iterable("byte string"), False)
        assert_equal(util.iterable(u"unicode string"), False)


class TestFilename(object):
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
        ns = util.NamedStream(StringIO(), self.filename)
        fn = util.filename(ns, ext=self.ext)
        # assert_equal replace by this if loop to avoid segfault on some systems
        if fn != ns:
            raise AssertionError("fn and ns are different")
        assert_equal(str(fn), self.filename2)
        assert_equal(ns.name, self.filename2)


class TestGeometryFunctions(object):
    def setUp(self):
        self.e1 = np.array([1., 0, 0])
        self.e2 = np.array([0, 1., 0])
        self.e3 = np.array([0, 0, 1.])
        self.a = np.array([np.cos(np.pi / 3), np.sin(np.pi / 3), 0])
        self.null = np.zeros(3)

    def test_AngleUnitvectors(self):
        assert_equal(mdamath.angle(self.e1, self.e2), np.pi / 2)
        assert_equal(mdamath.angle(self.e1, self.a), np.pi / 3)

    def test_AngleVectors(self):
        assert_equal(mdamath.angle(2 * self.e1, self.e2), np.pi / 2)
        assert_equal(mdamath.angle(-2 * self.e1, self.e2), np.pi - np.pi / 2)
        assert_equal(mdamath.angle(23.3 * self.e1, self.a), np.pi / 3)

    def test_AngleNullVector(self):
        assert_equal(mdamath.angle(self.e1, self.null), np.nan)

    def test_AngleColinear(self):
        assert_equal(mdamath.angle(self.a, self.a), 0.0)

    def test_AnglePi(self):
        assert_almost_equal(mdamath.angle(-2.3456e7 * self.e1, 3.4567e-6 * self.e1), np.pi)
        assert_almost_equal(mdamath.angle(2.3456e7 * self.e1, 3.4567e-6 * self.e1), 0.0)

    def _check_AngleRange(self, x):
        r = 1000.
        v = r * np.array([np.cos(x), np.sin(x), 0])
        assert_almost_equal(mdamath.angle(self.e1, v), x, 6)

    def test_AngleRange(self):
        for x in np.linspace(0, np.pi, 20):
            yield self._check_AngleRange, x

    def test_Norm(self):
        assert_equal(mdamath.norm(self.e3), 1)
        assert_equal(mdamath.norm(self.a), np.linalg.norm(self.a))

    def test_NormNullVector(self):
        assert_equal(mdamath.norm(self.null), 0.0)

    def _check_NormRange(self, x):
        r = 1000.
        v = r * np.array([np.cos(x), np.sin(x), 0])
        assert_almost_equal(mdamath.norm(v), r, 6)

    def test_NormRange(self):
        for x in np.linspace(0, np.pi, 20):
            yield self._check_NormRange, x

    def test_Normal(self):
        assert_equal(mdamath.normal(self.e1, self.e2), self.e3)
        # add more non-trivial tests

    def testNormalNullVector(self):
        assert_equal(mdamath.normal(self.e1, self.null), 0.0)

    def testStp(self):
        assert_equal(mdamath.stp(self.e1, self.e2, self.e3), 1.0)
        # add more non-trivial tests

    def testDihedral(self):
        ab = self.e1
        bc = ab + self.e2
        cd = bc + self.e3
        assert_almost_equal(mdamath.dihedral(ab, bc, cd), -np.pi / 2)

class TestMakeWhole(object):
    """Set up a simple system:

    +-----------+
    |           |
    | 6       3 | 6
    | !       ! | !
    |-5-8   1-2-|-5-8
    | !       ! | !
    | 7       4 | 7
    |           |
    +-----------+
    """
    def setUp(self):
        self.u = mda.Universe(Make_Whole)
        self.ag = self.u.residues[0].atoms

    def tearDown(self):
        del self.u

    def test_no_bonds(self):
        # NoData caused by no bonds
        assert_raises(NoDataError, mdamath.make_whole, self.ag)

    def test_not_orthogonal(self):
        # Not an orthogonal unit cell
        self._load_bonds()

        self.u.dimensions = [10., 10., 10., 80., 80., 80]
        assert_raises(ValueError, mdamath.make_whole, self.ag)

    def test_zero_box_size(self):
        self._load_bonds()

        self.u.dimensions = [0., 0., 0., 90., 90., 90.]
        assert_raises(ValueError, mdamath.make_whole, self.ag)

    def test_too_small_box_size(self):
        self._load_bonds()

        # Set the z dimensions to 0.5, which is small compared to the
        # bonds (1-2)
        self.u.dimensions = [100.0, 100.0, 0.5, 90., 90., 90.]
        assert_raises(ValueError, mdamath.make_whole, self.ag)

    def _load_bonds(self):
        # Load some bonds into Universe, not required for all tests
        bondlist = [(0, 1), (1, 2), (1, 3), (1, 4), (4, 5), (4, 6), (4, 7)]
        self.u.add_TopologyAttr(Bonds(bondlist))

    def test_wrong_reference_atom(self):
        # Reference atom not in atomgroup
        self._load_bonds()
        assert_raises(ValueError, mdamath.make_whole, self.ag,
                      reference_atom=self.u.atoms[-1])

    def test_impossible_solve(self):
        # check that the algorithm sees the bad walk
        self._load_bonds()
        assert_raises(ValueError, mdamath.make_whole, self.u.atoms)

    def test_walk_1(self):
        self._load_bonds()
        # self.ag is contiguous
        assert_(mdamath._is_contiguous(
            self.ag, self.u.residues[0].atoms[0]))

    def test_walk_2(self):
        self._load_bonds()
        # u.atoms isnt all contiguous
        assert_(not mdamath._is_contiguous(
            self.u.atoms, self.u.residues[0].atoms[0]))

    def test_solve_1(self):
        # regular usage of function
        self._load_bonds()

        refpos = self.u.atoms[:4].positions.copy()

        mdamath.make_whole(self.ag)

        assert_array_almost_equal(self.u.atoms[:4].positions, refpos)
        assert_array_almost_equal(self.u.atoms[4].position,
                                  np.array([110.0, 50.0, 0.0]))
        assert_array_almost_equal(self.u.atoms[5].position,
                                  np.array([110.0, 60.0, 0.0]))
        assert_array_almost_equal(self.u.atoms[6].position,
                                  np.array([110.0, 40.0, 0.0]))
        assert_array_almost_equal(self.u.atoms[7].position,
                                  np.array([120.0, 50.0, 0.0]))

    def test_solve_2(self):
        # use but specify the center atom
        self._load_bonds()

        refpos = self.u.atoms[4:8].positions.copy()

        mdamath.make_whole(self.ag,
                           reference_atom=self.u.residues[0].atoms[4])

        assert_array_almost_equal(self.u.atoms[4:8].positions, refpos)
        assert_array_almost_equal(self.u.atoms[0].position,
                                  np.array([-20.0, 50.0, 0.0]))
        assert_array_almost_equal(self.u.atoms[1].position,
                                  np.array([-10.0, 50.0, 0.0]))
        assert_array_almost_equal(self.u.atoms[2].position,
                                  np.array([-10.0, 60.0, 0.0]))
        assert_array_almost_equal(self.u.atoms[3].position,
                                  np.array([-10.0, 40.0, 0.0]))

    def test_solve_3(self):
        # put in a chunk that doesn't need any work
        self._load_bonds()

        refpos = self.u.atoms[:1].positions.copy()

        mdamath.make_whole(self.u.atoms[:1])

        assert_array_almost_equal(self.u.atoms[:1].positions, refpos)

    def test_solve_4(self):
        # Put in only some of a fragment,
        # check that not everything gets moved
        self._load_bonds()

        chunk = self.u.atoms[:7]
        refpos = self.u.atoms[7].position.copy()

        mdamath.make_whole(chunk)

        assert_array_almost_equal(self.u.atoms[7].position, refpos)
        assert_array_almost_equal(self.u.atoms[4].position,
                                  np.array([110.0, 50.0, 0.0]))
        assert_array_almost_equal(self.u.atoms[5].position,
                                  np.array([110.0, 60.0, 0.0]))
        assert_array_almost_equal(self.u.atoms[6].position,
                                  np.array([110.0, 40.0, 0.0]))

    def test_double_frag_short_bonds(self):
        # previous bug where if two fragments are given
        # but all bonds were short, the algorithm didn't
        # complain
        self._load_bonds()
        mdamath.make_whole(self.ag)
        assert_raises(ValueError, mdamath.make_whole, self.u.atoms)


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


class TestCachedDecorator(object):
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


class TestConvFloat(object):
    def test_float_1(self):
        assert_equal(util.conv_float('0.45'), 0.45)

    def test_float_2(self):
        assert_equal(util.conv_float('.45'), 0.45)

    def test_str(self):
        assert_equal(util.conv_float('a.b'), 'a.b')

    def test_map_1(self):
        ret = [util.conv_float(el) for el in ('0.45', '0.56', '6.7')]
        assert_equal(ret, [0.45, 0.56, 6.7])

    def test_map_2(self):
        ret = [util.conv_float(el) for el in ('0.45', 'a.b', '!!')]
        assert_equal(ret, [0.45, 'a.b', '!!'])

class TestFixedwidthBins(object):
    def test_keys(self):
        ret = util.fixedwidth_bins(0.5, 1.0, 2.0)
        for k in ['Nbins', 'delta', 'min', 'max']:
            assert_(k in ret)

    def test_VE(self):
        assert_raises(ValueError, util.fixedwidth_bins, 0.1, 5.0, 4.0)

    def test_usage_1(self):
        ret = util.fixedwidth_bins(0.1, 4.0, 5.0)
        assert_equal(ret['Nbins'], 10)
        assert_equal(ret['delta'], 0.1)
        assert_equal(ret['min'], 4.0)
        assert_equal(ret['max'], 5.0)

    def test_usage_2(self):
        ret = util.fixedwidth_bins(0.4, 4.0, 5.0)
        assert_equal(ret['Nbins'], 3)
        assert_equal(ret['delta'], 0.4)
        assert_almost_equal(ret['min'], 3.9)
        assert_almost_equal(ret['max'], 5.1)

class TestGuessFormat(object):
    """Test guessing of format from filenames

    Tests also getting the appropriate Parser and Reader from a
    given filename
    """
    # list of known formats, followed by the desired Parser and Reader
    # None indicates that there isn't a Parser/Reader for this format
    formats = [
        ('CHAIN', None, mda.coordinates.chain.ChainReader),
        ('CONFIG', mda.topology.DLPolyParser.ConfigParser, mda.coordinates.DLPoly.ConfigReader),
        ('CRD', mda.topology.CRDParser.CRDParser, mda.coordinates.CRD.CRDReader),
        ('DMS', mda.topology.DMSParser.DMSParser, mda.coordinates.DMS.DMSReader),
        ('GMS', mda.topology.GMSParser.GMSParser, mda.coordinates.GMS.GMSReader),
        ('GRO', mda.topology.GROParser.GROParser, mda.coordinates.GRO.GROReader),
        ('HISTORY', mda.topology.DLPolyParser.HistoryParser, mda.coordinates.DLPoly.HistoryReader),
        ('INPCRD', None, mda.coordinates.INPCRD.INPReader),
        ('MDCRD', None, mda.coordinates.TRJ.TRJReader),
        ('MMTF', mda.topology.MMTFParser.MMTFParser, mda.coordinates.MMTF.MMTFReader),
        ('MOL2', mda.topology.MOL2Parser.MOL2Parser, mda.coordinates.MOL2.MOL2Reader),
        ('NC', None, mda.coordinates.TRJ.NCDFReader),
        ('NCDF', None, mda.coordinates.TRJ.NCDFReader),
        ('PDB', mda.topology.PDBParser.PDBParser, mda.coordinates.PDB.PDBReader),
        ('PDBQT', mda.topology.PDBQTParser.PDBQTParser, mda.coordinates.PDBQT.PDBQTReader),
        ('PRMTOP', mda.topology.TOPParser.TOPParser, None),
        ('PQR', mda.topology.PQRParser.PQRParser, mda.coordinates.PQR.PQRReader),
        ('PSF', mda.topology.PSFParser.PSFParser, None),
        ('RESTRT', None, mda.coordinates.INPCRD.INPReader),
        ('TOP', mda.topology.TOPParser.TOPParser, None),
        ('TPR', mda.topology.TPRParser.TPRParser, None),
        ('TRJ', None, mda.coordinates.TRJ.TRJReader),
        ('TRR', None, mda.coordinates.TRR.TRRReader),
        ('XML', mda.topology.HoomdXMLParser.HoomdXMLParser, None),
        ('XPDB', mda.topology.ExtendedPDBParser.ExtendedPDBParser, mda.coordinates.PDB.ExtendedPDBReader),
        ('XTC', None, mda.coordinates.XTC.XTCReader),
        ('XYZ', mda.topology.XYZParser.XYZParser, mda.coordinates.XYZ.XYZReader),
    ]
    if six.PY2:
        # DCD, LAMMPS, and TRZ are not supported on Python 3 yet
        formats += [
            ('DATA', mda.topology.LAMMPSParser.DATAParser,
             mda.coordinates.LAMMPS.DATAReader),
            ('DCD', None, mda.coordinates.DCD.DCDReader),
            ('LAMMPS', None, mda.coordinates.LAMMPS.DCDReader),
            ('TRZ', None, mda.coordinates.TRZ.TRZReader),
        ]
    # list of possible compressed extensions
    # include no extension too!
    compressed_extensions = ['.bz2', '.gz']

    def _check_get_ext(self, f, fn):
        """Check that get_ext works"""
        a, b = util.get_ext(fn)

        assert_equal(a, 'file')
        assert_equal(b, f.lower())

    def _check_compressed(self, f, fn):
        """Check that format suffixed by compressed extension works"""
        a = util.format_from_filename_extension(fn)
        # expect answer to always be uppercase
        assert_equal(a, f.upper())

    def _check_guess_format(self, f, fn):
        a = util.guess_format(fn)
        # expect answer to always be uppercase
        assert_equal(a, f.upper())

    def _check_get_parser(self, fn, P):
        a = mda.topology.core.get_parser_for(fn)

        assert_equal(a, P)

    def _check_get_parser_invalid(self, fn):
        assert_raises(ValueError, mda.topology.core.get_parser_for, fn)

    def _check_get_reader(self, fn, R):
        a = mda.coordinates.core.get_reader_for(fn)

        assert_equal(a, R)

    def _check_get_reader_invalid(self, fn):
        assert_raises(ValueError, mda.coordinates.core.get_reader_for, fn)

    def test_formats(self):
        # f - format extension
        # P - parser class or None
        # R - reader class or None
        for form, P, R in self.formats:
            # should work with either lower or upper case extension
            for f in [form.upper(), form.lower()]:
                fn = 'file.{0}'.format(f)
                # check f doesn't trip up get_ext or guess_format
                yield self._check_get_ext, f, fn
                yield self._check_guess_format, f, fn

                # check adding extension to f
                # also checks f without extension
                yield self._check_compressed, f, fn
                for e in self.compressed_extensions:
                    yield self._check_compressed, f, fn + e
                    yield self._check_guess_format, f, fn + e

            # Check that expected parser is returned
            if P is not None:
                yield self._check_get_parser, fn, P
                for e in self.compressed_extensions:
                    yield self._check_get_parser, fn + e, P
            else:
                yield self._check_get_parser_invalid, fn

            # Check that expected reader is returned
            if R is not None:
                yield self._check_get_reader, fn, R
                for e in self.compressed_extensions:
                    yield self._check_get_reader, fn + e, R
            else:
                yield self._check_get_reader_invalid, fn

    def test_check_compressed_format_TE(self):
        assert_raises(TypeError, util.check_compressed_format, 1234, 'bz2')

    def test_format_from_filename_TE(self):
        assert_raises(TypeError, util.format_from_filename_extension, 1234)

    def test_guess_format_stream_VE(self):
        # This stream has no name, so can't guess format
        s = StringIO('this is a very fun file')

        assert_raises(ValueError, util.guess_format, s)

    def test_from_ndarray(self):
        fn = np.zeros((3, 3))
        rd = mda.coordinates.core.get_reader_for(fn)
        assert_equal(rd, mda.coordinates.memory.MemoryReader)


class TestUniqueRows(object):
    def test_unique_rows_2(self):
        a = np.array([[0, 1], [1, 2], [2, 1], [0, 1], [0, 1], [2, 1]])

        assert_array_equal(util.unique_rows(a),
                           np.array([[0, 1], [1, 2], [2, 1]]))

    def test_unique_rows_3(self):
        a = np.array([[0, 1, 2], [0, 1, 2], [2, 3, 4], [0, 1, 2]])

        assert_array_equal(util.unique_rows(a),
                           np.array([[0, 1, 2], [2, 3, 4]]))

    def test_unique_rows_with_view(self):
        # unique_rows doesn't work when flags['OWNDATA'] is False,
        # happens when second dimension is created through broadcast
        a = np.array([1, 2])

        assert_array_equal(util.unique_rows(a[None, :]),
                           np.array([[1, 2]]))


class TestGetWriterFor(object):
    def test_no_filename_argument(self):
        assert_raises(TypeError, mda.coordinates.core.get_writer_for)
        # Does ``get_writer_for`` fails as expected when provided no
        # filename arguments

    def test_precedence(self):
        writer = mda.coordinates.core.get_writer_for('test.pdb', 'GRO')
        assert_equal(writer, mda.coordinates.GRO.GROWriter)
        # Make sure ``get_writer_for`` uses *format* if provided

    def test_missing_extension(self):
        assert_raises(TypeError, mda.coordinates.core.get_writer_for,
                      filename='test', format=None)
        # Make sure ``get_writer_for`` behave as expected if *filename*
        # has no extension

    def test_wrong_format(self):
        assert_raises(TypeError, mda.coordinates.core.get_writer_for,
                      filename="fail_me", format='UNK')
        # Make sure ``get_writer_for`` fails if the format is unknown

    def test_compressed_extension(self):
        for ext in ('.gz', '.bz2'):
            fn = 'test.gro' + ext
            writter = mda.coordinates.core.get_writer_for(filename=fn)
            assert_equal(writter, mda.coordinates.GRO.GROWriter)
        # Make sure ``get_writer_for`` works with compressed file file names

    def test_compressed_extension_fail(self):
        for ext in ('.gz', '.bz2'):
            fn = 'test.unk' + ext
            assert_raises(TypeError, mda.coordinates.core.get_writer_for,
                          filename=fn)
        # Make sure ``get_writer_for`` fails if an unknown format is compressed

    def test_non_string_filename(self):
        assert_raises(ValueError, mda.coordinates.core.get_writer_for,
                      filename=StringIO(), format=None)
        # Does ``get_writer_for`` fails with non string filename, no format

    def test_multiframe_failure(self):
        assert_raises(TypeError, mda.coordinates.core.get_writer_for,
                      filename="fail_me", format='UNK', multiframe=True)
        assert_raises(TypeError, mda.coordinates.core.get_writer_for,
                      filename="fail_me", format='UNK', multiframe=False)
        # does ``get_writer_for`` fail with invalid format and multiframe not None

    def test_multiframe_nonsense(self):
        assert_raises(ValueError, mda.coordinates.core.get_writer_for,
                      filename='this.gro', multiframe='sandwich')

    @staticmethod
    def _check_singleframe(fmt, cls):
        assert_equal(mda.coordinates.core.get_writer_for('this', format=fmt, multiframe=False),
                     cls)

    @staticmethod
    def _check_singleframe_fails(fmt):
        assert_raises(TypeError,
                      mda.coordinates.core.get_writer_for,
                      'this', format=fmt, multiframe=False)

    @staticmethod
    def _check_multiframe(fmt, cls):
        assert_equal(mda.coordinates.core.get_writer_for('this', format=fmt, multiframe=True),
                     cls)

    @staticmethod
    def _check_multiframe_fails(fmt):
        assert_raises(TypeError,
                      mda.coordinates.core.get_writer_for,
                      'this', format=fmt, multiframe=True)

    formats = [
        # format name, related class, singleframe, multiframe
        ('CRD', mda.coordinates.CRD.CRDWriter, True, False),
        #('ENT', mda.coordinates.PDB.PDBWriter, True, False),
        ('GRO', mda.coordinates.GRO.GROWriter, True, False),
        ('MOL2', mda.coordinates.MOL2.MOL2Writer, True, True),
        ('NCDF', mda.coordinates.TRJ.NCDFWriter, True, True),
        ('NULL', mda.coordinates.null.NullWriter, True, True),
        # ('PDB', mda.coordinates.PDB.PDBWriter, True, True), special case, done separately
        ('PDBQT', mda.coordinates.PDBQT.PDBQTWriter, True, False),
        ('PQR', mda.coordinates.PQR.PQRWriter, True, False),
        ('TRR', mda.coordinates.TRR.TRRWriter, True, True),
        ('XTC', mda.coordinates.XTC.XTCWriter, True, True),
        ('XYZ', mda.coordinates.XYZ.XYZWriter, True, True),
    ]
    if six.PY2:
        formats += [
        ('DATA', mda.coordinates.LAMMPS.DATAWriter, True, False),
        ('DCD', mda.coordinates.DCD.DCDWriter, True, True),
        ('LAMMPS', mda.coordinates.LAMMPS.DCDWriter, True, True),
        ('TRZ', mda.coordinates.TRZ.TRZWriter, True, True),
    ]
    def test_get_writer_for(self):
        for fmt, cls, singleframe, multiframe in self.formats:
            for f in [fmt.upper(), fmt.lower()]:
                if singleframe:
                    yield self._check_singleframe, f, cls
                else:
                    yield self._check_singleframe_fails, f
                if multiframe:
                    yield self._check_multiframe, f, cls
                else:
                    yield self._check_multiframe_fails, f

    def test_get_writer_for_pdb(self):
        assert_equal(mda.coordinates.core.get_writer_for('this', format='PDB', multiframe=False),
                     mda.coordinates.PDB.PDBWriter)
        assert_equal(mda.coordinates.core.get_writer_for('this', format='PDB', multiframe=True),
                     mda.coordinates.PDB.MultiPDBWriter)
        assert_equal(mda.coordinates.core.get_writer_for('this', format='ENT', multiframe=False),
                     mda.coordinates.PDB.PDBWriter)
        assert_equal(mda.coordinates.core.get_writer_for('this', format='ENT', multiframe=True),
                     mda.coordinates.PDB.MultiPDBWriter)

class TestBlocksOf(object):
    def test_blocks_of_1(self):
        arr = np.arange(16).reshape(4, 4)

        view = util.blocks_of(arr, 1, 1)

        # should return a (4, 1, 1) view
        # ie 4 lots of 1x1
        assert_(view.shape == (4, 1, 1))
        assert_array_almost_equal(view, np.array([[[0]], [[5]], [[10]], [[15]]]))

        # Change my view, check changes are reflected in arr
        view[:] = 1001

        assert_array_almost_equal(arr,
                                  np.array([[1001, 1, 2, 3],
                                            [4, 1001, 6, 7],
                                            [8, 9, 1001, 11],
                                            [12, 13, 14, 1001]]))

    def test_blocks_of_2(self):
        arr = np.arange(16).reshape(4, 4)

        view = util.blocks_of(arr, 2, 2)

        # should return (2, 2, 2)
        assert_(view.shape == (2, 2, 2))
        assert_array_almost_equal(view, np.array([[[0, 1], [4, 5]],
                                                  [[10, 11], [14, 15]]]))

        view[0] = 100
        view[1] = 200

        assert_array_almost_equal(arr,
                                  np.array([[100, 100, 2, 3],
                                            [100, 100, 6, 7],
                                            [8, 9, 200, 200],
                                            [12, 13, 200, 200]]))

    def test_blocks_of_3(self):
        # testing non square array
        arr = np.arange(32).reshape(8, 4)

        view = util.blocks_of(arr, 2, 1)

        assert_(view.shape == (4, 2, 1))

    def test_blocks_of_VE(self):
        arr = np.arange(16).reshape(4, 4)

        assert_raises(ValueError, util.blocks_of, arr, 2, 1)


class TestNamespace(object):
    def setUp(self):
        self.ns = util.Namespace()

    def tearDown(self):
        del self.ns

    def test_getitem(self):
        self.ns.this = 42

        assert_(self.ns['this'] == 42)

    def test_getitem_KE(self):
        assert_raises(KeyError, dict.__getitem__, self.ns, 'this')

    def test_setitem(self):
        self.ns['this'] = 42

        assert_(self.ns['this'] == 42)

    def test_delitem(self):
        self.ns['this'] = 42
        assert_('this' in self.ns)
        del self.ns['this']
        assert_(not ('this' in self.ns))

    def test_delitem_AE(self):
        def deller():
            del self.ns.this
        assert_raises(AttributeError, deller)

    def test_setattr(self):
        self.ns.this = 42

        assert_(self.ns.this == 42)

    def test_getattr(self):
        self.ns['this'] = 42

        assert_(self.ns.this == 42)

    def test_getattr_AE(self):
        assert_raises(AttributeError, getattr, self.ns, 'this')

    def test_delattr(self):
        self.ns['this'] = 42

        assert_('this' in self.ns)
        del self.ns.this
        assert_(not ('this' in self.ns))

    def test_eq(self):
        self.ns['this'] = 42

        ns2 = util.Namespace()
        ns2['this'] = 42

        assert_(self.ns == ns2)

    def test_len(self):
        assert_(len(self.ns) == 0)
        self.ns['this'] = 1
        self.ns['that'] = 2
        assert_(len(self.ns) == 2)

    def test_iter(self):
        self.ns['this'] = 12
        self.ns['that'] = 24
        self.ns['other'] = 48

        seen = []
        for val in self.ns:
            seen.append(val)
        for val in ['this', 'that', 'other']:
            assert_(val in seen)


class TestTruncateInteger(object):
    @staticmethod
    def _check_vals(a, b):
        assert_(util.ltruncate_int(*a) == b)

    def test_ltruncate_int(self):
        for vals, exp in (
                ((1234, 1), 4),
                ((1234, 2), 34),
                ((1234, 3), 234),
                ((1234, 4), 1234),
                ((1234, 5), 1234),
        ):
            yield self._check_vals, vals, exp
