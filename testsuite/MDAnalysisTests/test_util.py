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
import numpy as np
import numpy.random
from numpy.testing import *


from MDAnalysisTests.datafiles import PSF, Make_Whole

import MDAnalysis.lib.util as util
import MDAnalysis.lib.mdamath as mdamath
from MDAnalysis.lib.util import cached
from MDAnalysis.core.topologyobjects import TopologyGroup, Bond
from MDAnalysis.exceptions import NoDataError
import MDAnalysis as mda

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
        assert_equal(util.iterable(np.array([1, 2, 3])), True)

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
        self.e1 = np.array([1., 0, 0])
        self.e2 = np.array([0, 1., 0])
        self.e3 = np.array([0, 0, 1.])
        self.a = np.array([np.cos(np.pi / 3), np.sin(np.pi / 3), 0])
        self.null = np.zeros(3)

    def testAngleUnitvectors(self):
        assert_equal(mdamath.angle(self.e1, self.e2), np.pi / 2)
        assert_equal(mdamath.angle(self.e1, self.a), np.pi / 3)

    def testAngleVectors(self):
        assert_equal(mdamath.angle(2 * self.e1, self.e2), np.pi / 2)
        assert_equal(mdamath.angle(-2 * self.e1, self.e2), np.pi - np.pi / 2)
        assert_equal(mdamath.angle(23.3 * self.e1, self.a), np.pi / 3)

    def testAngleNullVector(self):
        assert_equal(mdamath.angle(self.e1, self.null), np.nan)

    def testAngleColinear(self):
        assert_equal(mdamath.angle(self.a, self.a), 0.0)

    def testAnglePi(self):
        assert_almost_equal(mdamath.angle(-2.3456e7 * self.e1, 3.4567e-6 * self.e1), np.pi)
        assert_almost_equal(mdamath.angle(2.3456e7 * self.e1, 3.4567e-6 * self.e1), 0.0)

    def testAngleRandom(self):
        for x in np.random.uniform(0, np.pi, 20):
            r = np.random.uniform(0, 1000)
            v = r * np.array([np.cos(x), np.sin(x), 0])
            assert_almost_equal(mdamath.angle(self.e1, v), x, 6)

    def testNorm(self):
        assert_equal(mdamath.norm(self.e3), 1)
        assert_equal(mdamath.norm(self.a), np.linalg.norm(self.a))

    def testNormNullVector(self):
        assert_equal(mdamath.norm(self.null), 0.0)

    def testNormRandom(self):
        for x in np.random.uniform(0, np.pi, 20):
            r = np.random.uniform(0, 1000)
            v = r * np.array([np.cos(x), np.sin(x), 0])
            assert_almost_equal(mdamath.norm(v), r, 6)

    def testNormal(self):
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

    def tearDown(self):
        del self.u

    def test_no_bonds(self):
        # NoData caused by no bonds
        assert_raises(NoDataError, mdamath.make_whole, self.u.residues[0])

    def test_not_orthogonal(self):
        # Not an orthogonal unit cell
        self._load_bonds()

        self.u.dimensions = [10., 10., 10., 80., 80., 80]
        assert_raises(ValueError, mdamath.make_whole, self.u.residues[0])

    def test_zero_box_size(self):
        self._load_bonds()

        self.u.dimensions = [0., 0., 0., 90., 90., 90.]
        assert_raises(ValueError, mdamath.make_whole, self.u.residues[0])

    def test_too_small_box_size(self):
        self._load_bonds()

        # Set the z dimensions to 0.5, which is small compared to the
        # bonds (1-2)
        self.u.dimensions = [100.0, 100.0, 0.5, 90., 90., 90.]
        assert_raises(ValueError, mdamath.make_whole, self.u.residues[0])

    def _load_bonds(self):
        # Load some bonds into Universe, not required for all tests
        bondlist = [(0, 1), (1, 2), (1, 3), (1, 4), (4, 5), (4, 6), (4, 7)]
        tg = TopologyGroup.from_indices(bondlist, self.u.atoms, bondclass=Bond)
        self.u.bonds = tg

    def test_wrong_reference_atom(self):
        # Reference atom not in atomgroup
        self._load_bonds()
        assert_raises(ValueError, mdamath.make_whole, self.u.residues[0],
                      reference_atom=self.u.atoms[-1])

    def test_impossible_solve(self):
        # check that the algorithm sees the bad walk
        self._load_bonds()
        assert_raises(ValueError, mdamath.make_whole, self.u.atoms)

    def test_walk_1(self):
        self._load_bonds()

        assert_equal(mdamath._is_contiguous(self.u.residues[0], self.u.residues[0][0]), True)

    def test_walk_2(self):
        self._load_bonds()

        assert_equal(mdamath._is_contiguous(self.u.atoms, self.u.residues[0][0]), False)

    def test_solve_1(self):
        # regular usage of function
        self._load_bonds()

        refpos = self.u.atoms[:4].positions.copy()

        mdamath.make_whole(self.u.residues[0])

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

        mdamath.make_whole(self.u.residues[0], reference_atom=self.u.residues[0][4])

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


class TestConvFloat(object):
    def test_float_1(self):
        assert_equal(util.conv_float('0.45'), 0.45)

    def test_float_2(self):
        assert_equal(util.conv_float('.45'), 0.45)

    def test_str(self):
        assert_equal(util.conv_float('a.b'), 'a.b')

    def test_map_1(self):
        ret = map(util.conv_float, ['0.45', '0.56', '6.7'])
        assert_equal(ret, [0.45, 0.56, 6.7])

    def test_map_2(self):
        ret = map(util.conv_float, ['0.45', 'a.b', '!!'])
        assert_equal(ret, [0.45, 'a.b', '!!'])

class TestFixedwidthBins(object):
    def test_keys(self):
        ret = util.fixedwidth_bins(0.5, 1.0, 2.0)
        for k in ['Nbins', 'delta', 'min', 'max']:
            assert k in ret

    def test_VE(self):
        assert_raises(ValueError, util.fixedwidth_bins, 0.1, 5.0, 4.0)

    def test_usage_1(self):
        ret = util.fixedwidth_bins(0.1, 4.0, 5.0)
        assert ret['Nbins'] == 10
        assert ret['delta'] == 0.1
        assert ret['min'] == 4.0
        assert ret['max'] == 5.0

    def test_usage_2(self):
        ret = util.fixedwidth_bins(0.4, 4.0, 5.0)
        assert ret['Nbins'] == 3
        assert ret['delta'] == 0.4
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
        ('CHAIN', None, mda.coordinates.base.ChainReader),
        ('CONFIG', mda.topology.DLPolyParser.ConfigParser, mda.coordinates.DLPoly.ConfigReader),
        ('CRD', mda.topology.CRDParser.CRDParser, mda.coordinates.CRD.CRDReader),
        ('DATA', mda.topology.LAMMPSParser.DATAParser, mda.coordinates.LAMMPS.DATAReader),
        ('DCD', None, mda.coordinates.DCD.DCDReader),
        ('DMS', mda.topology.DMSParser.DMSParser, mda.coordinates.DMS.DMSReader),
        ('GMS', mda.topology.GMSParser.GMSParser, mda.coordinates.GMS.GMSReader),
        ('GRO', mda.topology.GROParser.GROParser, mda.coordinates.GRO.GROReader),
        ('HISTORY', mda.topology.DLPolyParser.HistoryParser, mda.coordinates.DLPoly.HistoryReader),
        ('INPCRD', None, mda.coordinates.INPCRD.INPReader),
        ('LAMMPS', None, mda.coordinates.LAMMPS.DCDReader),
        ('MDCRD', None, mda.coordinates.TRJ.TRJReader),
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
        ('TRZ', None, mda.coordinates.TRZ.TRZReader),
        ('XML', mda.topology.HoomdXMLParser.HoomdXMLParser, None),
        ('XPDB', mda.topology.ExtendedPDBParser.ExtendedPDBParser, mda.coordinates.PDB.ExtendedPDBReader),
        ('XTC', None, mda.coordinates.XTC.XTCReader),
        ('XYZ', mda.topology.XYZParser.XYZParser, mda.coordinates.XYZ.XYZReader),
    ]
    # list of possible compressed extensions
    # include no extension too!
    compressed_extensions = ['.bz2', '.gz']

    def _check_get_ext(self, f, fn):
        """Check that get ext works"""
        a, b = util.get_ext(fn)

        assert a == 'file'
        assert b == f.lower()

    def _check_compressed(self, f, fn):
        """Check that format suffixed by compressed extension works"""
        a = util.format_from_filename_extension(fn)

        assert a == f

    def _check_guess_format(self, f, fn):
        a = util.guess_format(fn)

        assert a == f

    def _check_get_parser(self, fn, P):
        a = mda.topology.core.get_parser_for(fn)

        assert a == P

    def _check_get_parser_invalid(self, fn):
        assert_raises(ValueError, mda.topology.core.get_parser_for, fn)

    def _check_get_reader(self, fn, R):
        a = mda.coordinates.core.get_reader_for(fn)

        assert a == R

    def _check_get_reader_invalid(self, fn):
        assert_raises(ValueError, mda.coordinates.core.get_reader_for, fn)

    def test_formats(self):
        # f - format extension
        # P - parser class or None
        # R - reader class or None
        for f, P, R in self.formats:
            fn = 'file.{}'.format(f)
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
        s = cStringIO.StringIO('this is a very fun file')

        assert_raises(ValueError, util.guess_format, s)
