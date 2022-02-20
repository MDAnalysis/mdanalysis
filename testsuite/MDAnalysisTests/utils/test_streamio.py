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
from os.path import abspath, basename, dirname, expanduser, normpath, relpath, split, splitext

from io import StringIO

import pytest
import numpy as np
from numpy.testing import assert_equal, assert_almost_equal, assert_array_almost_equal


import MDAnalysis
import MDAnalysis.lib.util as util
import MDAnalysis.tests.datafiles as datafiles
from MDAnalysisTests.coordinates.reference import RefAdKSmall

import os


class TestIsstream(object):
    def test_hasmethod(self):
        obj = "random string"
        assert_equal(util.hasmethod(obj, "rfind"), True)
        assert_equal(util.hasmethod(obj, "bogusXXX"), False)

    def test_string(self):
        obj = datafiles.PSF  # filename
        assert_equal(util.isstream(obj), False)

    def test_list(self):
        obj = [1, 2, 3]
        assert_equal(util.isstream(obj), False)

    def test_iterator(self):
        obj = (i for i in range(3))
        assert_equal(util.isstream(obj), False)

    def test_file(self):
        with open(datafiles.PSF) as obj:
            assert_equal(util.isstream(obj), True)

    def test_StringIO_read(self):
        with open(datafiles.PSF, "r") as f:
            obj = StringIO(f.read())
        assert_equal(util.isstream(obj), True)
        obj.close()

    def test_StringIO_write(self):
        obj = StringIO()
        assert_equal(util.isstream(obj), True)
        obj.close()


class TestNamedStream(object):
    filename = datafiles.PSF
    numlines = 12326  # len(open(self.filename).readlines())
    text = [
        "The Jabberwock, with eyes of flame,\n",
        "Came whiffling through the tulgey wood,\n",
        "And burbled as it came!"
    ]
    textname = "jabberwock.txt"

    def test_closing(self):
        obj = StringIO("".join(self.text))
        ns = util.NamedStream(obj, self.textname, close=True)
        assert_equal(ns.closed, False)
        ns.close()
        assert_equal(ns.closed, True)

    def test_closing_force(self):
        obj = StringIO("".join(self.text))
        ns = util.NamedStream(obj, self.textname)
        assert_equal(ns.closed, False)
        ns.close()
        assert_equal(ns.closed, False)
        ns.close(force=True)
        assert_equal(ns.closed, True)
        # Issue 3386 - calling close again shouldn't raise an error
        ns.close()

    def test_StringIO_read(self):
        obj = StringIO("".join(self.text))
        ns = util.NamedStream(obj, self.textname)
        assert_equal(ns.name, self.textname)
        assert_equal(str(ns), self.textname)
        assert_equal(len(ns.readlines()), len(self.text))
        ns.reset()
        assert_equal(len(ns.readlines()), len(self.text))
        ns.close(force=True)

    def test_File_read(self):
        obj = open(self.filename, 'r')
        ns = util.NamedStream(obj, self.filename)
        assert_equal(ns.name, self.filename)
        assert_equal(str(ns), self.filename)
        assert_equal(len(ns.readlines()), self.numlines)
        ns.reset()
        assert_equal(len(ns.readlines()), self.numlines)
        ns.close(force=True)

    def test_StringIO_write(self):
        obj = StringIO()
        ns = util.NamedStream(obj, self.textname)
        ns.writelines(self.text)
        assert_equal(ns.name, self.textname)
        assert_equal(str(ns), self.textname)
        ns.reset()
        assert_equal(len(ns.readlines()), len(self.text))
        ns.reset()
        assert_equal(ns.read(20), "".join(self.text)[:20])
        ns.close(force=True)

    def test_File_write(self, tmpdir):
        with tmpdir.as_cwd():
            outfile = "lookingglas.txt"
            with open(outfile, "w") as obj:
                with util.NamedStream(obj, outfile, close=True) as ns:
                    assert_equal(ns.name, outfile)
                    assert_equal(str(ns), outfile)
                    ns.writelines(self.text)

            with open(outfile) as fh:
                text = fh.readlines()

            assert_equal(len(text), len(self.text))
            assert_equal("".join(text), "".join(self.text))

    def test_matryoshka(self):
        obj = StringIO()
        ns = util.NamedStream(obj, 'r')
        with pytest.warns(RuntimeWarning):
            ns2 = util.NamedStream(ns, 'f')
        assert not isinstance(ns2.stream, util.NamedStream)
        assert ns2.name == 'f'


class TestNamedStream_filename_behavior(object):
    textname = os.path.join("~", "stories", "jabberwock.txt") # with tilde ~ to test regular expanduser()
    # note: no setUp() because classes with generators would run it
    #       *for each generated test* and we need it for the generator method

    def create_NamedStream(self, name=None):
        if name is None:
            name = self.textname
        obj = StringIO()
        return util.NamedStream(obj, name)

    @pytest.mark.parametrize('func', (
            abspath,
            basename,
            dirname,
            expanduser,
            normpath,
            relpath,
            split,
            splitext
    ))
    def test_func(self, func):
        # - "expandvars" gave Segmentation fault (OS X 10.6, Python 2.7.11 -- orbeckst)
        # - "expanduser" will either return a string if it carried out interpolation
        #   or "will do nothing" and return the NamedStream (see extra test below).
        #   On systems without a user or HOME, it will also do nothing and the test
        #   below will fail.
        ns = self.create_NamedStream()
        fn = self.textname
        reference = func(fn)
        value = func(ns)
        assert_equal(value, reference,
                     err_msg=("os.path.{0}() does not work with "
                              "NamedStream").format(func.__name__))

    def test_join(self, tmpdir, funcname="join"):
        # join not included because of different call signature
        # but added first argument for the sake of it showing up in the verbose
        # nose output
        ns = self.create_NamedStream()
        fn = self.textname
        reference = str(tmpdir.join(fn))
        value = os.path.join(str(tmpdir), ns)
        assert_equal(value, reference,
                     err_msg=("os.path.{0}() does not work with "
                              "NamedStream").format(funcname))

    def test_expanduser_noexpansion_returns_NamedStream(self):
        ns = self.create_NamedStream("de/zipferlack.txt")  # no tilde ~ in name!
        reference = ns.name
        value = os.path.expanduser(ns)
        assert_equal(value, reference,
                     err_msg=("os.path.expanduser() without '~' did not "
                              "return NamedStream --- weird!!"))

    @pytest.mark.skipif("HOME" not in os.environ, reason='It is needed')
    def test_expandvars(self):
        name = "${HOME}/stories/jabberwock.txt"
        ns = self.create_NamedStream(name)
        reference = os.path.expandvars(name)
        value = os.path.expandvars(ns)
        assert_equal(value, reference,
                     err_msg="os.path.expandvars() did not expand HOME")

    def test_expandvars_noexpansion_returns_NamedStream(self):
        ns = self.create_NamedStream() # no $VAR constructs
        reference = ns.name
        value = os.path.expandvars(ns)
        assert_equal(value, reference,
                     err_msg=("os.path.expandvars() without '$VARS' did not "
                              "return NamedStream --- weird!!"))

    def test_add(self):
        ns = self.create_NamedStream()
        try:
            assert_equal(ns + "foo", self.textname + "foo")
        except TypeError:
            raise pytest.fail("NamedStream does not support  "
                                 "string concatenation, NamedStream + str")

    def test_radd(self):
        ns = self.create_NamedStream()
        try:
            assert_equal("foo" + ns, "foo" + self.textname)
        except TypeError:
            raise pytest.fail("NamedStream does not support right "
                                 "string concatenation, str + NamedStream")


class _StreamData(object):
    """Data for StreamIO functions."""
    filenames = {
        'PSF': datafiles.PSF,
        'CRD': datafiles.CRD,
        'PDB': datafiles.PDB_small,
        'PQR': datafiles.PQR,
        'GRO': datafiles.GRO_velocity,
        'MOL2': datafiles.mol2_molecules,
        'PDBQT': datafiles.PDBQT_input,
    }

    def __init__(self):
        self.buffers = {}
        for name, fn in self.filenames.items():
            with open(fn) as filed:
                self.buffers[name] = "".join(filed.readlines())
        self.filenames['XYZ_PSF'] = u"bogus/path/mini.psf"
        self.buffers['XYZ_PSF'] = u"""\
PSF CMAP

      1 !NTITLE
Mini PSF for in memory XYZ

       8 !NATOM
       1 A    380  THR  N    NH1   -0.470000       14.0070           0
       2 A    380  THR  HN   H      0.310000        1.0080           0
       3 A    380  THR  CA   CT1    0.070000       12.0110           0
       4 A    380  THR  CB   CT1    0.140000       12.0110           0
       5 A    380  THR  OG1  OH1   -0.660000       15.9990           0
       6 A    380  THR  CG2  CT3   -0.270000       12.0110           0
       7 A    380  THR  C    C      0.510000       12.0110           0
       8 A    380  THR  O    O     -0.510000       15.9990           0
"""
        self.filenames['XYZ'] = "bogus/path/mini.xyz"
        self.buffers['XYZ'] = """\
8
frame 1
       N     0.93100   17.31800   16.42300
      HN     1.86100   17.06500   16.17100
      CA     0.48600   18.66500   16.14300
      CB     1.65900   19.66600   15.88700
     OG1     2.53100   19.43000   14.75700
     CG2     2.56700   19.70400   17.04500
       C    -0.38500   18.72400   14.93500
       O    -0.22300   17.81000   14.13400
8
frame 2
       N     1.00200   17.11400   16.52100
      HN     1.85100   16.93900   16.02800
      CA     0.45600   18.48700   16.26500
      CB     1.49700   19.58900   16.08900
     OG1     2.38300   19.42200   14.96500
     CG2     2.47300   19.54600   17.26500
       C    -0.31500   18.63800   14.99300
       O    -0.23100   17.83800   14.10800
8
frame 3
       N     0.94000   16.97600   16.44500
      HN     1.85800   16.71700   16.15500
      CA     0.53300   18.34800   16.17400
      CB     1.79500   19.24700   15.93000
     OG1     2.61400   18.84000   14.91900
     CG2     2.54700   19.25800   17.26500
       C    -0.27300   18.58100   14.94400
       O    -0.23800   17.82300   13.97300

"""

    def as_StringIO(self, name):
        return StringIO(self.buffers[name])

    def as_NamedStream(self, name):
        return util.NamedStream(self.as_StringIO(name), self.filenames[name])


@pytest.fixture(scope='module')
def streamData():
    return _StreamData()


# possibly add tests to individual readers instead?
class TestStreamIO(RefAdKSmall):
    def test_PrimitivePDBReader(self, streamData):
        u = MDAnalysis.Universe(streamData.as_NamedStream('PDB'))
        assert_equal(u.atoms.n_atoms, self.ref_n_atoms)

    def test_PDBReader(self, streamData):
        try:
            u = MDAnalysis.Universe(streamData.as_NamedStream('PDB'))
        except Exception as err:
            raise pytest.fail("StreamIO not supported:\n>>>>> {0}".format(err))
        assert_equal(u.atoms.n_atoms, self.ref_n_atoms)

    def test_CRDReader(self, streamData):
        u = MDAnalysis.Universe(streamData.as_NamedStream('CRD'))
        assert_equal(u.atoms.n_atoms, self.ref_n_atoms)

    def test_PSFParser(self, streamData):
        u = MDAnalysis.Universe(streamData.as_NamedStream('PSF'))
        assert_equal(u.atoms.n_atoms, self.ref_n_atoms)

    def test_PSF_CRD(self, streamData):
        u = MDAnalysis.Universe(streamData.as_NamedStream('PSF'),
                                streamData.as_NamedStream('CRD'))
        assert_equal(u.atoms.n_atoms, self.ref_n_atoms)

    def test_PQRReader(self, streamData):
        u = MDAnalysis.Universe(streamData.as_NamedStream('PQR'))
        assert_equal(u.atoms.n_atoms, self.ref_n_atoms)
        assert_almost_equal(u.atoms.total_charge(), self.ref_charmm_totalcharge, 3,
                            "Total charge (in CHARMM) does not match expected value.")
        assert_almost_equal(u.atoms.select_atoms('name H').charges, self.ref_charmm_Hcharges, 3,
                            "Charges for H atoms do not match.")

    def test_PDBQTReader(self, streamData):
        u = MDAnalysis.Universe(streamData.as_NamedStream('PDBQT'))
        sel = u.select_atoms('backbone')
        assert_equal(sel.n_atoms, 796)
        sel = u.select_atoms('segid A')
        assert_equal(sel.n_atoms, 909, "failed to select segment A")
        sel = u.select_atoms('segid B')
        assert_equal(sel.n_atoms, 896, "failed to select segment B")

    def test_GROReader(self, streamData):
        u = MDAnalysis.Universe(streamData.as_NamedStream('GRO'))
        assert_equal(u.atoms.n_atoms, 6)
        assert_almost_equal(u.atoms[3].position,
                            10. * np.array([1.275, 0.053, 0.622]), 3,  # manually convert nm -> A
                            err_msg="wrong coordinates for water 2 OW")
        assert_almost_equal(u.atoms[3].velocity,
                            10. * np.array([0.2519, 0.3140, -0.1734]), 3,  # manually convert nm/ps -> A/ps
                            err_msg="wrong velocity for water 2 OW")

    def test_MOL2Reader(self, streamData):
        u = MDAnalysis.Universe(streamData.as_NamedStream('MOL2'))
        assert_equal(len(u.atoms), 49)
        assert_equal(u.trajectory.n_frames, 200)
        u.trajectory[199]
        assert_array_almost_equal(u.atoms.positions[0], [1.7240, 11.2730, 14.1200])

    def test_XYZReader(self, streamData):
        u = MDAnalysis.Universe(streamData.as_NamedStream('XYZ_PSF'),
                                streamData.as_NamedStream('XYZ'))
        assert_equal(len(u.atoms), 8)
        assert_equal(u.trajectory.n_frames, 3)
        assert_equal(u.trajectory.frame, 0)  # weird, something odd with XYZ reader
        u.trajectory.next()  # (should really only need one next()... )
        assert_equal(u.trajectory.frame, 1)  # !!!! ???
        u.trajectory.next()  # frame 2
        assert_equal(u.trajectory.frame, 2)
        assert_almost_equal(u.atoms[2].position, np.array([0.45600, 18.48700, 16.26500]), 3,
                            err_msg="wrong coordinates for atom CA at frame 2")
