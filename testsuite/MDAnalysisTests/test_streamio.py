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
from numpy.testing import *

import MDAnalysis
import MDAnalysis.lib.util as util
import MDAnalysis.tests.datafiles as datafiles
from MDAnalysisTests.test_coordinates import RefAdKSmall
from MDAnalysisTests.plugins.knownfailure import knownfailure

import StringIO
import cStringIO
import tempfile
import os


class TestIsstream(TestCase):
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
        obj = (i for i in xrange(3))
        assert_equal(util.isstream(obj), False)

    def test_file(self):
        with open(datafiles.PSF) as obj:
            assert_equal(util.isstream(obj), True)

    def test_cStringIO_read(self):
        with open(datafiles.PSF, "r") as f:
            obj = cStringIO.StringIO(f.read())
        assert_equal(util.isstream(obj), True)
        obj.close()

    def test_cStringIO_write(self):
        obj = cStringIO.StringIO()
        assert_equal(util.isstream(obj), True)
        obj.close()

    def test_StringIO_read(self):
        with open(datafiles.PSF, "r") as f:
            obj = StringIO.StringIO(f)
        assert_equal(util.isstream(obj), True)
        obj.close()

    def test_StringIO_write(self):
        obj = StringIO.StringIO()
        assert_equal(util.isstream(obj), True)
        obj.close()


class TestNamedStream(TestCase):
    def setUp(self):
        self.filename = datafiles.PSF
        self.numlines = 12326  # len(open(self.filename).readlines())
        self.text = [
            "The Jabberwock, with eyes of flame,\n",
            "Came whiffling through the tulgey wood,\n",
            "And burbled as it came!"]
        self.textname = "jabberwock.txt"
        self.numtextlines = len(self.text)

    def testClosing(self):
        obj = cStringIO.StringIO("".join(self.text))
        ns = util.NamedStream(obj, self.textname, close=True)
        assert_equal(ns.closed, False)
        ns.close()
        assert_equal(ns.closed, True)

    def testClosingForce(self):
        obj = cStringIO.StringIO("".join(self.text))
        ns = util.NamedStream(obj, self.textname)
        assert_equal(ns.closed, False)
        ns.close()
        assert_equal(ns.closed, False)
        ns.close(force=True)
        assert_equal(ns.closed, True)

    def testcStringIO_read(self):
        obj = cStringIO.StringIO("".join(self.text))
        ns = util.NamedStream(obj, self.textname)
        assert_equal(ns.name, self.textname)
        assert_equal(str(ns), self.textname)
        assert_equal(len(ns.readlines()), self.numtextlines)
        ns.reset()
        assert_equal(len(ns.readlines()), self.numtextlines)
        ns.close(force=True)

    def testFile_read(self):
        obj = open(self.filename, 'r')
        ns = util.NamedStream(obj, self.filename)
        assert_equal(ns.name, self.filename)
        assert_equal(str(ns), self.filename)
        assert_equal(len(ns.readlines()), self.numlines)
        ns.reset()
        assert_equal(len(ns.readlines()), self.numlines)
        ns.close(force=True)

    def testcStringIO_write(self):
        obj = cStringIO.StringIO()
        ns = util.NamedStream(obj, self.textname)
        ns.writelines(self.text)
        assert_equal(ns.name, self.textname)
        assert_equal(str(ns), self.textname)
        ns.reset()
        assert_equal(len(ns.readlines()), len(self.text))
        ns.reset()
        assert_equal(ns.read(20), "".join(self.text)[:20])
        ns.close(force=True)

    def testFile_write(self):
        fd, outfile = tempfile.mkstemp(suffix=".txt")
        os.close(fd)
        try:
            obj = open(outfile, "w")
            ns = util.NamedStream(obj, outfile, close=True)
            ns.writelines(self.text)
            ns.close()
            text = open(outfile).readlines()

            assert_equal(ns.name, outfile)
            assert_equal(str(ns), outfile)
            assert_equal(len(text), len(self.text))
            assert_equal("".join(text), "".join(self.text))
        finally:
            ns.close()
            obj.close()
            try:
                os.unlink(outfile)
            except OSError:
                pass


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
        self.buffers = dict(
            (name, "".join(open(fn).readlines())) for name, fn in self.filenames.iteritems())
        self.filenames['XYZ_PSF'] = "bogus/path/mini.psf"
        self.buffers['XYZ_PSF'] = """\
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
        return StringIO.StringIO(self.buffers[name])

    def as_cStringIO(self, name):
        return cStringIO.StringIO(self.buffers[name])

    def as_NamedStream(self, name):
        return util.NamedStream(self.as_cStringIO(name), self.filenames[name])


streamData = _StreamData()
del _StreamData


# possibly add tests to individual readers instead?
class TestStreamIO(TestCase, RefAdKSmall):
    def test_PrimitivePDBReader(self):
        u = MDAnalysis.Universe(streamData.as_NamedStream('PDB'))
        assert_equal(u.atoms.n_atoms, self.ref_n_atoms)

    @knownfailure()
    def test_PDBReader(self):
        try:
            u = MDAnalysis.Universe(streamData.as_NamedStream('PDB'), permissive=False)
        except Exception as err:
            raise AssertionError("StreamIO not supported:\n>>>>> {0}".format(err))
        assert_equal(u.atoms.n_atoms, self.ref_n_atoms)

    def test_CRDReader(self):
        u = MDAnalysis.Universe(streamData.as_NamedStream('CRD'))
        assert_equal(u.atoms.n_atoms, self.ref_n_atoms)

    def test_PSFParser(self):
        u = MDAnalysis.Universe(streamData.as_NamedStream('PSF'))
        assert_equal(u.atoms.n_atoms, self.ref_n_atoms)

    def test_PSF_CRD(self):
        u = MDAnalysis.Universe(streamData.as_NamedStream('PSF'),
                                streamData.as_NamedStream('CRD'))
        assert_equal(u.atoms.n_atoms, self.ref_n_atoms)

    def test_PQRReader(self):
        u = MDAnalysis.Universe(streamData.as_NamedStream('PQR'))
        assert_equal(u.atoms.n_atoms, self.ref_n_atoms)
        assert_almost_equal(u.atoms.total_charge(), self.ref_charmm_totalcharge, 3,
                            "Total charge (in CHARMM) does not match expected value.")
        assert_almost_equal(u.atoms.H.charges, self.ref_charmm_Hcharges, 3,
                            "Charges for H atoms do not match.")

    def test_PDBQTReader(self):
        u = MDAnalysis.Universe(streamData.as_NamedStream('PDBQT'))
        sel = u.select_atoms('backbone')
        assert_equal(sel.n_atoms, 796)
        sel = u.select_atoms('segid A')
        assert_equal(sel.n_atoms, 909, "failed to select segment A")
        sel = u.select_atoms('segid B')
        assert_equal(sel.n_atoms, 896, "failed to select segment B")

    def test_GROReader(self):
        u = MDAnalysis.Universe(streamData.as_NamedStream('GRO'))
        assert_equal(u.atoms.n_atoms, 6)
        assert_almost_equal(u.atoms[3].position,
                            10. * np.array([1.275, 0.053, 0.622]), 3,  # manually convert nm -> A
                            err_msg="wrong coordinates for water 2 OW")
        assert_almost_equal(u.atoms[3].velocity,
                            10. * np.array([0.2519, 0.3140, -0.1734]), 3,  # manually convert nm/ps -> A/ps
                            err_msg="wrong velocity for water 2 OW")

    def test_MOL2Reader(self):
        u = MDAnalysis.Universe(streamData.as_NamedStream('MOL2'))
        assert_equal(len(u.atoms), 49)
        assert_equal(u.trajectory.n_frames, 200)
        u.trajectory[199]
        assert_array_almost_equal(u.atoms.positions[0], [1.7240, 11.2730, 14.1200])

    def test_XYZReader(self):
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
