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

import MDAnalysis
from numpy.testing import *
from MDAnalysis.tests.datafiles import XTC, TRR

import MDAnalysis.coordinates.xdrfile.libxdrfile as xdr

# FIXES: test_xdropen: error because assert_ not found in numpy < 1.3
# maybe move this into separate module together with
# from numpy.testing import * ?
try:
    from numpy.testing import assert_
except ImportError:
    def assert_(val, msg='') :
        """
        Assert that works in release mode.

        The Python built-in ``assert`` does not work when executing code in
        optimized mode (the ``-O`` flag) - no byte-code is generated for it.

        For documentation on usage, refer to the Python documentation.

        (Code taken from numpy.testing 1.4)
        """
        if not val :
            raise AssertionError(msg)


class TestLib(TestCase):
    def test_constants(self):
        assert_equal(xdr.DIM, 3, "xdr library not compiled for DIM=3 ?!?")

    def test_xdropen(self):
        XDR = xdr.xdrfile_open(XTC, 'r')
        assert_(XDR != None, "Failed top open xtc file")
        rc = xdr.xdrfile_close(XDR)
        assert_equal(rc, 0, "Failed to close xtc file")  # this can segfault


class TestXTC(TestCase):
    def test_numatoms(self):
        natoms = xdr.read_xtc_natoms(XTC)
        assert_equal(natoms, 47681, "Number of atoms in XTC frame")

    def test_numframes(self):
        numframes, offsets = xdr.read_xtc_numframes(XTC)
        assert_equal(numframes, 10, "Number of frames in XTC trajectory")


class TestTRR(TestCase):
    def test_numatoms(self):
        natoms = xdr.read_trr_natoms(TRR)
        assert_equal(natoms, 47681, "Number of atoms in TRR frame")

    def test_numframes(self):
        numframes, offsets = xdr.read_trr_numframes(TRR)
        assert_equal(numframes, 10, "Number of frames in TRR trajectory")

