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

from numpy.testing import TestCase, assert_, assert_equal, assert_array_almost_equal

from MDAnalysis.tests.datafiles import XTC, TRR
import MDAnalysis.coordinates.xdrfile.libxdrfile2 as xdr

class TestLib(TestCase):
    def test_constants(self):
        assert_equal(xdr.DIM, 3, "xdr library not compiled for DIM=3 ?!?")

    def test_xdropen(self):
        XDR = xdr.xdrfile_open(XTC, 'r')
        assert_(XDR is not None, "Failed to open xtc file")
        rc = xdr.xdrfile_close(XDR)
        assert_equal(rc, 0, "Failed to close xtc file")  # this can segfault


class TestXTC(TestCase):
    def test_n_atoms(self):
        natoms = xdr.read_xtc_natoms(XTC)
        assert_equal(natoms, 47681, "Number of atoms in XTC frame")

    def test_n_frames_offsets(self):
        n_frames, offsets = xdr.read_xtc_n_frames(XTC)
        desired_offsets = [0,  165188,  330364,  495520,  660708,  825872,  991044,
                           1156212, 1321384, 1486544]

        assert_equal(n_frames, len(desired_offsets), "Number of frames in XTC trajectory")

        assert_array_almost_equal(offsets, desired_offsets, err_msg="wrong xtc frame offsets")


class TestTRR(TestCase):
    def test_n_atoms(self):
        natoms = xdr.read_trr_natoms(TRR)
        assert_equal(natoms, 47681, "Number of atoms in TRR frame")

    def test_n_frames_offsets(self):
        n_frames, offsets = xdr.read_trr_n_frames(TRR)
        desired_offsets = [0,  1144464,  2288928,  3433392,  4577856,  5722320,
                           6866784,  8011248,  9155712, 10300176]

        assert_equal(n_frames, len(desired_offsets), "Number of frames in TRR trajectory")

        assert_array_almost_equal(offsets, desired_offsets, err_msg="wrong trr frame offsets")

