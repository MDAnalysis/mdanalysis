# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import print_function

import MDAnalysis
import MDAnalysis.analysis.align as align
import MDAnalysis.analysis.rms as rms
from MDAnalysis import SelectionError

from numpy.testing import (TestCase, dec,
                           assert_almost_equal, assert_raises, assert_equal)
import numpy as np
from nose.plugins.attrib import attr

from os import path

from MDAnalysisTests.datafiles import PSF, DCD, FASTA
from MDAnalysisTests import executable_not_found, parser_not_found, tempdir


class TestRotationMatrix(object):
    def __init__(self):
        self.a = np.array([[0.1, 0.2, 0.3],
                           [1.1, 1.1, 1.1]])
        self.b = np.array([[0.1, 0.1, 0.1],
                           [1.1, 1.1, 1.1]])
        self.w = np.array([1.3, 2.3])

    def test_no_solution_no_weights(self):
        rot, rmsd = align.rotation_matrix(self.a, self.b)
        assert_equal(rot, np.eye(3))
        assert_equal(rmsd, None)

    def test_no_solution_with_weights(self):
        rot, rmsd = align.rotation_matrix(self.a, self.b, self.w)
        assert_equal(rot, np.eye(3))
        assert_equal(rmsd, None)

    def test_wrong_dtype(self):
        rot, rmsd = align.rotation_matrix(self.a.astype(np.int),
                                          self.b.astype(np.int),
                                          self.w.astype(np.float32))
        assert_equal(rot, np.eye(3))
        assert_equal(rmsd, None)

    @staticmethod
    def test_list_args():
        a = [[0.1, 0.2, 0.3],
             [1.1, 1.1, 1.1]]
        b = [[0.1, 0.1, 0.1],
             [1.1, 1.1, 1.1]]
        w = [1.3, 2.3]
        rot, rmsd = align.rotation_matrix(a, b, w)
        assert_equal(rot, np.eye(3))
        assert_equal(rmsd, None)


class TestAlign(TestCase):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.reference = MDAnalysis.Universe(PSF, DCD)
        # output is always same as input (=DCD)
        self.tempdir = tempdir.TempDir()
        self.outfile = path.join(self.tempdir.name, 'align_test.dcd')

    def tearDown(self):
        del self.tempdir
        del self.universe
        del self.reference

    def test_rmsd(self):
        self.universe.trajectory[0]  # ensure first frame
        bb = self.universe.select_atoms('backbone')
        first_frame = bb.positions
        self.universe.trajectory[-1]
        last_frame = bb.positions
        assert_almost_equal(rms.rmsd(first_frame, first_frame), 0.0, 5,
                            err_msg="error: rmsd(X,X) should be 0")
        # rmsd(A,B) = rmsd(B,A) should be exact but spurious failures in the
        # 9th decimal have been observed (see Issue 57 comment #1) so we relax
        # the test to 6 decimals.
        rmsd = rms.rmsd(first_frame, last_frame, superposition=True)
        assert_almost_equal(rms.rmsd(last_frame, first_frame,
                                     superposition=True),
                            rmsd, 6,
                            err_msg="error: rmsd() is not symmetric")
        assert_almost_equal(rmsd, 6.820321761927005, 5,
                            err_msg="RMSD calculation between 1st and last "
                            "AdK frame gave wrong answer")

    @dec.slow
    @attr('issue')
    def test_rms_fit_trj(self):
        """Testing align.rms_fit_trj() for all atoms (Issue 58)"""
        # align to *last frame* in target... just for the heck of it
        self.reference.trajectory[-1]
        align.rms_fit_trj(self.universe, self.reference, select="all",
                          filename=self.outfile, quiet=True)
        fitted = MDAnalysis.Universe(PSF, self.outfile)
        # RMSD against the reference frame
        # calculated on Mac OS X x86 with MDA 0.7.2 r689
        # VMD: 6.9378711
        self._assert_rmsd(fitted, 0, 6.929083044751061)
        self._assert_rmsd(fitted, -1, 0.0)

    def _assert_rmsd(self, fitted, frame, desired):
        fitted.trajectory[frame]
        rmsd = rms.rmsd(self.reference.atoms.positions,
                        fitted.atoms.positions, superposition=True)
        assert_almost_equal(rmsd, desired, decimal=5,
                            err_msg="frame {0:d} of fit does not have "
                            "expected RMSD".format(frame))

    @attr('issue')
    def test_alignto_checks_selections(self):
        """Testing that alignto() fails if selections do not
        match (Issue 143)"""
        u = self.universe

        def different_size():
            a = u.atoms[10:100]
            b = u.atoms[10:101]
            return MDAnalysis.analysis.align.alignto(a, b)

        assert_raises(SelectionError, different_size)

        def different_atoms():
            a = u.atoms[10:20]
            b = u.atoms[10:17] + u.atoms[18:21]
            return MDAnalysis.analysis.align.alignto(a, b)

        assert_raises(SelectionError, different_atoms)


class TestAlignmentProcessing(TestCase):
    def setUp(self):
        self.seq = FASTA
        self.tempdir = tempdir.TempDir()
        self.alnfile = path.join(self.tempdir.name, 'alignmentprocessing.aln')
        self.treefile = path.join(self.tempdir.name, 'alignmentprocessing.dnd')

    def tearDown(self):
        del self.tempdir

    @attr('issue')
    def test_fasta2select_aligned(self):
        """test align.fasta2select() on aligned FASTA (Issue 112)"""
        from MDAnalysis.analysis.align import fasta2select

        sel = fasta2select(self.seq, is_aligned=True)
        # length of the output strings, not residues or anything real...
        assert_equal(len(sel['reference']), 30623,
                     err_msg="selection string has unexpected length")
        assert_equal(len(sel['mobile']), 30623,
                     err_msg="selection string has unexpected length")

    @attr('issue')
    @dec.skipif(executable_not_found("clustalw2"),
                msg="Test skipped because clustalw2 executable not found")
    def test_fasta2select_ClustalW(self):
        """MDAnalysis.analysis.align: test fasta2select() with calling
        ClustalW (Issue 113)"""
        # note: will not be run if clustalw is not installed
        from MDAnalysis.analysis.align import fasta2select

        sel = fasta2select(self.seq, is_aligned=False,
                           alnfilename=self.alnfile,
                           treefilename=self.treefile)
        # numbers computed from alignment with clustalw 2.1 on Mac OS X
        # [orbeckst] length of the output strings, not residues or anything
        # real...
        assert_equal(len(sel['reference']), 23080,
                     err_msg="selection string has unexpected length")
        assert_equal(len(sel['mobile']), 23090,
                     err_msg="selection string has unexpected length")
