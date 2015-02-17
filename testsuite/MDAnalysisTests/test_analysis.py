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
from __future__ import print_function
import MDAnalysis
import MDAnalysis.analysis.distances
import MDAnalysis.analysis.align
import MDAnalysis.analysis.hbonds
import MDAnalysis.analysis.helanal
from MDAnalysis import SelectionError, FinishTimeException

from numpy.testing import *
from nose.plugins.attrib import attr

import os
import errno
import tempfile

from MDAnalysis.tests.datafiles import PSF, DCD, FASTA, PDB_helix, PDB_HOLE, XTC_HOLE, GRO, XTC
from . import executable_not_found_runtime


class TestContactMatrix(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.dcd = self.universe.trajectory
        # reasonable precision so that tests succeed on 32 and 64 bit machines
        # (the reference values were obtained on 64 bit)
        # Example:
        #   Items are not equal: wrong maximum distance value
        #   ACTUAL: 52.470254967456412
        #   DESIRED: 52.470257062419059
        self.prec = 5

    def tearDown(self):
        del self.universe
        del self.dcd

    def test_numpy(self):
        U = self.universe
        self.dcd.rewind()
        self.dcd[10]
        # small cutoff value as the input file is a protein
        contacts = MDAnalysis.analysis.distances.contact_matrix(U.atoms.coordinates(), cutoff=1.5, returntype="numpy")
        assert_equal(contacts.shape, (3341, 3341), "wrong shape (should be (Natoms,Natoms))")
        assert_equal(contacts[0][0], True, "first entry should be a contact")
        assert_equal(contacts[0][-1], False, "last entry for first atom should be a non-contact")

    def test_sparse(self):
        U = self.universe
        self.dcd.rewind()
        self.dcd[10]
        # Just taking first 50 atoms as the sparse method is slow
        selection = U.selectAtoms('bynum 1:50')
        # small cutoff value as the input file is a protein
        # High progress_meter_freq so progress meter is not printed during test
        contacts = MDAnalysis.analysis.distances.contact_matrix(selection.coordinates(), cutoff=1.07,
                                                                returntype="sparse", suppress_progmet=True)
        assert_equal(contacts.shape, (50, 50), "wrong shape (should be (50,50))")
        assert_equal(contacts[0, 0], False, "entry (0,0) should be a non-contact")
        assert_equal(contacts[0, 2], True, "entry (0,2) should be a contact")
        assert_equal(contacts[0, 3], True, "entry (0,3) should be a contact")
        assert_equal(contacts[0, 4], False, "entry (0,3) should be a contact")


class TestAlign(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.reference = MDAnalysis.Universe(PSF, DCD)
        fd, self.outfile = tempfile.mkstemp(suffix=".dcd")  # output is always same as input (=DCD)
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.universe
        del self.reference

    def test_rmsd(self):
        self.universe.trajectory[0]  # ensure first frame
        bb = self.universe.selectAtoms('backbone')
        A = bb.coordinates(copy=True)  # coordinates of first frame (copy=True just in case)
        self.universe.trajectory[-1]  # forward to last frame
        B = bb.coordinates()  # coordinates of last frame
        rmsd = MDAnalysis.analysis.align.rmsd(A, B)
        assert_almost_equal(MDAnalysis.analysis.align.rmsd(A, A), 0.0, 5,
                            err_msg="error: rmsd(X,X) should be 0")
        # rmsd(A,B) = rmsd(B,A)  should be exact but spurious failures in the 9th decimal have been
        # observed (see Issue 57 comment #1) so we relax the test to 6 decimals.
        assert_almost_equal(MDAnalysis.analysis.align.rmsd(B, A), rmsd, 6,
                            err_msg="error: rmsd() is not symmetric")
        assert_almost_equal(rmsd, 6.8342494129169804, 5,
                            err_msg="RMSD calculation between 1st and last AdK frame gave wrong answer")

    @dec.slow
    @attr('issue')
    def test_rms_fit_trj(self):
        """Testing align.rms_fit_trj() for all atoms (Issue 58)"""
        # align to *last frame* in target... just for the heck of it
        self.reference.trajectory[-1]
        MDAnalysis.analysis.align.rms_fit_trj(self.universe, self.reference, select="all",
                                              filename=self.outfile, quiet=True)
        fitted = MDAnalysis.Universe(PSF, self.outfile)
        # RMSD against the reference frame
        # calculated on Mac OS X x86 with MDA 0.7.2 r689
        # VMD: 6.9378711
        self._assert_rmsd(fitted, 0, 6.92913674516568)
        self._assert_rmsd(fitted, -1, 0.0)

    def _assert_rmsd(self, fitted, frame, desired):
        fitted.trajectory[frame]
        rmsd = MDAnalysis.analysis.align.rmsd(self.reference.atoms.coordinates(), fitted.atoms.coordinates())
        assert_almost_equal(rmsd, desired, decimal=5,
                            err_msg="frame %d of fit does not have expected RMSD" % frame)

    @attr('issue')
    def test_alignto_checks_selections(self):
        """Testing that alignto() fails if selections do not match (Issue 143)"""
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


class TestHydrogenBondAnalysis(TestCase):
    def setUp(self):
        self.universe = u = MDAnalysis.Universe(PDB_helix)
        self.kwargs = {
            'detect_hydrogens': "distance",
            'distance': 3.0,
            'angle': 150.0,
        }
        # ideal helix with 1 proline:
        self.num_bb_hbonds = u.atoms.numberOfResidues() - u.SYSTEM.PRO.numberOfResidues() - 4

    def _run(self):
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(self.universe, 'protein', 'protein', **self.kwargs)
        h.run()
        return h

    def test_helix_backbone(self):
        h = self._run()
        assert_equal(len(h.timeseries[0]), self.num_bb_hbonds, "wrong number of backbone hydrogen bonds")
        assert_equal(h.timesteps, [0.0])

    # TODO: Expand tests because the following ones are a bit superficial
    #       because we should really run them on a trajectory

    def test_count_by_time(self):
        h = self._run()
        c = h.count_by_time()
        assert_equal(c.tolist(), [(0.0, self.num_bb_hbonds)])

    def test_count_by_type(self):
        h = self._run()
        c = h.count_by_type()
        assert_equal(c.frequency, self.num_bb_hbonds * [1.0])

    def test_count_by_type(self):
        h = self._run()
        t = h.timesteps_by_type()
        assert_equal(t.time, self.num_bb_hbonds * [0.0])

    def tearDown(self):
        del self.universe


class TestHydrogenBondAnalysisHeuristic(TestHydrogenBondAnalysis):
    def setUp(self):
        super(TestHydrogenBondAnalysisHeuristic, self).setUp()
        self.kwargs['detect_hydrogens'] = "heuristic"


class TestHydrogenBondAnalysisHeavy(TestHydrogenBondAnalysis):
    def setUp(self):
        super(TestHydrogenBondAnalysisHeavy, self).setUp()
        self.kwargs['distance_type'] = "heavy"
        self.kwargs["distance"] = 3.5


class TestHydrogenBondAnalysisHeavyFail(TestHydrogenBondAnalysisHeavy):
    def setUp(self):
        super(TestHydrogenBondAnalysisHeavyFail, self).setUp()
        self.kwargs["distance"] = 3.0
        self.num_bb_hbonds = 0  # no H-bonds with a D-A distance < 3.0 A (they start at 3.05 A)


class TestAlignmentProcessing(TestCase):
    def setUp(self):
        self.seq = FASTA
        fd, self.alnfile = tempfile.mkstemp(suffix=".aln")
        os.close(fd)
        fd, self.treefile = tempfile.mkstemp(suffix=".dnd")
        os.close(fd)

    def tearDown(self):
        for f in self.alnfile, self.treefile:
            try:
                os.unlink(f)
            except OSError:
                pass

    @attr('issue')
    def test_fasta2select_aligned(self):
        """test align.fasta2select() on aligned FASTA (Issue 112)"""
        from MDAnalysis.analysis.align import fasta2select

        sel = fasta2select(self.seq, is_aligned=True)
        # length of the output strings, not residues or anything real...
        assert_equal(len(sel['reference']), 30623, err_msg="selection string has unexpected length")
        assert_equal(len(sel['mobile']), 30623, err_msg="selection string has unexpected length")

    @attr('issue')
    @dec.skipif(executable_not_found_runtime("clustalw2"),
                msg="Test skipped because clustalw executable not found")
    def test_fasta2select_ClustalW(self):
        """test align.fasta2select() with calling ClustalW (Issue 113)"""
        # note: will not be run if clustalw is not installed
        from MDAnalysis.analysis.align import fasta2select

        sel = fasta2select(self.seq, is_aligned=False, alnfilename=self.alnfile, treefilename=self.treefile)
        # numbers computed from alignment with clustalw 2.1 on Mac OS X [orbeckst]
        # length of the output strings, not residues or anything real...
        assert_equal(len(sel['reference']), 23080, err_msg="selection string has unexpected length")
        assert_equal(len(sel['mobile']), 23090, err_msg="selection string has unexpected length")


def rlimits_missing():
    # return True if resources module not accesible (ie setting of rlimits)
    try:
        # on Unix we can manipulate our limits: http://docs.python.org/2/library/resource.html
        import resource

        soft_max_open_files, hard_max_open_files = resource.getrlimit(resource.RLIMIT_NOFILE)
    except ImportError:
        return True
    return False


class TestHoleModule(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(PDB_HOLE, XTC_HOLE)
        self.dir_name = tempfile.mkdtemp()
        try:
            # on Unix we can manipulate our limits: http://docs.python.org/2/library/resource.html
            import resource

            self.soft_max_open_files, self.hard_max_open_files = resource.getrlimit(resource.RLIMIT_NOFILE)
        except ImportError:
            pass

    @attr('slow')
    @attr('issue')
    @dec.skipif(rlimits_missing, msg="Test skipped because platform does not allow setting rlimits")
    @dec.skipif(executable_not_found_runtime("hole"), msg="Test skipped because HOLE not found")
    def test_hole_module_fd_closure(self):
        """Issue 129: ensure low level file descriptors to PDB files used by Hole program are properly closed"""
        # If Issue 129 isn't resolved, this function will produce an OSError on
        # the system, and cause many other tests to fail as well.
        #
        # Successful test takes ~10 s, failure ~2 s.
        try:
            # Hasten failure by setting "ulimit -n 64" (can't go too low because of open modules etc...)
            import resource

            resource.setrlimit(resource.RLIMIT_NOFILE, (64, self.hard_max_open_files))
        except ImportError:
            raise NotImplementedError("Test cannot be run without the resource module.")
        import errno
        from MDAnalysis.analysis.hole import HOLEtraj

        os.chdir(self.dir_name)
        try:
            # will need to have the 'hole' command available in the path
            H = HOLEtraj(self.universe, cvect=[0, 1, 0], sample=20.0)
        except OSError as err:
            if err.errno == errno.ENOENT:
                raise OSError(errno.ENOENT, "HOLE binary not found")
            raise
        finally:
            self._restore_rlimits()

        # pretty unlikely that the code will get through 2 rounds if the MDA
        # issue 129 isn't fixed, although this depends on the file descriptor
        # open limit for the machine in question
        try:
            for i in xrange(2):
                # will typically get an OSError for too many files being open after
                # about 2 seconds if issue 129 isn't resolved
                H.run()
        except OSError as err:
            if err.errno == errno.EMFILE:
                raise AssertionError("HOLEtraj does not close file descriptors (Issue 129)")
            elif err.errno == errno.ENOENT:
                raise OSError(errno.ENOENT, "HOLE binary not found")
            raise
        finally:
            # make sure to restore open file limit !!
            self._restore_rlimits()

    def _restore_rlimits(self):
        try:
            import resource

            resource.setrlimit(resource.RLIMIT_NOFILE, (self.soft_max_open_files, self.hard_max_open_files))
        except ImportError:
            pass

    def tearDown(self):
        self._restore_rlimits()
        del self.universe
        import shutil

        shutil.rmtree(self.dir_name, ignore_errors=True)


class Test_Helanal(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(GRO, XTC)
        self.selection = 'name CA'
        self.tempdir = tempfile.mkdtemp()
        os.chdir(self.tempdir)

    def tearDown(self):
        del self.universe
        del self.selection

    def test_xtc_striding(self):
        """Check for sustained resolution of Issue 188."""
        u = self.universe
        sel = self.selection
        with assert_raises(FinishTimeException):
            try:
                MDAnalysis.analysis.helanal.helanal_trajectory(u, selection=sel, finish=5)
            except IndexError:
                self.fail("IndexError consistent with Issue 188.")
