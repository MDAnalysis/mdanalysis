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
from __future__ import print_function
import MDAnalysis
import MDAnalysis.analysis.distances
import MDAnalysis.analysis.align
import MDAnalysis.analysis.hbonds
import MDAnalysis.analysis.helanal
import MDAnalysis.analysis.rms
import MDAnalysis.analysis.waterdynamics
from MDAnalysis import SelectionError, SelectionWarning, FinishTimeException

from numpy.testing import *
import numpy as np
import nose
from nose.plugins.attrib import attr

import os
import errno
import tempfile
import itertools
import warnings

from MDAnalysisTests.datafiles import PSF, DCD, FASTA, PDB_helix, PDB_HOLE, XTC_HOLE, GRO, XTC, waterDCD, waterPSF, rmsfArray
from MDAnalysisTests import executable_not_found_runtime


class TestContactMatrix(TestCase):
    def setUp(self):
        self.coord = np.array([[1, 1, 1],
                                  [5, 5, 5],
                                  [1.1, 1.1, 1.1],
                                  [11, 11, 11],  # neighboring image with pbc
                                  [21, 21, 21]],  # non neighboring image with pbc
                                 dtype=np.float32)
        self.box = np.array([10, 10, 10], dtype=np.float32)
        self.shape = (5, 5)
        self.res_no_pbc = np.array([[1, 0, 1, 0, 0],
                                       [0, 1, 0, 0, 0],
                                       [1, 0, 1, 0, 0],
                                       [0, 0, 0, 1, 0],
                                       [0, 0, 0, 0, 1]], dtype=np.bool)
        self.res_pbc = np.array([[1, 0, 1, 1, 1],
                                    [0, 1, 0, 0, 0],
                                    [1, 0, 1, 1, 1],
                                    [1, 0, 1, 1, 1],
                                    [1, 0, 1, 1, 1]], dtype=np.bool)

    def test_np(self):
        contacts = MDAnalysis.analysis.distances.contact_matrix(
            self.coord, cutoff=1, returntype="numpy")
        assert_equal(contacts.shape, self.shape,
                     "wrong shape (should be {})".format(self.shape))
        assert_equal(contacts, self.res_no_pbc)

    def test_sparse(self):
        contacts = MDAnalysis.analysis.distances.contact_matrix(
            self.coord, cutoff=1.5, returntype="sparse")
        assert_equal(contacts.shape, self.shape,
                     "wrong shape (should be {})".format(self.shape))
        assert_equal(contacts.toarray(), self.res_no_pbc)

    def test_box_numpy(self):
        contacts = MDAnalysis.analysis.distances.contact_matrix(
            self.coord, box=self.box, cutoff=1)
        assert_equal(contacts.shape, self.shape,
                     "wrong shape (should be {})".format(self.shape))
        assert_equal(contacts, self.res_pbc)

    def test_box_sparse(self):
        contacts = MDAnalysis.analysis.distances.contact_matrix(
            self.coord, box=self.box, cutoff=1, returntype='sparse')
        assert_equal(contacts.shape, self.shape,
                     "wrong shape (should be {})".format(self.shape))
        assert_equal(contacts.toarray(), self.res_pbc)


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
        bb = self.universe.select_atoms('backbone')
        A = bb.coordinates(copy=True)  # coordinates of first frame (copy=True just in case)
        self.universe.trajectory[-1]  # forward to last frame
        B = bb.coordinates()  # coordinates of last frame
        rmsd = MDAnalysis.analysis.rms.rmsd(A, B)
        assert_almost_equal(MDAnalysis.analysis.rms.rmsd(A, A), 0.0, 5,
                            err_msg="error: rmsd(X,X) should be 0")
        # rmsd(A,B) = rmsd(B,A)  should be exact but spurious failures in the 9th decimal have been
        # observed (see Issue 57 comment #1) so we relax the test to 6 decimals.
        assert_almost_equal(MDAnalysis.analysis.rms.rmsd(B, A), rmsd, 6,
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
        rmsd = MDAnalysis.analysis.rms.rmsd(self.reference.atoms.coordinates(), fitted.atoms.coordinates())
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


class TestRMSF(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(GRO, XTC)
        fd, self.outfile = tempfile.mkstemp(suffix=".xtc")  # output is always same as input (=XTC)
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass

        del self.universe

    def test_rmsf(self):
        rmsfs = MDAnalysis.analysis.rms.RMSF(self.universe.select_atoms('name CA'))
        rmsfs.run(quiet=True)
        test_rmsfs = np.load(rmsfArray)

        assert_almost_equal(rmsfs.rmsf, test_rmsfs, 5,
                            err_msg="error: rmsf profile should match test " +
                            "values")

    def test_rmsf_single_frame(self):
        rmsfs = MDAnalysis.analysis.rms.RMSF(self.universe.select_atoms('name CA'))
        rmsfs.run(start=5, stop=6, quiet=True)

        assert_almost_equal(rmsfs.rmsf, 0, 5,
                            err_msg="error: rmsfs should all be zero")

    def test_rmsf_identical_frames(self):
        # write a dummy trajectory of all the same frame
        with MDAnalysis.Writer(self.outfile, self.universe.atoms.n_atoms) as W:
            for i in xrange(self.universe.trajectory.n_frames):
                W.write(self.universe)

        self.universe = MDAnalysis.Universe(GRO, self.outfile)
        rmsfs = MDAnalysis.analysis.rms.RMSF(self.universe.select_atoms('name CA'))
        rmsfs.run(quiet=True)

        assert_almost_equal(rmsfs.rmsf, 0, 5,
                            err_msg="error: rmsfs should all be 0")


class TestHydrogenBondAnalysis(TestCase):
    def setUp(self):
        self.universe = u = MDAnalysis.Universe(PDB_helix)
        self.kwargs = {
            'selection1': 'protein',
            'selection2': 'protein',
            'detect_hydrogens': "distance",
            'distance': 3.0,
            'angle': 150.0,
        }
        # ideal helix with 1 proline:
        self.values = {
            'num_bb_hbonds':  u.atoms.n_residues - u.SYSTEM.PRO.n_residues - 4,
            'donor_resid': np.array([5,  6,  8,  9, 10, 11, 12, 13]),
            'acceptor_resnm': np.array(['ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'PRO', 'ALA', 'ALA']),
            }

    def _run(self, **kwargs):
        kw = self.kwargs.copy()
        kw.update(kwargs)
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(self.universe, **kw)
        h.run(quiet=True)
        return h

    def test_helix_backbone(self):
        h = self._run()
        assert_equal(len(h.timeseries[0]),
                     self.values['num_bb_hbonds'], "wrong number of backbone hydrogen bonds")
        assert_equal(h.timesteps, [0.0])

    def test_generate_table(self):
        h = self._run()
        h.generate_table()
        assert_equal(len(h.table),
                     self.values['num_bb_hbonds'], "wrong number of backbone hydrogen bonds in table")
        assert_(isinstance(h.table, np.core.records.recarray))
        assert_array_equal(h.table.donor_resid, self.values['donor_resid'])
        assert_array_equal(h.table.acceptor_resnm, self.values['acceptor_resnm'])

    # TODO: Expand tests because the following ones are a bit superficial
    #       because we should really run them on a trajectory

    def test_count_by_time(self):
        h = self._run()
        c = h.count_by_time()
        assert_equal(c.tolist(), [(0.0, self.values['num_bb_hbonds'])])

    def test_count_by_type(self):
        h = self._run()
        c = h.count_by_type()
        assert_equal(c.frequency, self.values['num_bb_hbonds'] * [1.0])

    def test_count_by_type(self):
        h = self._run()
        t = h.timesteps_by_type()
        assert_equal(t.time, self.values['num_bb_hbonds'] * [0.0])

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
        self.values['num_bb_hbonds'] = 0  # no H-bonds with a D-A distance < 3.0 A (they start at 3.05 A)
        self.values['donor_resid'] = np.array([])
        self.values['acceptor_resnm'] = np.array([], dtype="|S3")


class TestHydrogenBondAnalysisChecking(object):
    def _setUp(self):
        self.universe = u = MDAnalysis.Universe(PDB_helix)
        self.kwargs = {
            'selection1': 'protein',
            'selection2': 'protein',
            'detect_hydrogens': "distance",
            'distance': 3.0,
            'angle': 150.0,
        }

    def _tearDown(self):
        del self.universe

    def _run(self, **kwargs):
        kw = self.kwargs.copy()
        kw.update(kwargs)
        with warnings.catch_warnings():
            # ignore SelectionWarning
            warnings.simplefilter("ignore")
            h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(self.universe, **kw)
            h.run(quiet=True)
        return h

    def test_check_static_selections(self):
        self._setUp()  # manually set up (because with yield cannot use TestCase)
        try:
            def run_HBA(s1, s2, s1type):
                """test that HydrogenBondAnalysis() raises SelectionError for missing donors/acceptors"""
                # no donors/acceptors; only raises error if no updates
                return self._run(selection1=s1, selection2=s2,
                                 update_selection1=False, update_selection2=False,
                                 selection1_type=s1type,
                                 )
            protein = "protein"
            nothing = "resname ALA and not backbone"
            for s1, s2, s1type in itertools.product((protein, nothing),
                                                    (protein, nothing),
                                                    ("donor", "acceptor", "both")):
                if s1 == s2 == protein:
                    def runOK():
                        """test that HydrogenBondAnalysis() works for protein/protein"""
                        try:
                            h = run_HBA(s1, s2, s1type)
                        except:
                            raise AssertionError("HydrogenBondAnalysis protein/protein failed")
                        else:
                            return True
                    yield runOK
                else:
                    yield assert_raises, SelectionError, run_HBA, s1, s2, s1type
        finally:
            self._tearDown() # manually tear down (because with yield cannot use TestCase)


    def test_run_empty_selections(self):
        self._setUp()  # manually set up (because with yield cannot use TestCase)
        try:
            def run_HBA(s1, s2, s1type):
                # no donors/acceptors; should not raise error because updates=True
                return self._run(selection1=s1, selection2=s2,
                                 update_selection1=True, update_selection2=True,
                                 selection1_type=s1type,
                                 )
            protein = "protein"
            nothing = "resname ALA and not backbone"
            for s1, s2, s1type in itertools.product((protein, nothing),
                                                    (protein, nothing),
                                                    ("donor", "acceptor", "both")):
                def run_HBA_dynamic_selections(*args):
                    try:
                        h = run_HBA(*args)
                    except:
                        raise AssertionError("HydrogenBondAnalysis with update=True failed")
                    else:
                        return True
                yield run_HBA_dynamic_selections, s1, s2, s1type
        finally:
            self._tearDown() # manually tear down (because with yield cannot use TestCase)


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
        u.trajectory[1]

        assert_raises(FinishTimeException, MDAnalysis.analysis.helanal.helanal_trajectory,
                      u, selection=sel, finish=5)
        #with assert_raises(FinishTimeException):
        #    try:
        #        MDAnalysis.analysis.helanal.helanal_trajectory(u, selection=sel, finish=5)
         #   except IndexError:
         #       self.fail("IndexError consistent with Issue 188.")


class TestWaterdynamics(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(waterPSF, waterDCD)
        self.selection1 = "byres name OH2"
        self.selection2 = self.selection1

    def test_HydrogenBondLifetimes(self):
        hbl = MDAnalysis.analysis.waterdynamics.HydrogenBondLifetimes(self.universe, self.selection1, self.selection2, 0, 5, 3)
        hbl.run(quiet=True)
        assert_equal(round(hbl.timeseries[2][1],5), 0.75)

    def test_WaterOrientationalRelaxation(self):
        wor = MDAnalysis.analysis.waterdynamics.WaterOrientationalRelaxation(self.universe, self.selection1, 0, 5, 2)
        wor.run(quiet=True)
        assert_equal(round(wor.timeseries[1][2],5), 0.35887)

    def test_AngularDistribution(self):
        ad = MDAnalysis.analysis.waterdynamics.AngularDistribution(self.universe,self.selection1,40)
        ad.run(quiet=True)
        assert_equal(str(ad.graph[0][39]), str("0.951172947884 0.48313682125") )

    def test_MeanSquareDisplacement(self):
        msd = MDAnalysis.analysis.waterdynamics.MeanSquareDisplacement(self.universe, self.selection1, 0, 10, 2)
        msd.run(quiet=True)
        assert_equal(round(msd.timeseries[1],5), 0.03984)

    def test_SurvivalProbability(self):
        sp = MDAnalysis.analysis.waterdynamics.SurvivalProbability(self.universe, self.selection1, 0, 6, 3)
        sp.run(quiet=True)
        assert_equal(round(sp.timeseries[1],5), 1.0)
