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

from six.moves import range

import MDAnalysis
import MDAnalysis.analysis.hole
from MDAnalysis.analysis.hole import HOLEtraj, HOLE

from numpy.testing import (TestCase, dec, 
                           assert_equal, assert_almost_equal, 
                           assert_array_equal,
                           assert_array_almost_equal)
import numpy as np
import nose
from nose.plugins.attrib import attr

import errno

from MDAnalysisTests.datafiles import PDB_HOLE, MULTIPDB_HOLE
from MDAnalysisTests import executable_not_found, tempdir

def rlimits_missing():
    # return True if resources module not accesible (ie setting of rlimits)
    try:
        # on Unix we can manipulate our limits: http://docs.python.org/2/library/resource.html
        import resource

        soft_max_open_files, hard_max_open_files = resource.getrlimit(resource.RLIMIT_NOFILE)
    except ImportError:
        return True
    return False

@attr('slow')
@dec.skipif(executable_not_found("hole"), msg="Test skipped because HOLE not found")
def test_HOLE(filename=PDB_HOLE):
    with tempdir.in_tempdir():
        H = HOLE(filename, raseed=31415)
        H.run()
        H.collect()
    profiles = H.profiles.values()
    assert_equal(len(profiles), 1, 
                 err_msg="HOLE.profile should contain exactly 1 profile")

    p = profiles[0]
    
    assert_equal(len(p), 425, 
                 err_msg="wrong number of points in HOLE profile")
    assert_almost_equal(p.rxncoord.mean(), -1.41225,
                        err_msg="wrong mean HOLE rxncoord")
    assert_almost_equal(p.radius.min(), 1.19707,
                        err_msg="wrong min HOLE radius")


class TestHOLEtraj(TestCase):
    filename = MULTIPDB_HOLE
    start = 5
    stop = 7

    # HOLE is so slow so we only run it once and keep it in
    # the class; note that you may not change universe.trajectory
    # (eg iteration) because this is not safe in parallel
    @classmethod
    def setUpClass(cls):
        cls.universe = MDAnalysis.Universe(cls.filename)
        if not executable_not_found("hole"):
            with tempdir.in_tempdir():
                H = HOLEtraj(cls.universe, start=cls.start, 
                             stop=cls.stop, raseed=31415)
                H.run()
            cls.H = H
        else:
            cls.H = None

        cls.frames = [ts.frame 
                      for ts in cls.universe.trajectory[cls.start:cls.stop]]

    @classmethod
    def tearDownClass(cls):
        del cls.H
        del cls.universe

    # This is VERY slow on 11 frames so we just take 2
    @attr('slow')
    @dec.skipif(executable_not_found("hole"), msg="Test skipped because HOLE not found")
    def test_HOLEtraj(self):
        assert_array_equal(sorted(self.H.profiles.keys()), self.frames,
                           err_msg="H.profiles.keys() should contain the frame numbers")

        data = np.transpose([(len(p), p.rxncoord.mean(), p.radius.min()) 
                             for p in self.H.profiles.values()])

        assert_array_equal(data[0], [401, 399], 
                           err_msg="incorrect profile lengths")
        assert_array_almost_equal(data[1], [1.98767,  0.0878],
                                  err_msg="wrong mean HOLE rxncoord")
        assert_array_almost_equal(data[2], [1.19819, 1.29628],
                                  err_msg="wrong minimum radius")

    @attr('slow')
    @dec.skipif(executable_not_found("hole"), msg="Test skipped because HOLE not found")
    def test_min_radius(self):
        assert_array_almost_equal(self.H.min_radius(),
                                  np.array([[ 5.     ,  1.19819],
                                            [ 6.     ,  1.29628]]),
                                  err_msg="min_radius() array not correct")                                  


class TestHoleModule(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(MULTIPDB_HOLE)
        try:
            # on Unix we can manipulate our limits: http://docs.python.org/2/library/resource.html
            import resource
            self.soft_max_open_files, self.hard_max_open_files = resource.getrlimit(resource.RLIMIT_NOFILE)
        except ImportError:
            pass

    @attr('slow')
    @attr('issue')
    @dec.skipif(rlimits_missing, msg="Test skipped because platform does not allow setting rlimits")
    @dec.skipif(executable_not_found("hole"), msg="Test skipped because HOLE not found")
    def test_hole_module_fd_closure(self):
        """MDAnalysis.analysis.hole: Issue 129: ensure low level file descriptors to PDB files used by Hole program are properly closed"""
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

        with tempdir.in_tempdir():
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
                for i in range(2):
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
