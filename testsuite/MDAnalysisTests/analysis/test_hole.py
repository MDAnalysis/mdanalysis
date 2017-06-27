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
from __future__ import print_function, absolute_import

from six.moves import range

import MDAnalysis
import MDAnalysis.analysis.hole
from MDAnalysis.analysis.hole import HOLEtraj, HOLE

from numpy.testing import (TestCase, dec,
                           assert_equal, assert_almost_equal,
                           assert_array_equal,
                           assert_array_almost_equal, assert_)
import numpy as np
import matplotlib
import mpl_toolkits.mplot3d

import nose
from nose.plugins.attrib import attr

import os
import errno

from MDAnalysisTests.datafiles import PDB_HOLE, MULTIPDB_HOLE
from MDAnalysisTests import (executable_not_found, module_not_found,
                             tempdir, in_dir)

def rlimits_missing():
    # return True if resources module not accesible (ie setting of rlimits)
    try:
        # on Unix we can manipulate our limits: http://docs.python.org/2/library/resource.html
        import resource

        soft_max_open_files, hard_max_open_files = resource.getrlimit(resource.RLIMIT_NOFILE)
    except ImportError:
        return True
    return False


class TestHOLE(TestCase):
    filename = PDB_HOLE

    @dec.skipif(executable_not_found("hole"), msg="Test skipped because HOLE not found")
    def setUp(self):
        # keep tempdir around for the whole lifetime of the class
        self.tempdir = tempdir.TempDir()
        with in_dir(self.tempdir.name):
            H = HOLE(self.filename, raseed=31415)
            H.run()
            H.collect()
        self.H = H

    def tearDown(self):
        del self.H
        del self.tempdir

    @attr('slow')
    @dec.skipif(executable_not_found("hole"), msg="Test skipped because HOLE not found")
    def test_HOLE(self):
        profiles = self.H.profiles.values()
        assert_equal(len(profiles), 1,
                     err_msg="HOLE.profile should contain exactly 1 profile")

        p = profiles[0]

        assert_equal(len(p), 425,
                     err_msg="wrong number of points in HOLE profile")
        assert_almost_equal(p.rxncoord.mean(), -1.41225,
                            err_msg="wrong mean HOLE rxncoord")
        assert_almost_equal(p.radius.min(), 1.19707,
                            err_msg="wrong min HOLE radius")

    @attr('slow')
    @dec.skipif(executable_not_found("hole"), msg="Test skipped because HOLE not found")
    def test_vmd_surface(self):
        with in_dir(self.tempdir.name):
            filename = self.H.create_vmd_surface(filename="hole.vmd")
            assert_equal(len(open(filename).readlines()), 6504,
                         err_msg="HOLE VMD surface file is incomplete")

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

    @attr('slow')
    @dec.skipif(executable_not_found("hole"), msg="Test skipped because HOLE not found")
    def test_plot(self):
        ax = self.H.plot(label=True)
        assert_(isinstance(ax, matplotlib.axes.Axes),
                msg="H.plot() did not produce an Axes instance")

    @attr('slow')
    @dec.skipif(executable_not_found("hole"), msg="Test skipped because HOLE not found")
    def test_plot3D(self):
        ax = self.H.plot3D()
        assert_(isinstance(ax, mpl_toolkits.mplot3d.Axes3D),
                msg="H.plot3D() did not produce an Axes3D instance")

    @attr('slow')
    @dec.skipif(executable_not_found("hole"), msg="Test skipped because HOLE not found")
    def test_plot3D_rmax(self):
        ax = self.H.plot3D(rmax=2.5)
        assert_(isinstance(ax, mpl_toolkits.mplot3d.Axes3D),
                msg="H.plot3D(rmax=float) did not produce an Axes3D instance")


class TestHoleModule(TestCase):
    @dec.skipif(rlimits_missing, msg="Test skipped because platform does not allow setting rlimits")
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
        """test open file descriptors are closed (MDAnalysisTests.analysis.test_hole.TestHoleModule): Issue 129"""
        # If Issue 129 isn't resolved, this function will produce an OSError on
        # the system, and cause many other tests to fail as well.
        #
        # Successful test takes ~10 s, failure ~2 s.

        # Hasten failure by setting "ulimit -n 64" (can't go too low because of open modules etc...)
        import resource

        # ----- temporary hack -----
        # on Mac OS X (on Travis) we run out of open file descriptors
        # before even starting this test (see
        # https://github.com/MDAnalysis/mdanalysis/pull/901#issuecomment-231938093);
        # if this issue is solved by #363 then revert the following
        # hack:
        #
        import platform
        if platform.platform() == "Darwin":
            max_open_files = 512
        else:
            max_open_files = 64
        #
        # --------------------------

        resource.setrlimit(resource.RLIMIT_NOFILE,
                           (max_open_files, self.hard_max_open_files))

        with tempdir.in_tempdir():
            try:
                H = HOLEtraj(self.universe, cvect=[0, 1, 0], sample=20.0)
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
                raise
            finally:
                # make sure to restore open file limit !!
                self._restore_rlimits()

    def _restore_rlimits(self):
        try:
            import resource
            resource.setrlimit(resource.RLIMIT_NOFILE,
                               (self.soft_max_open_files, self.hard_max_open_files))
        except ImportError:
            pass

    def tearDown(self):
        self._restore_rlimits()
        del self.universe
