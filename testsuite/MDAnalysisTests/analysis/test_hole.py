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
from __future__ import print_function, absolute_import

from six.moves import range
import pytest

import MDAnalysis
import MDAnalysis.analysis.hole
from MDAnalysis.analysis.hole import HOLEtraj, HOLE

from numpy.testing import (
    assert_almost_equal,
    assert_array_equal,
    assert_array_almost_equal
)
import numpy as np
import matplotlib
import mpl_toolkits.mplot3d

import errno

from MDAnalysisTests.datafiles import PDB_HOLE, MULTIPDB_HOLE
from MDAnalysisTests import executable_not_found


def rlimits_missing():
    # return True if resources module not accesible (ie setting of rlimits)
    try:
        # on Unix we can manipulate our limits: http://docs.python.org/2/library/resource.html
        import resource

        soft_max_open_files, hard_max_open_files = resource.getrlimit(resource.RLIMIT_NOFILE)
    except ImportError:
        return True
    return False


@pytest.mark.skipif(executable_not_found("hole"), reason="Test skipped because HOLE not found")
class TestHOLE(object):
    @staticmethod
    @pytest.fixture()
    def H(tmpdir):
        # keep tempdir around for the whole lifetime of the class
        with tmpdir.as_cwd():
            filename = PDB_HOLE
            H = HOLE(filename, raseed=31415)
            H.run()
            H.collect()
        return H

    def test_HOLE(self, H):
        profiles_values = list(H.profiles.values())
        assert len(profiles_values) == 1, "HOLE.profile should contain exactly 1 profile"

        p = profiles_values[0]

        assert len(p) == 425, "wrong number of points in HOLE profile"
        assert_almost_equal(p.rxncoord.mean(), -1.41225, err_msg="wrong mean HOLE rxncoord")
        assert_almost_equal(p.radius.min(), 1.19707, err_msg="wrong min HOLE radius")

    def test_vmd_surface(self, tmpdir, H):
        with tmpdir.as_cwd():
            filename = H.create_vmd_surface(filename="hole.vmd")
            assert len(open(filename).readlines()) == 6504, "HOLE VMD surface file is incomplete"


@pytest.mark.skipif(executable_not_found("hole"), reason="Test skipped because HOLE not found")
class TestHOLEtraj(object):
    filename = MULTIPDB_HOLE
    start = 5
    stop = 7

    # HOLE is so slow so we only run it once and keep it in
    # the class; note that you may not change universe.trajectory
    # (eg iteration) because this is not safe in parallel
    @pytest.fixture()
    def universe(self):
        return MDAnalysis.Universe(self.filename)

    @pytest.fixture()
    def H(self, universe, tmpdir):
        with tmpdir.as_cwd():
            H = HOLEtraj(universe, raseed=31415)
            H.run(start=self.start, stop=self.stop)
        return H

    @pytest.fixture()
    def frames(self, universe):
        return [ts.frame for ts in universe.trajectory[self.start:self.stop]]

    # This is VERY slow on 11 frames so we just take 2
    def test_HOLEtraj(self, H, frames):
        assert_array_equal(sorted(H.profiles.keys()), frames,
                           err_msg="H.profiles.keys() should contain the frame numbers")

        data = np.transpose([(len(p), p.rxncoord.mean(), p.radius.min())
                             for p in H.profiles.values()])

        assert_array_equal(data[0], [401, 399], err_msg="incorrect profile lengths")
        assert_array_almost_equal(data[1], [1.98767,  0.0878], err_msg="wrong mean HOLE rxncoord")
        assert_array_almost_equal(data[2], [1.19819, 1.29628], err_msg="wrong minimum radius")

    def test_min_radius(self, H):
        assert_array_almost_equal(
            H.min_radius(),
            np.array([[5., 1.19819],
                      [6., 1.29628]]),
            err_msg="min_radius() array not correct"
        )

    def test_plot(self, H):
        ax = H.plot(label=True)
        assert isinstance(ax, matplotlib.axes.Axes), "H.plot() did not produce an Axes instance"

    def test_plot3D(self, H):
        ax = H.plot3D()
        assert isinstance(ax, mpl_toolkits.mplot3d.Axes3D), "H.plot3D() did not produce an Axes3D instance"

    def test_plot3D_rmax(self, H):
        ax = H.plot3D(rmax=2.5)
        assert isinstance(ax, mpl_toolkits.mplot3d.Axes3D), "H.plot3D(rmax=float) did not produce an Axes3D instance"


@pytest.mark.skipif(executable_not_found("hole"), reason="Test skipped because HOLE not found")
class TestHoleModule(object):
    try:
        # on Unix we can manipulate our limits: http://docs.python.org/2/library/resource.html
        import resource
        soft_max_open_files, hard_max_open_files = resource.getrlimit(resource.RLIMIT_NOFILE)
    except ImportError:
        pass

    @staticmethod
    @pytest.fixture()
    def universe():
        return MDAnalysis.Universe(MULTIPDB_HOLE)

    @pytest.mark.skipif(rlimits_missing, reason="Test skipped because platform does not allow setting rlimits")
    def test_hole_module_fd_closure(self, universe, tmpdir):
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

        with tmpdir.as_cwd():
            try:
                H = HOLEtraj(universe, cvect=[0, 1, 0], sample=20.0)
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
                    raise pytest.fail("HOLEtraj does not close file descriptors (Issue 129)")
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
