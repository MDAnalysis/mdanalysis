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
import numpy as np
import pytest

from numpy.testing import (
    assert_equal,
    assert_almost_equal
)

import MDAnalysis as mda
from MDAnalysisTests.coordinates.reference import RefACHE, RefCappedAla
from MDAnalysisTests.datafiles import (PRM, TRJ, TRJ_bz2, PRMpbc, TRJpbc_bz2)


class _TRJReaderTest(object):
    @pytest.fixture(scope='class')
    def universe(self):
        return mda.Universe(self.topology_file, self.trajectory_file)

    def test_load_prm(self, universe):
        assert_equal(len(universe.atoms), self.ref_n_atoms,
                     "load Universe from PRM and TRJ")

    def test_n_atoms(self, universe):
        assert_equal(universe.trajectory.n_atoms, self.ref_n_atoms,
                     "wrong number of atoms")

    def test_n_frames(self, universe):
        assert_equal(universe.trajectory.n_frames, self.ref_n_frames,
                     "wrong number of frames in xyz")

    def test_periodic(self, universe):
        assert_equal(universe.trajectory.periodic, self.ref_periodic)

    def test_amber_proteinselection(self, universe):
        protein = universe.select_atoms('protein')
        assert_equal(protein.n_atoms, self.ref_proteinatoms,
                     "error in protein selection (HIS or termini?)")

    def test_sum_centres_of_geometry(self, universe):
        protein = universe.select_atoms('protein')
        total = np.sum([protein.center_of_geometry() for ts in
                        universe.trajectory])
        assert_almost_equal(total, self.ref_sum_centre_of_geometry, self.prec,
                            err_msg="sum of centers of geometry over the "
                                    "trajectory do not match")

    def test_initial_frame_is_0(self, universe):
        assert_equal(universe.trajectory.ts.frame, 0,
                     "initial frame is not 0 but {0}".format(
                         universe.trajectory.ts.frame))

    def test_starts_with_first_frame(self, universe):
        """Test that coordinate arrays are filled as soon as the trajectory
        has been opened."""
        assert np.any(universe.atoms.positions > 0), "Reader does not " \
                                                     "populate positions right away."

    def test_rewind(self, universe):
        trj = universe.trajectory
        trj.next()
        trj.next()  # for readers that do not support indexing
        assert_equal(trj.ts.frame, 2,
                     "failed to forward to frame 2 (frameindex 2)")
        trj.rewind()
        assert_equal(trj.ts.frame, 0, "failed to rewind to first frame")
        assert np.any(universe.atoms.positions > 0), "Reader does not " \
                                                     "populate positions after rewinding."

    def test_full_slice(self, universe):
        trj_iter = universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(universe.trajectory.n_frames))

    def test_random_access(self, universe):
        pos1 = universe.atoms[0].position
        universe.trajectory.next()
        universe.trajectory.next()
        pos3 = universe.atoms[0].position

        universe.trajectory[0]

        assert_equal(universe.atoms[0].position, pos1)

        universe.trajectory[2]

        assert_equal(universe.atoms[0].position, pos3)


class TestTRJReader(_TRJReaderTest, RefACHE):
    topology_file = PRM
    trajectory_file = TRJ
    prec = 3

    def test_read_frame_reopens(self, universe):
        # should automatically reopen
        u = universe
        u.trajectory.close()
        u.trajectory[2]
        assert u.trajectory.ts.frame == 2


class TestBzippedTRJReader(TestTRJReader):
    topology_file = PRM
    trajectory_file = TRJ_bz2
    prec = 3


class TestBzippedTRJReaderPBC(_TRJReaderTest, RefCappedAla):
    topology_file = PRMpbc
    trajectory_file = TRJpbc_bz2
    prec = 3


def test_trj_no_natoms():
    with pytest.raises(ValueError):
        mda.coordinates.TRJ.TRJReader('somefile.txt')
