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
import sys
from unittest.mock import patch
import pytest
import numpy as np
from numpy.testing import assert_equal, assert_almost_equal, assert_allclose
import itertools
from itertools import combinations_with_replacement as comb

import MDAnalysis
from MDAnalysis.lib import distances
from MDAnalysis.lib.distances import HAS_DISTOPIA
from MDAnalysis.lib import mdamath
from MDAnalysis.tests.datafiles import PSF, DCD, TRIC


def distopia_conditional_backend():
    # functions that allow distopia acceleration need to be tested with
    # distopia backend argument but distopia is an optional dep.
    if HAS_DISTOPIA:
        return ["serial", "openmp", "distopia"]
    else:
        return ["serial", "openmp"]



class TestCheckResultArray(object):

    ref = np.zeros(1, dtype=np.float64)

    def test_check_result_array_pass(self):
        # Assert input array is returned if it has correct shape and dtype:
        res = distances._check_result_array(self.ref, self.ref.shape)
        assert res is self.ref
        # Assert correct array is returned if input is None:
        res = distances._check_result_array(None, self.ref.shape)
        assert_equal(res, self.ref)
        assert res.dtype == np.float64

    def test_check_result_array_wrong_shape(self):
        wrong_shape = (1,) + self.ref.shape
        with pytest.raises(ValueError) as err:
            res = distances._check_result_array(self.ref, wrong_shape)
            assert err.msg == ("Result array has incorrect shape, should be "
                               "{0}, got {1}.".format(self.ref.shape,
                                                      wrong_shape))

    def test_check_result_array_wrong_dtype(self):
        wrong_dtype = np.int64
        ref_wrong_dtype = self.ref.astype(wrong_dtype)
        with pytest.raises(TypeError) as err:
            res = distances._check_result_array(ref_wrong_dtype, self.ref.shape)
            assert err.msg == ("Result array must be of type numpy.float64, "
                               "got {}.".format(wrong_dtype))


@pytest.mark.parametrize('coord_dtype', (np.float32, np.float64))
def test_transform_StoR_pass(coord_dtype):
    box = np.array([10, 7, 3, 45, 60, 90], dtype=np.float32)
    s = np.array([[0.5, -0.1, 0.5]], dtype=coord_dtype)

    original_r = np.array([[ 5.75,  0.36066014, 0.75]], dtype=np.float32)

    test_r = distances.transform_StoR(s, box)

    assert_allclose(original_r, test_r)


class TestCappedDistances(object):

    npoints_1 = (1, 100)

    boxes_1 = (np.array([10, 20, 30, 90, 90, 90], dtype=np.float32),  # ortho
               np.array([10, 20, 30, 30, 45, 60], dtype=np.float32),  # tri_box
               None,  # Non Periodic
               )

    @pytest.fixture()
    def query_1(self):
        q1 = np.array([0.1, 0.1, 0.1], dtype=np.float32)
        return q1

    @pytest.fixture()
    def query_2(self):
        q2 = np.array([[0.1, 0.1, 0.1], [0.2, 0.1, 0.1]], dtype=np.float32)
        return q2

    @pytest.fixture()
    def query_1_atomgroup(self, query_1):
        q1 = query_1
        u = MDAnalysis.Universe.empty(len(q1), trajectory=True)
        u.atoms.positions = q1
        return u.atoms

    @pytest.fixture()
    def query_2_atomgroup(self, query_2):
        q2 = query_2
        u = MDAnalysis.Universe.empty(len(q2), trajectory=True)
        u.atoms.positions = q2
        return u.atoms

    method_1 = ('bruteforce', 'pkdtree', 'nsgrid')

    min_cutoff_1 = (None, 0.1)

    def test_capped_distance_noresults(self):
        point1 = np.array([0.1, 0.1, 0.1], dtype=np.float32)
        point2 = np.array([0.95, 0.1, 0.1], dtype=np.float32)

        pairs, dists = distances.capped_distance(point1,
                                                 point2, max_cutoff=0.2)

        assert_equal(len(pairs), 0)

    @pytest.mark.parametrize('query', ['query_1', 'query_2',
                             'query_1_atomgroup', 'query_2_atomgroup'])
    @pytest.mark.parametrize('npoints', npoints_1)
    @pytest.mark.parametrize('box', boxes_1)
    @pytest.mark.parametrize('method', method_1)
    @pytest.mark.parametrize('min_cutoff', min_cutoff_1)
    def test_capped_distance_checkbrute(self, npoints, box, method,
                                        min_cutoff, query, request):
        q = request.getfixturevalue(query)
        np.random.seed(90003)
        points = (np.random.uniform(low=0, high=1.0,
                  size=(npoints, 3))*(self.boxes_1[0][:3])).astype(np.float32)
        max_cutoff = 2.5
        # capped distance should be able to handle array of vectors
        # as well as single vectors.
        pairs, dist = distances.capped_distance(q, points, max_cutoff,
                                                min_cutoff=min_cutoff, box=box,
                                                method=method)

        if pairs.shape != (0, ):
            found_pairs = pairs[:, 1]
        else:
            found_pairs = list()

        if isinstance(q, np.ndarray):
            if(q.shape[0] == 3):
                q = q.reshape((1, 3))

        dists = distances.distance_array(q, points, box=box)

        if min_cutoff is None:
            min_cutoff = 0.
        indices = np.where((dists <= max_cutoff) & (dists > min_cutoff))

        assert_equal(np.sort(found_pairs, axis=0), np.sort(indices[1], axis=0))

    # for coverage
    @pytest.mark.parametrize('query', ['query_1', 'query_2',
                             'query_1_atomgroup', 'query_2_atomgroup'])
    @pytest.mark.parametrize('npoints', npoints_1)
    @pytest.mark.parametrize('box', boxes_1)
    @pytest.mark.parametrize('method', method_1)
    @pytest.mark.parametrize('min_cutoff', min_cutoff_1)
    def test_capped_distance_return(self, npoints, box, query, request,
                                    method, min_cutoff):
        q = request.getfixturevalue(query)
        np.random.seed(90003)
        points = (np.random.uniform(low=0, high=1.0,
                  size=(npoints, 3))*(self.boxes_1[0][:3])).astype(np.float32)
        max_cutoff = 0.3
        # capped distance should be able to handle array of vectors
        # as well as single vectors.
        pairs = distances.capped_distance(q, points, max_cutoff,
                                          min_cutoff=min_cutoff, box=box,
                                          method=method,
                                          return_distances=False)

        if pairs.shape != (0, ):
            found_pairs = pairs[:, 1]
        else:
            found_pairs = list()

        if isinstance(q, np.ndarray):
            if(q.shape[0] == 3):
                q = q.reshape((1, 3))

        dists = distances.distance_array(q, points, box=box)

        if min_cutoff is None:
            min_cutoff = 0.
        indices = np.where((dists <= max_cutoff) & (dists > min_cutoff))

        assert_equal(np.sort(found_pairs, axis=0),
                     np.sort(indices[1], axis=0))

    def points_or_ag_self_capped(self, npoints, atomgroup=False):
        np.random.seed(90003)
        points = (np.random.uniform(low=0, high=1.0,
                  size=(npoints, 3))*(self.boxes_1[0][:3])).astype(np.float32)
        if atomgroup:
            u = MDAnalysis.Universe.empty(points.shape[0], trajectory=True)
            u.atoms.positions = points
            return u.atoms
        else:
            return points

    @pytest.mark.parametrize('npoints', npoints_1)
    @pytest.mark.parametrize('box', boxes_1)
    @pytest.mark.parametrize('method', method_1)
    @pytest.mark.parametrize('min_cutoff', min_cutoff_1)
    @pytest.mark.parametrize('ret_dist', (False, True))
    @pytest.mark.parametrize('atomgroup', (False, True))
    def test_self_capped_distance(self, npoints, box, method, min_cutoff,
                                  ret_dist, atomgroup):
        points = self.points_or_ag_self_capped(npoints, atomgroup=atomgroup)
        max_cutoff = 0.2
        result = distances.self_capped_distance(points, max_cutoff,
                                                min_cutoff=min_cutoff, box=box,
                                                method=method,
                                                return_distances=ret_dist)
        if ret_dist:
            pairs, cdists = result
        else:
            pairs = result

        # Check we found all hits
        ref = distances.self_distance_array(points, box)
        ref_d = ref[ref < 0.2]
        if min_cutoff is not None:
            ref_d = ref_d[ref_d > min_cutoff]
        assert len(ref_d) == len(pairs)

        # Go through hit by hit and check we got the indices correct too
        ref = distances.distance_array(points, points, box)
        if ret_dist:
            for (i, j), d in zip(pairs, cdists):
                d_ref = ref[i, j]
                assert d_ref < 0.2
                if min_cutoff is not None:
                    assert d_ref > min_cutoff
                assert_almost_equal(d, d_ref, decimal=6)
        else:
            for i, j in pairs:
                d_ref = ref[i, j]
                assert d_ref < 0.2
                if min_cutoff is not None:
                    assert d_ref > min_cutoff

    @pytest.mark.parametrize('box', (None,
                                     np.array([1, 1, 1,  90, 90, 90],
                                              dtype=np.float32),
                                     np.array([1, 1, 1, 60, 75, 80],
                                              dtype=np.float32)))
    @pytest.mark.parametrize('npoints,cutoff,meth',
                             [(1, 0.02, '_bruteforce_capped_self'),
                              (1, 0.2, '_bruteforce_capped_self'),
                              (600, 0.02, '_pkdtree_capped_self'),
                              (600, 0.2, '_nsgrid_capped_self')])
    def test_method_selfselection(self, box, npoints, cutoff, meth):
        np.random.seed(90003)
        points = (np.random.uniform(low=0, high=1.0,
                  size=(npoints, 3))).astype(np.float32)
        method = distances._determine_method_self(points, cutoff, box=box)
        assert_equal(method.__name__, meth)

    @pytest.mark.parametrize('box', (None,
                                     np.array([1, 1, 1,  90, 90, 90],
                                              dtype=np.float32),
                                     np.array([1, 1, 1, 60, 75, 80],
                                              dtype=np.float32)))
    @pytest.mark.parametrize('npoints,cutoff,meth',
                             [(1, 0.02, '_bruteforce_capped'),
                              (1, 0.2, '_bruteforce_capped'),
                              (200, 0.02, '_nsgrid_capped'),
                              (200, 0.35, '_bruteforce_capped'),
                              (10000, 0.35, '_nsgrid_capped')])
    def test_method_selection(self, box, npoints, cutoff, meth):
        np.random.seed(90003)
        points = (np.random.uniform(low=0, high=1.0,
                  size=(npoints, 3)).astype(np.float32))
        method = distances._determine_method(points, points, cutoff, box=box)
        assert_equal(method.__name__, meth)


@pytest.fixture()
def ref_system():
    box = np.array([1., 1., 2., 90., 90., 90], dtype=np.float32)
    points = np.array(
        [
            [0, 0, 0], [1, 1, 2], [1, 0, 2],  # identical under PBC
            [0.5, 0.5, 1.5],
        ], dtype=np.float32)
    ref = points[0:1]
    conf = points[1:]

    return box, points, ref, conf


@pytest.fixture()
def ref_system_universe(ref_system):
    box, points, ref, conf = ref_system
    u = MDAnalysis.Universe.empty(points.shape[0], trajectory=True)
    u.atoms.positions = points
    u.trajectory.ts.dimensions = box
    return (box, u.atoms, u.select_atoms("index 0"),
            u.select_atoms("index 1 to 3"))


@pytest.mark.parametrize("backend", distopia_conditional_backend())
class TestDistanceArray(object):
    @staticmethod
    def _dist(x, ref):
        ref = np.asarray(ref, dtype=np.float32)
        r = x - ref
        return np.sqrt(np.dot(r, r))

    # test both AtomGroup and numpy array
    @pytest.mark.parametrize('pos', ['ref_system', 'ref_system_universe'])
    def test_noPBC(self, backend, ref_system, pos, request):
        _, points, reference, _ = ref_system  # reference values
        _, all, ref, _ = request.getfixturevalue(pos)

        d = distances.distance_array(ref, all, backend=backend)
        assert_almost_equal(d, np.array([[
            self._dist(points[0], reference[0]),
            self._dist(points[1], reference[0]),
            self._dist(points[2], reference[0]),
            self._dist(points[3], reference[0])]
        ]))

    # cycle through combinations of numpy array and AtomGroup
    @pytest.mark.parametrize('pos0', ['ref_system', 'ref_system_universe'])
    @pytest.mark.parametrize('pos1', ['ref_system', 'ref_system_universe'])
    def test_noPBC_mixed_combinations(self, backend, ref_system, pos0, pos1,
                                      request):
        _, points, reference, _ = ref_system  # reference values
        _, _, ref_val, _ = request.getfixturevalue(pos0)
        _, points_val, _, _ = request.getfixturevalue(pos1)
        d = distances.distance_array(ref_val, points_val,
                                     backend=backend)
        assert_almost_equal(d, np.array([[
            self._dist(points[0], reference[0]),
            self._dist(points[1], reference[0]),
            self._dist(points[2], reference[0]),
            self._dist(points[3], reference[0])]
        ]))

    # test both AtomGroup and numpy array
    @pytest.mark.parametrize('pos', ['ref_system', 'ref_system_universe'])
    def test_PBC(self, backend, ref_system, pos, request):
        box, points, _, _ = ref_system
        _, all, ref, _ = request.getfixturevalue(pos)

        d = distances.distance_array(ref, all, box=box, backend=backend)

        assert_almost_equal(d, np.array([[0., 0., 0., self._dist(points[3],
                            ref=[1, 1, 2])]]))

    # cycle through combinations of numpy array and AtomGroup
    @pytest.mark.parametrize('pos0', ['ref_system', 'ref_system_universe'])
    @pytest.mark.parametrize('pos1', ['ref_system', 'ref_system_universe'])
    def test_PBC_mixed_combinations(self, backend, ref_system, pos0, pos1,
                                    request):
        box, points, _, _ = ref_system
        _, _, ref_val, _ = request.getfixturevalue(pos0)
        _, points_val, _, _ = request.getfixturevalue(pos1)
        d = distances.distance_array(ref_val, points_val,
                                     box=box,
                                     backend=backend)
        assert_almost_equal(
            d, np.array([[0., 0., 0., self._dist(points[3], ref=[1, 1, 2])]]))

    def test_PBC2(self, backend):
        a = np.array([7.90146923, -13.72858524, 3.75326586], dtype=np.float32)
        b = np.array([-1.36250901, 13.45423985, -0.36317623], dtype=np.float32)
        box = np.array([5.5457325, 5.5457325, 5.5457325, 90., 90., 90.], dtype=np.float32)

        def mindist(a, b, box):
            x = a - b
            return np.linalg.norm(x - np.rint(x / box) * box)

        ref = mindist(a, b, box[:3])
        val = distances.distance_array(a, b, box=box, backend=backend)[0, 0]

        assert_almost_equal(val, ref, decimal=6,
                            err_msg="Issue 151 not correct (PBC in distance array)")

def test_distance_array_overflow_exception():
    class FakeArray(np.ndarray):
        shape = (4294967296, 3)  # upper limit is sqrt(UINT64_MAX)
        ndim = 2
    dummy_array = FakeArray([1, 2, 3])
    box = np.array([100, 100, 100, 90., 90., 90.], dtype=np.float32)
    with pytest.raises(ValueError, match="Size of resulting array"):
        distances.distance_array.__wrapped__(dummy_array, dummy_array, box=box)

def test_self_distance_array_overflow_exception():
    class FakeArray(np.ndarray):
        shape = (6074001001, 3)  # solution of x**2 -x = 2*UINT64_MAX
        ndim = 2
    dummy_array = FakeArray([1, 2, 3])
    box = np.array([100, 100, 100, 90., 90., 90.], dtype=np.float32)
    with pytest.raises(ValueError, match="Size of resulting array"):
        distances.self_distance_array.__wrapped__(dummy_array, box=box)


@pytest.fixture()
def DCD_Universe():
    universe = MDAnalysis.Universe(PSF, DCD)
    return universe


# second independent universe required for
# TestDistanceArrayDCD_TRIC.test_atomgroup_simple
@pytest.fixture()
def DCD_Universe2():
    universe = MDAnalysis.Universe(PSF, DCD)
    return universe


@pytest.fixture()
def Triclinic_Universe():
    universe = MDAnalysis.Universe(TRIC)
    return universe

@pytest.mark.parametrize('backend', distopia_conditional_backend())
class TestDistanceArrayDCD_TRIC(object):
    # reasonable precision so that tests succeed on 32 and 64 bit machines
    # (the reference values were obtained on 64 bit)
    # Example:
    #   Items are not equal: wrong maximum distance value
    #   ACTUAL: 52.470254967456412
    #   DESIRED: 52.470257062419059
    prec = 5

    def test_simple(self, DCD_Universe, backend):
        U = DCD_Universe
        trajectory = U.trajectory
        trajectory.rewind()
        x0 = U.atoms.positions
        trajectory[10]
        x1 = U.atoms.positions
        d = distances.distance_array(x0, x1, backend=backend)
        assert_equal(d.shape, (3341, 3341), "wrong shape (should be"
                     "(Natoms,Natoms))")
        assert_almost_equal(d.min(), 0.11981228170520701, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 53.572192429459619, self.prec,
                            err_msg="wrong maximum distance value")

    def test_outarray(self, DCD_Universe, backend):
        U = DCD_Universe
        trajectory = DCD_Universe.trajectory
        trajectory.rewind()
        x0 = U.atoms.positions
        trajectory[10]
        x1 = U.atoms.positions
        natoms = len(U.atoms)
        d = np.zeros((natoms, natoms), np.float64)
        distances.distance_array(x0, x1, result=d, backend=backend)
        print("AFTER")
        print(d)
        assert_equal(d.shape, (natoms, natoms), "wrong shape, should be"
                     " (Natoms,Natoms) entries")
        assert_almost_equal(d.min(), 0.11981228170520701, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 53.572192429459619, self.prec,
                            err_msg="wrong maximum distance value")

    def test_periodic(self, DCD_Universe, backend):
        # boring with the current dcd as that has no PBC
        U = DCD_Universe
        trajectory = DCD_Universe.trajectory
        trajectory.rewind()
        x0 = U.atoms.positions
        trajectory[10]
        x1 = U.atoms.positions
        d = distances.distance_array(x0, x1, box=U.coord.dimensions,
                                     backend=backend)
        assert_equal(d.shape, (3341, 3341), "should be square matrix with"
                     " Natoms entries")
        assert_almost_equal(d.min(), 0.11981228170520701, self.prec,
                            err_msg="wrong minimum distance value with PBC")
        assert_almost_equal(d.max(), 53.572192429459619, self.prec,
                            err_msg="wrong maximum distance value with PBC")

    def test_atomgroup_simple(self, DCD_Universe, DCD_Universe2, backend):
        # need two copies as moving ts updates underlying array on atomgroup
        U1 = DCD_Universe
        U2 = DCD_Universe2
        trajectory1 = DCD_Universe.trajectory
        trajectory2 = DCD_Universe2.trajectory
        trajectory1.rewind()
        trajectory2.rewind()
        x0 = U1.select_atoms("all")
        trajectory2[10]
        x1 = U2.select_atoms("all")
        d = distances.distance_array(x0, x1, backend=backend)
        assert_equal(d.shape, (3341, 3341), "wrong shape (should be"
                     " (Natoms,Natoms))")
        assert_almost_equal(d.min(), 0.11981228170520701, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 53.572192429459619, self.prec,
                            err_msg="wrong maximum distance value")

    # check no box and ortho box types and some slices
    @pytest.mark.parametrize('box', [None, [50, 50, 50, 90, 90, 90]])
    @pytest.mark.parametrize("sel, np_slice", [("all", np.s_[:, :]),
                             ("index 0 to 8 ", np.s_[0:9, :]),
                             ("index 9", np.s_[8, :])])
    def test_atomgroup_matches_numpy(self, DCD_Universe, backend, sel,
                                     np_slice, box):
        U = DCD_Universe
        x0_ag = U.select_atoms(sel)
        x0_arr = U.atoms.positions[np_slice]
        x1_ag = U.select_atoms(sel)
        x1_arr = U.atoms.positions[np_slice]
        d_ag = distances.distance_array(x0_ag, x1_ag, box=box,
                                        backend=backend)
        d_arr = distances.distance_array(x0_arr, x1_arr, box=box,
                                         backend=backend)
        assert_allclose(d_ag, d_arr,
                        err_msg="AtomGroup and NumPy distances do not match")

    # check triclinic box and some slices
    @pytest.mark.parametrize("sel, np_slice", [("all", np.s_[:, :]),
                             ("index 0 to 8 ", np.s_[0:9, :]),
                             ("index 9", np.s_[8, :])])
    def test_atomgroup_matches_numpy_tric(self, Triclinic_Universe, backend,
                                          sel, np_slice):
        U = Triclinic_Universe
        x0_ag = U.select_atoms(sel)
        x0_arr = U.atoms.positions[np_slice]
        x1_ag = U.select_atoms(sel)
        x1_arr = U.atoms.positions[np_slice]
        d_ag = distances.distance_array(x0_ag, x1_ag, box=U.coord.dimensions,
                                        backend=backend)
        d_arr = distances.distance_array(x0_arr, x1_arr,
                                         box=U.coord.dimensions,
                                         backend=backend)
        assert_allclose(d_ag, d_arr,
                        err_msg="AtomGroup and NumPy distances do not match")


@pytest.mark.parametrize('backend', distopia_conditional_backend())
class TestSelfDistanceArrayDCD_TRIC(object):
    prec = 5

    def test_simple(self, DCD_Universe, backend):
        U = DCD_Universe
        trajectory = DCD_Universe.trajectory
        trajectory.rewind()
        x0 = U.atoms.positions
        d = distances.self_distance_array(x0, backend=backend)
        N = 3341 * (3341 - 1) / 2
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 52.4702570624190590, self.prec,
                            err_msg="wrong maximum distance value")

    def test_outarray(self, DCD_Universe, backend):
        U = DCD_Universe
        trajectory = DCD_Universe.trajectory
        trajectory.rewind()
        x0 = U.atoms.positions
        natoms = len(U.atoms)
        N = natoms * (natoms - 1) // 2
        d = np.zeros((N,), np.float64)
        distances.self_distance_array(x0, result=d, backend=backend)
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 52.4702570624190590, self.prec,
                            err_msg="wrong maximum distance value")

    def test_periodic(self, DCD_Universe, backend):
        # boring with the current dcd as that has no PBC
        U = DCD_Universe
        trajectory = DCD_Universe.trajectory
        trajectory.rewind()
        x0 = U.atoms.positions
        natoms = len(U.atoms)
        N = natoms * (natoms - 1) / 2
        d = distances.self_distance_array(x0, box=U.coord.dimensions,
                                          backend=backend)
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, self.prec,
                            err_msg="wrong minimum distance value with PBC")
        assert_almost_equal(d.max(), 52.4702570624190590, self.prec,
                            err_msg="wrong maximum distance value with PBC")

    def test_atomgroup_simple(self, DCD_Universe, backend):
        U = DCD_Universe
        trajectory = DCD_Universe.trajectory
        trajectory.rewind()
        x0 = U.select_atoms("all")
        d = distances.self_distance_array(x0, backend=backend)
        N = 3341 * (3341 - 1) / 2
        assert_equal(d.shape, (N,), "wrong shape (should be"
                     " (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 52.4702570624190590, self.prec,
                            err_msg="wrong maximum distance value")

    # check no box and ortho box types and some slices
    @pytest.mark.parametrize('box', [None, [50, 50, 50, 90, 90, 90]])
    @pytest.mark.parametrize("sel, np_slice", [("all", np.s_[:, :]),
                             ("index 0 to 8 ", np.s_[0:9, :]),
                             ("index 9", np.s_[8, :])])
    def test_atomgroup_matches_numpy(self, DCD_Universe, backend,
                                     sel, np_slice, box):
        U = DCD_Universe

        x0_ag = U.select_atoms(sel)
        x0_arr = U.atoms.positions[np_slice]
        d_ag = distances.self_distance_array(x0_ag, box=box,
                                             backend=backend)
        d_arr = distances.self_distance_array(x0_arr, box=box,
                                              backend=backend)
        assert_allclose(d_ag, d_arr,
                        err_msg="AtomGroup and NumPy distances do not match")

    # check triclinic box and some slices
    @pytest.mark.parametrize("sel, np_slice", [
                            ("index 0 to 8 ", np.s_[0:9, :]),
                            ("index 9", np.s_[8, :])])
    def test_atomgroup_matches_numpy_tric(self, Triclinic_Universe, backend,
                                          sel, np_slice):
        U = Triclinic_Universe
        x0_ag = U.select_atoms(sel)
        x0_arr = U.atoms.positions[np_slice]
        d_ag = distances.self_distance_array(x0_ag, box=U.coord.dimensions,
                                             backend=backend)
        d_arr = distances.self_distance_array(x0_arr, box=U.coord.dimensions,
                                              backend=backend)
        assert_allclose(d_ag, d_arr,
                        err_msg="AtomGroup and NumPy distances do not match")


class TestTriclinicDistances(object):
    """Unit tests for the Triclinic PBC functions.
    Tests:
      # transforming to and from S space (fractional coords)
      MDAnalysis.lib.distances.transform_StoR
      MDAnalysis.lib.distances.transform_RtoS
      # distance calculations with PBC
      MDAnalysis.lib.distances.self_distance_array
      MDAnalysis.lib.distances.distance_array
    """

    prec = 2

    @staticmethod
    @pytest.fixture()
    def TRIC():
        return MDAnalysis.Universe(TRIC)

    @staticmethod
    @pytest.fixture()
    def tri_vec_box(TRIC):
        return MDAnalysis.coordinates.core.triclinic_vectors(TRIC.dimensions)

    @staticmethod
    @pytest.fixture()
    def box(TRIC):
        return TRIC.dimensions

    @staticmethod
    @pytest.fixture()
    def S_mol(TRIC):
        S_mol1 = np.array([TRIC.atoms[383].position])
        S_mol2 = np.array([TRIC.atoms[390].position])

        return S_mol1, S_mol2

    @staticmethod
    @pytest.fixture()
    def S_mol_single(TRIC):
        S_mol1 = TRIC.atoms[383].position
        S_mol2 = TRIC.atoms[390].position
        return S_mol1, S_mol2

    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    @pytest.mark.parametrize('S_mol', [S_mol, S_mol_single], indirect=True)
    def test_transforms(self, S_mol, tri_vec_box, box, backend):
        # To check the cython coordinate transform, the same operation is done in numpy
        # Is a matrix multiplication of Coords x tri_vec_box = NewCoords, so can use np.dot
        S_mol1, S_mol2 = S_mol
        # Test transformation
        R_mol1 = distances.transform_StoR(S_mol1, box, backend=backend)
        R_np1 = np.dot(S_mol1, tri_vec_box)
        R_mol2 = distances.transform_StoR(S_mol2, box, backend=backend)
        R_np2 = np.dot(S_mol2, tri_vec_box)

        assert_almost_equal(R_mol1, R_np1, self.prec, err_msg="StoR transform failed for S_mol1")
        assert_almost_equal(R_mol2, R_np2, self.prec, err_msg="StoR transform failed for S_mol2")

        # Round trip test
        S_test1 = distances.transform_RtoS(R_mol1, box, backend=backend)
        S_test2 = distances.transform_RtoS(R_mol2, box, backend=backend)

        assert_almost_equal(S_test1, S_mol1, self.prec, err_msg="Round trip 1 failed in transform")
        assert_almost_equal(S_test2, S_mol2, self.prec, err_msg="Round trip 2 failed in transform")


    @pytest.mark.parametrize('backend', distopia_conditional_backend())
    def test_selfdist(self, S_mol, box, tri_vec_box, backend):
        S_mol1, S_mol2 = S_mol
        R_coords = distances.transform_StoR(S_mol1, box, backend="serial")
        # Transform functions are tested elsewhere so taken as working here
        dists = distances.self_distance_array(R_coords, box=box, backend=backend)
        # Manually calculate self_distance_array
        manual = np.zeros(len(dists), dtype=np.float64)
        distpos = 0
        for i, Ri in enumerate(R_coords):
            for Rj in R_coords[i + 1:]:
                Rij = Rj - Ri
                Rij -= round(Rij[2] / tri_vec_box[2][2]) * tri_vec_box[2]
                Rij -= round(Rij[1] / tri_vec_box[1][1]) * tri_vec_box[1]
                Rij -= round(Rij[0] / tri_vec_box[0][0]) * tri_vec_box[0]
                Rij = np.linalg.norm(Rij)  # find norm of Rij vector
                manual[distpos] = Rij  # and done, phew
                distpos += 1

        assert_almost_equal(dists, manual, self.prec,
                            err_msg="self_distance_array failed with input 1")

        # Do it again for input 2 (has wider separation in points)
        R_coords = distances.transform_StoR(S_mol2, box, backend="serial")
        # Transform functions are tested elsewhere so taken as working here
        dists = distances.self_distance_array(R_coords, box=box, backend=backend)
        # Manually calculate self_distance_array
        manual = np.zeros(len(dists), dtype=np.float64)
        distpos = 0
        for i, Ri in enumerate(R_coords):
            for Rj in R_coords[i + 1:]:
                Rij = Rj - Ri
                Rij -= round(Rij[2] / tri_vec_box[2][2]) * tri_vec_box[2]
                Rij -= round(Rij[1] / tri_vec_box[1][1]) * tri_vec_box[1]
                Rij -= round(Rij[0] / tri_vec_box[0][0]) * tri_vec_box[0]
                Rij = np.linalg.norm(Rij)  # find norm of Rij vector
                manual[distpos] = Rij  # and done, phew
                distpos += 1

        assert_almost_equal(dists, manual, self.prec,
                            err_msg="self_distance_array failed with input 2")


    @pytest.mark.parametrize('backend', distopia_conditional_backend())
    def test_distarray(self, S_mol, tri_vec_box, box, backend):
        S_mol1, S_mol2 = S_mol

        R_mol1 = distances.transform_StoR(S_mol1, box, backend="serial")
        R_mol2 = distances.transform_StoR(S_mol2, box, backend="serial")

        # Try with box
        dists = distances.distance_array(R_mol1, R_mol2, box=box, backend=backend)
        # Manually calculate distance_array
        manual = np.zeros((len(R_mol1), len(R_mol2)))
        for i, Ri in enumerate(R_mol1):
            for j, Rj in enumerate(R_mol2):
                Rij = Rj - Ri
                Rij -= round(Rij[2] / tri_vec_box[2][2]) * tri_vec_box[2]
                Rij -= round(Rij[1] / tri_vec_box[1][1]) * tri_vec_box[1]
                Rij -= round(Rij[0] / tri_vec_box[0][0]) * tri_vec_box[0]
                Rij = np.linalg.norm(Rij)  # find norm of Rij vector
                manual[i][j] = Rij

        assert_almost_equal(dists, manual, self.prec,
                            err_msg="distance_array failed with box")


    @pytest.mark.parametrize('backend', distopia_conditional_backend())
    def test_pbc_dist(self, S_mol, box, backend):
        S_mol1, S_mol2 = S_mol
        results = np.array([[37.629944]])
        dists = distances.distance_array(S_mol1, S_mol2, box=box, backend=backend)

        assert_almost_equal(dists, results, self.prec,
                            err_msg="distance_array failed to retrieve PBC distance")


    @pytest.mark.parametrize('backend', distopia_conditional_backend())
    def test_pbc_wrong_wassenaar_distance(self, backend):
        box = [2, 2, 2, 60, 60, 60]
        tri_vec_box = mdamath.triclinic_vectors(box)
        a, b, c = tri_vec_box
        point_a = a + b
        point_b = .5 * point_a
        dist = distances.distance_array(point_a, point_b, box=box, backend=backend)
        assert_almost_equal(dist[0, 0], 1)
        # check that our distance is different from the wassenaar distance as
        # expected.
        assert np.linalg.norm(point_a - point_b) != dist[0, 0]


@pytest.mark.parametrize("box",
    [
        None,
        np.array([10., 15., 20., 90., 90., 90.]), # otrho
        np.array([10., 15., 20., 70.53571, 109.48542, 70.518196]), # TRIC
    ]
)
def test_issue_3725(box):
    """
    Code from @hmacdope
    https://github.com/MDAnalysis/mdanalysis/issues/3725
    """
    random_coords = np.random.uniform(-50, 50, (1000, 3))

    self_da_serial = distances.self_distance_array(
        random_coords, box=box, backend='serial'
    )
    self_da_openmp = distances.self_distance_array(
        random_coords, box=box, backend='openmp'
    )

    np.testing.assert_allclose(self_da_serial, self_da_openmp)


def conv_dtype_if_ndarr(a, dtype):
    if isinstance(a, np.ndarray):
        return a.astype(dtype)
    else:
        return a


def convert_position_dtype_if_ndarray(a, b, c, d, dtype):
    return (conv_dtype_if_ndarr(a, dtype),
            conv_dtype_if_ndarr(b, dtype),
            conv_dtype_if_ndarr(c, dtype),
            conv_dtype_if_ndarr(d, dtype))




def test_HAS_DISTOPIA_distopia_too_old():
    # mock a version of distopia that is too old
    sys.modules.pop("distopia", None)
    sys.modules.pop("MDAnalysis.lib._distopia", None)   
    if HAS_DISTOPIA:
        with patch('distopia.__version__', '0.1.0'):
            with pytest.warns(RuntimeWarning,
                              match="distopia will NOT be used"):
                import MDAnalysis.lib._distopia
                assert not MDAnalysis.lib._distopia.HAS_DISTOPIA



class TestCythonFunctions(object):
    # Unit tests for calc_bonds calc_angles and calc_dihedrals in lib.distances
    # Tests both numerical results as well as input types as Cython will silently
    # produce nonsensical results if given wrong data types otherwise.
    prec = 5

    @staticmethod
    @pytest.fixture()
    def box():
        return np.array([10., 10., 10., 90., 90., 90.], dtype=np.float32)

    @staticmethod
    @pytest.fixture()
    def triclinic_box():
        box_vecs = np.array([[10., 0., 0.], [1., 10., 0., ], [1., 0., 10.]],
                            dtype=np.float32)
        return mdamath.triclinic_box(box_vecs[0], box_vecs[1], box_vecs[2])

    @staticmethod
    @pytest.fixture()
    def positions():
        # dummy atom data
        a = np.array([[0., 0., 0.], [0., 0., 0.], [0., 11., 0.], [1., 1., 1.]], dtype=np.float32)
        b = np.array([[0., 0., 0.], [1., 1., 1.], [0., 0., 0.], [29., -21., 99.]], dtype=np.float32)
        c = np.array([[0., 0., 0.], [2., 2., 2.], [11., 0., 0.], [1., 9., 9.]], dtype=np.float32)
        d = np.array([[0., 0., 0.], [3., 3., 3.], [11., -11., 0.], [65., -65., 65.]], dtype=np.float32)
        return a, b, c, d

    @staticmethod
    @pytest.fixture()
    def positions_atomgroups(positions):
        a, b, c, d = positions
        arrs = [a, b, c, d]
        universes = [MDAnalysis.Universe.empty(arr.shape[0],
                     trajectory=True) for arr in arrs]
        for u, a in zip(universes, arrs):
            u.atoms.positions = a
        return tuple([u.atoms for u in universes])

    @staticmethod
    @pytest.fixture()
    def wronglength():
        # has a different length to other inputs and should raise ValueError
        return np.array([[0., 0., 0.], [3., 3., 3.]],
                        dtype=np.float32)

    # coordinate shifts for single coord tests
    shifts = [((0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)),  # no shifting
              ((1, 0, 0), (0, 1, 1), (0, 0, 1), (1, 1, 0)),  # single box lengths
              ((-1, 0, 1), (0, -1, 0), (1, 0, 1), (-1, -1, -1)),  # negative single
              ((4, 3, -2), (-2, 2, 2), (-5, 2, 2), (0, 2, 2))]  # multiple boxlengths

    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    @pytest.mark.parametrize("pos", ["positions", "positions_atomgroups"])
    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_bonds(self, box, backend, dtype, pos, request):
        a, b, c, d = request.getfixturevalue(pos)
        a, b, c, d = convert_position_dtype_if_ndarray(a, b, c, d, dtype)
        dists = distances.calc_bonds(a, b, backend=backend)
        assert_equal(len(dists), 4, err_msg="calc_bonds results have wrong length")
        dists_pbc = distances.calc_bonds(a, b, box=box, backend=backend)
        #tests 0 length
        assert_almost_equal(dists[0], 0.0, self.prec, err_msg="Zero length calc_bonds fail")
        assert_almost_equal(dists[1], 1.7320508075688772, self.prec,
                            err_msg="Standard length calc_bonds fail")  # arbitrary length check
        # PBC checks, 2 without, 2 with
        assert_almost_equal(dists[2], 11.0, self.prec,
                            err_msg="PBC check #1 w/o box")  # pbc check 1, subtract single box length
        assert_almost_equal(dists_pbc[2], 1.0, self.prec,
                            err_msg="PBC check #1 with box")
        assert_almost_equal(dists[3], 104.26888318, self.prec,  # pbc check 2, subtract multiple box
                            err_msg="PBC check #2 w/o box")  # lengths in all directions
        assert_almost_equal(dists_pbc[3], 3.46410072, self.prec,
                            err_msg="PBC check #w with box")

    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_results_inplace_all_backends(self,  backend, dtype,):
        N = 10
        c0 = np.ones(3 * N, dtype=dtype).reshape(N, 3) * 2
        c1 = np.ones(3 * N, dtype=dtype).reshape(N, 3) * 3

        result = np.zeros(N, dtype=np.float64)
        distances.calc_bonds(c0, c1, result=result, backend=backend)
        expected = np.ones(N, dtype=dtype) * 3**(1/2)
        # test the result array is updated in place
        assert_almost_equal(result, expected, self.prec, err_msg="calc_bonds inplace failed")

    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_bonds_badbox(self, positions, backend):
        a, b, c, d = positions
        badbox1 = np.array([10., 10., 10.], dtype=np.float64)
        badbox2 = np.array([[10., 10.], [10., 10., ]], dtype=np.float32)

        with pytest.raises(ValueError):
            distances.calc_bonds(a, b, box=badbox1, backend=backend)

        with pytest.raises(ValueError):
            distances.calc_bonds(a, b, box=badbox2, backend=backend)

    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_bonds_badresult(self, positions, backend):
        a, b, c, d = positions
        badresult = np.zeros(len(a) - 1)  # Bad result array
        with pytest.raises(ValueError):
            distances.calc_bonds(a, b, result=badresult, backend=backend)

    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    @pytest.mark.parametrize("pos", ["positions", "positions_atomgroups"])
    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_bonds_triclinic(self, triclinic_box, backend, dtype, pos, request):
        a, b, c, d = request.getfixturevalue(pos)
        a, b, c, d = convert_position_dtype_if_ndarray(a, b, c, d, dtype)
        dists = distances.calc_bonds(a, b, box=triclinic_box, backend=backend)
        reference = np.array([0.0, 1.7320508, 1.4142136, 2.82842712])
        assert_almost_equal(dists, reference, self.prec, err_msg="calc_bonds with triclinic box failed")

    @pytest.mark.parametrize("shift", shifts)
    @pytest.mark.parametrize("periodic", [True, False])
    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_bonds_single_coords(self, shift, periodic, backend):
        box = np.array([10, 20, 30, 90., 90., 90.], dtype=np.float32)

        coords = np.array([[1, 1, 1], [3, 1, 1]], dtype=np.float32)

        shift1, shift2, _, _ = shift

        coords[0] += shift1 * box[:3]
        coords[1] += shift2 * box[:3]

        box = box if periodic else None
        result = distances.calc_bonds(coords[0], coords[1], box, backend=backend)

        reference = 2.0 if periodic else np.linalg.norm(coords[0] - coords[1])

        assert_almost_equal(result, reference, decimal=self.prec)

    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    @pytest.mark.parametrize("pos", ["positions", "positions_atomgroups"])
    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_angles(self, backend, dtype, pos, request):
        a, b, c, d = request.getfixturevalue(pos)
        a, b, c, d = convert_position_dtype_if_ndarray(a, b, c, d, dtype)
        angles = distances.calc_angles(a, b, c, backend=backend)
        # Check calculated values
        assert_equal(len(angles), 4, err_msg="calc_angles results have wrong length")
        #        assert_almost_equal(angles[0], 0.0, self.prec,
        #                           err_msg="Zero length angle calculation failed") # What should this be?
        assert_almost_equal(angles[1], np.pi, self.prec,
                            err_msg="180 degree angle calculation failed")
        assert_almost_equal(np.rad2deg(angles[2]), 90., self.prec,
                            err_msg="Ninety degree angle in calc_angles failed")
        assert_almost_equal(angles[3], 0.098174833, self.prec,
                            err_msg="Small angle failed in calc_angles")

    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_angles_bad_result(self, positions, backend):
        a, b, c, d = positions
        badresult = np.zeros(len(a) - 1)  # Bad result array
        with pytest.raises(ValueError):
            distances.calc_angles(a, b, c, result=badresult, backend=backend)

    @pytest.mark.parametrize(
        "case",
        [
            (
                np.array([[1, 1, 1], [1, 2, 1], [2, 2, 1]], dtype=np.float32),
                0.5 * np.pi,
            ),  # 90 degree angle
            (
                np.array([[1, 1, 1], [1, 2, 1], [1, 3, 1]], dtype=np.float32),
                np.pi,
            ),  # straight line / 180.
            (
                np.array([[1, 1, 1], [1, 2, 1], [2, 1, 1]], dtype=np.float32),
                0.25 * np.pi,
            ),  # 45
        ],
    )
    @pytest.mark.parametrize("shift", shifts)
    @pytest.mark.parametrize("periodic", [True, False])
    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_angles_single_coords(self, case, shift, periodic, backend):
        def manual_angle(x, y, z):
            return mdamath.angle(y - x, y - z)

        box = np.array([10, 20, 30, 90., 90., 90.], dtype=np.float32)
        (a, b, c), ref = case

        shift1, shift2, shift3, _ = shift

        a += shift1 * box[:3]
        b += shift2 * box[:3]
        c += shift3 * box[:3]

        box = box if periodic else None
        result = distances.calc_angles(a, b, c, box, backend=backend)
        reference = ref if periodic else manual_angle(a, b, c)
        assert_almost_equal(result, reference, decimal=4)

    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    @pytest.mark.parametrize("pos", ["positions", "positions_atomgroups"])
    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_dihedrals(self, backend, dtype, pos, request):
        a, b, c, d = request.getfixturevalue(pos)
        a, b, c, d = convert_position_dtype_if_ndarray(a, b, c, d, dtype)
        dihedrals = distances.calc_dihedrals(a, b, c, d, backend=backend)
        # Check calculated values
        assert_equal(len(dihedrals), 4, err_msg="calc_dihedrals results have wrong length")
        assert np.isnan(dihedrals[0]), "Zero length dihedral failed"
        assert np.isnan(dihedrals[1]), "Straight line dihedral failed"
        # 180 degree dihedral can be either pi or (-pi for distopia)
        assert_almost_equal(np.abs(dihedrals[2]), np.pi, self.prec, err_msg="180 degree dihedral failed")
        assert_almost_equal(dihedrals[3], -0.50714064, self.prec,
                            err_msg="arbitrary dihedral angle failed")

    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_dihedrals_wronglength(self, positions, wronglength, backend):
        a, b, c, d = positions
        with pytest.raises(ValueError):
            distances.calc_dihedrals(a, wronglength, c, d, backend=backend)

        with pytest.raises(ValueError):
            distances.calc_dihedrals(wronglength, b, c, d, backend=backend)

        with pytest.raises(ValueError):
            distances.calc_dihedrals(a, b, wronglength, d, backend=backend)

        with pytest.raises(ValueError):
            distances.calc_dihedrals(a, b, c, wronglength, backend=backend)

    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_dihedrals_bad_result(self, positions, backend):
        a, b, c, d = positions
        badresult = np.zeros(len(a) - 1)  # Bad result array

        with pytest.raises(ValueError):
            distances.calc_dihedrals(a, b, c, d, result=badresult, backend=backend)

    @pytest.mark.parametrize(
        "case",
        [
            (
                np.array(
                    [[1, 2, 1], [1, 1, 1], [2, 1, 1], [2, 2, 1]], dtype=np.float32
                ),
                0.0,
            ),  # 0 degree angle (cis)
            (
                np.array(
                    [[1, 2, 1], [1, 1, 1], [2, 1, 1], [2, 0, 1]], dtype=np.float32
                ),
                np.pi,
            ),  # 180 degree (trans)
            (
                np.array(
                    [[1, 2, 1], [1, 1, 1], [2, 1, 1], [2, 1, 2]], dtype=np.float32
                ),
                0.5 * np.pi,
            ),  # 90 degree
            (
                np.array(
                    [[1, 2, 1], [1, 1, 1], [2, 1, 1], [2, 1, 0]], dtype=np.float32
                ),
                0.5 * np.pi,
            ),  # other 90 degree
            (
                np.array(
                    [[1, 2, 1], [1, 1, 1], [2, 1, 1], [2, 2, 2]], dtype=np.float32
                ),
                0.25 * np.pi,
            ),  # 45 degree
            (
                np.array(
                    [[1, 2, 1], [1, 1, 1], [2, 1, 1], [2, 0, 2]], dtype=np.float32
                ),
                0.75 * np.pi,
            ),  # 135
        ],
    )
    @pytest.mark.parametrize("shift", shifts)
    @pytest.mark.parametrize("periodic", [True, False])
    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_dihedrals_single_coords(self, case, shift, periodic, backend):
        def manual_dihedral(a, b, c, d):
            return mdamath.dihedral(b - a, c - b, d - c)

        box = np.array([10., 10., 10., 90., 90., 90.], dtype=np.float32)

        (a, b, c, d), ref = case

        shift1, shift2, shift3, shift4 = shift

        a += shift1 * box[:3]
        b += shift2 * box[:3]
        c += shift3 * box[:3]
        d += shift4 * box[:3]

        box = box if periodic else None
        result = distances.calc_dihedrals(a, b, c, d, box, backend=backend)
        reference = ref if periodic else manual_dihedral(a, b, c, d)
        assert_almost_equal(abs(result), abs(reference), decimal=4)

    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_numpy_compliance_bonds(self, positions, backend):
        a, b, c, d = positions
        # Checks that the cython functions give identical results to the numpy versions
        bonds = distances.calc_bonds(a, b, backend=backend)
        bonds_numpy = np.array([mdamath.norm(y - x) for x, y in zip(a, b)])

        assert_almost_equal(
            bonds,
            bonds_numpy,
            self.prec,
            err_msg="Cython bonds didn't match numpy calculations",
        )

    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_numpy_compliance_angles(self, positions, backend):
        a, b, c, d = positions
        # Checks that the cython functions give identical results to the numpy versions
        angles = distances.calc_angles(a, b, c, backend=backend)
        vec1 = a - b
        vec2 = c - b
        angles_numpy = np.array([mdamath.angle(x, y) for x, y in zip(vec1, vec2)])
        # numpy 0 angle returns NaN rather than 0
        assert_almost_equal(
            angles[1:],
            angles_numpy[1:],
            self.prec,
            err_msg="Cython angles didn't match numpy calcuations",
        )

    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_numpy_compliance_dihedrals(self, positions, backend):
        a, b, c, d = positions
        # Checks that the cython functions give identical results to the numpy versions
        dihedrals = distances.calc_dihedrals(a, b, c, d, backend=backend)
        ab = a - b
        bc = b - c
        cd = c - d
        dihedrals_numpy = np.array([mdamath.dihedral(x, y, z) for x, y, z in zip(ab, bc, cd)])
        # 180 (and 360) degree dihedral can be either pi or (-pi for distopia)
        dihedrals[2] = np.abs(dihedrals[2]) 
        assert_almost_equal(dihedrals, dihedrals_numpy, self.prec,
                            err_msg="Cython dihedrals didn't match numpy calculations")


@pytest.mark.parametrize('backend', ['serial', 'openmp'])
class Test_apply_PBC(object):
    prec = 6

    @pytest.fixture()
    def DCD_universe_pos(self, DCD_Universe):
        U = DCD_Universe
        return U.atoms.positions

    @pytest.fixture()
    def DCD_universe_ag(self, DCD_Universe):
        return DCD_Universe.atoms

    @pytest.fixture()
    def Triclinic_universe_pos_box(self, Triclinic_Universe):
        U = Triclinic_Universe
        atoms = U.atoms.positions
        box = U.dimensions
        return atoms, box

    @pytest.fixture()
    def Triclinic_universe_pos_box(self, Triclinic_Universe):
        U = Triclinic_Universe
        atoms = U.atoms.positions
        box = U.dimensions
        return atoms, box

    @pytest.fixture()
    def Triclinic_universe_ag_box(self, Triclinic_Universe):
        U = Triclinic_Universe
        atoms = U.atoms
        box = U.dimensions
        return atoms, box

    @pytest.mark.parametrize('pos', ['DCD_universe_pos', 'DCD_universe_ag'])
    def test_ortho_PBC(self, backend, pos, request, DCD_universe_pos):
        positions = request.getfixturevalue(pos)
        box = np.array([2.5, 2.5, 3.5, 90., 90., 90.], dtype=np.float32)
        with pytest.raises(ValueError):
            cyth1 = distances.apply_PBC(positions, box[:3], backend=backend)
        cyth2 = distances.apply_PBC(positions, box, backend=backend)
        reference = (DCD_universe_pos -
                     np.floor(DCD_universe_pos / box[:3]) * box[:3])

        assert_almost_equal(cyth2, reference, self.prec,
                            err_msg="Ortho apply_PBC #2 failed comparison with np")

    @pytest.mark.parametrize('pos', ['Triclinic_universe_pos_box',
                             'Triclinic_universe_ag_box'])
    def test_tric_PBC(self, backend, pos, request):
        positions, box = request.getfixturevalue(pos)
        def numpy_PBC(coords, box):
            # need this to allow both AtomGroup and array
            if isinstance(coords, MDAnalysis.core.AtomGroup):
                coords = coords.positions
            else:
                pass
            # move to fractional coordinates
            fractional = distances.transform_RtoS(coords, box)
            # move fractional coordinates to central cell
            fractional -= np.floor(fractional)
            # move back to real coordinates
            return distances.transform_StoR(fractional, box)

        cyth1 = distances.apply_PBC(positions, box, backend=backend)

        reference = numpy_PBC(positions, box)

        assert_almost_equal(cyth1, reference, decimal=4,
                            err_msg="Triclinic apply_PBC failed comparison with np")

        box = np.array([10, 7, 3, 45, 60, 90], dtype=np.float32)
        r = np.array([5.75, 0.36066014, 0.75], dtype=np.float32)
        r_in_cell = distances.apply_PBC(r, box)

        assert_almost_equal([5.75, 7.3606596, 0.75], r_in_cell, self.prec)

    def test_coords_strictly_in_central_image_ortho(self, backend):
        box = np.array([10.1, 10.1, 10.1, 90.0, 90.0, 90.0], dtype=np.float32)
        # coordinates just below lower or exactly at the upper box boundaries:
        coords = np.array([[-1.0e-7, -1.0e-7, -1.0e-7],
                           [-1.0e-7, -1.0e-7,  box[2]],
                           [-1.0e-7,  box[1], -1.0e-7],
                           [ box[0], -1.0e-7, -1.0e-7],
                           [ box[0],  box[1], -1.0e-7],
                           [ box[0], -1.0e-7,  box[2]],
                           [-1.0e-7,  box[1],  box[2]],
                           [ box[0],  box[1],  box[2]]], dtype=np.float32)
        # Check that all test coordinates actually lie below the lower or
        # exactly at the upper box boundary:
        assert np.all((coords < 0.0) | (coords == box[:3]))
        res = distances.apply_PBC(coords, box, backend=backend)
        # Assert all result coordinates lie strictly within the primary image:
        assert np.all(res >= 0.0)
        assert np.all(res < box[:3])

    def test_coords_in_central_image_tric(self, backend):
        # Triclinic box corresponding to this box matrix:
        tbx = np.array([[10.1      ,  0.       ,  0.       ],
                        [ 1.0100002, 10.1      ,  0.       ],
                        [ 1.0100006,  1.0100021, 10.1      ]],
                       dtype=np.float32)
        box = mdamath.triclinic_box(*tbx)
        # coordinates just below lower or exactly at the upper box boundaries:
        coords = np.array([[  -1.0e-7,   -1.0e-7,   -1.0e-7],
                           [tbx[0, 0],   -1.0e-7,   -1.0e-7],
                           [   1.01  , tbx[1, 1],   -1.0e-7],
                           [   1.01  ,    1.01  , tbx[2, 2]],
                           [tbx[0, 0] + tbx[1, 0], tbx[1, 1], -1.0e-7],
                           [tbx[0, 0] + tbx[2, 0], 1.01, tbx[2, 2]],
                           [2.02, tbx[1, 1] + tbx[2, 1], tbx[2, 2]],
                           [tbx[0, 0] + tbx[1, 0] + tbx[2, 0],
                            tbx[1, 1] + tbx[2, 1], tbx[2, 2]]],
                          dtype=np.float32)
        relcoords = distances.transform_RtoS(coords, box)
        # Check that all test coordinates actually lie below the lower or
        # exactly at the upper box boundary:
        assert np.all((relcoords < 0.0) | (relcoords == 1.0))
        res = distances.apply_PBC(coords, box, backend=backend)
        relres = distances.transform_RtoS(res, box)
        # Assert all result coordinates lie strictly within the primary image:
        assert np.all(relres >= 0.0)
        assert np.all(relres < 1.0)


@pytest.mark.parametrize('backend', ['serial', 'openmp'])
class TestPeriodicAngles(object):
    """Test case for properly considering minimum image convention when calculating angles and dihedrals
    (Issue 172)
    """
    @staticmethod
    @pytest.fixture()
    def positions():
        a = np.array([[0.0, 1.0, 0.0]], dtype=np.float32)
        b = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        c = np.array([[1.0, 0.0, 0.0]], dtype=np.float32)
        d = np.array([[1.0, 0.0, 1.0]], dtype=np.float32)
        box = np.array([10.0, 10.0, 10.0], dtype=np.float32)
        return a, b, c, d, box

    prec = 5

    def test_angles(self, positions, backend):
        # Shift atom coordinates a few box lengths in random directions and see if we still get same results
        a, b, c, d, box = positions
        a2 = a + box * (-1, 0, 0)
        b2 = b + box * (1, 0, 1)
        c2 = c + box * (-2, 5, -7)

        ref = distances.calc_angles(a, b, c, backend=backend)

        box = np.append(box, [90, 90, 90])
        test1 = distances.calc_angles(a2, b, c, box=box, backend=backend)
        test2 = distances.calc_angles(a, b2, c, box=box, backend=backend)
        test3 = distances.calc_angles(a, b, c2, box=box, backend=backend)
        test4 = distances.calc_angles(a2, b2, c2, box=box, backend=backend)

        for val in [test1, test2, test3, test4]:
            assert_almost_equal(ref, val, self.prec, err_msg="Min image in angle calculation failed")

    def test_dihedrals(self, positions, backend):
        a, b, c, d, box = positions
        a2 = a + box * (-1, 0, 0)
        b2 = b + box * (1, 0, 1)
        c2 = c + box * (-2, 5, -7)
        d2 = d + box * (0, -5, 0)

        ref = distances.calc_dihedrals(a, b, c, d, backend=backend)

        box = np.append(box, [90, 90, 90])
        test1 = distances.calc_dihedrals(a2, b, c, d, box=box, backend=backend)
        test2 = distances.calc_dihedrals(a, b2, c, d, box=box, backend=backend)
        test3 = distances.calc_dihedrals(a, b, c2, d, box=box, backend=backend)
        test4 = distances.calc_dihedrals(a, b, c, d2, box=box, backend=backend)
        test5 = distances.calc_dihedrals(a2, b2, c2, d2, box=box, backend=backend)

        for val in [test1, test2, test3, test4, test5]:
            assert_almost_equal(ref, val, self.prec, err_msg="Min image in dihedral calculation failed")

class TestInputUnchanged(object):
    """Tests ensuring that the following functions in MDAnalysis.lib.distances
    do not alter their input coordinate arrays:
      * distance_array
      * self_distance_array
      * capped_distance
      * self_capped_distance
      * transform_RtoS
      * transform_StoR
      * calc_bonds
      * calc_angles
      * calc_dihedrals
      * apply_PBC
    """

    boxes = ([1.0, 1.0, 1.0, 90.0, 90.0, 90.0],  # orthorhombic
             [1.0, 1.0, 1.0, 80.0, 80.0, 80.0],  # triclinic
             None)  # no PBC

    @staticmethod
    @pytest.fixture()
    def coords():
        # input coordinates, some outside the [1, 1, 1] box:
        return [np.array([[0.1, 0.1, 0.1], [-0.9, -0.9, -0.9]], dtype=np.float32),
                np.array([[0.1, 0.1, 1.9], [-0.9, -0.9,  0.9]], dtype=np.float32),
                np.array([[0.1, 1.9, 1.9], [-0.9,  0.9,  0.9]], dtype=np.float32),
                np.array([[0.1, 1.9, 0.1], [-0.9,  0.9, -0.9]], dtype=np.float32)]

    @staticmethod
    @pytest.fixture()
    def coords_atomgroups(coords):
        universes = [MDAnalysis.Universe.empty(arr.shape[0], trajectory=True)
                     for arr in coords]
        for u, a in zip(universes, coords):
            u.atoms.positions = a
        return [u.atoms for u in universes]

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_input_unchanged_distance_array(self, coords, box, backend):
        crds = coords[:2]
        refs = [crd.copy() for crd in crds]
        res = distances.distance_array(crds[0], crds[1], box=box,
                                       backend=backend)
        assert_equal(crds, refs)

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_input_unchanged_distance_array_atomgroup(self, coords_atomgroups,
                                                      box, backend):
        crds = coords_atomgroups[:2]
        refs = [crd.positions.copy() for crd in crds]
        res = distances.distance_array(crds[0], crds[1], box=box,
                                       backend=backend)
        assert_equal([crd.positions for crd in crds], refs)

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_input_unchanged_self_distance_array(self, coords, box, backend):
        crd = coords[0]
        ref = crd.copy()
        res = distances.self_distance_array(crd, box=box, backend=backend)
        assert_equal(crd, ref)

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_input_unchanged_self_distance_array_atomgroup(self,
                                                           coords_atomgroups,
                                                           box, backend):
        crd = coords_atomgroups[0]
        ref = crd.positions.copy()
        res = distances.self_distance_array(crd, box=box, backend=backend)
        assert_equal(crd.positions, ref)

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('met', ["bruteforce", "pkdtree", "nsgrid", None])
    def test_input_unchanged_capped_distance(self, coords, box, met):
        crds = coords[:2]
        refs = [crd.copy() for crd in crds]
        res = distances.capped_distance(crds[0], crds[1], max_cutoff=0.3,
                                        box=box, method=met)
        assert_equal(crds, refs)

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('met', ["bruteforce", "pkdtree", "nsgrid", None])
    def test_input_unchanged_self_capped_distance(self, coords, box, met):
        crd = coords[0]
        ref = crd.copy()
        r_cut = 0.25
        res = distances.self_capped_distance(crd, max_cutoff=r_cut, box=box,
                                             method=met)
        assert_equal(crd, ref)

    @pytest.mark.parametrize('box', boxes[:2])
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_input_unchanged_transform_RtoS_and_StoR(self, coords, box, backend):
        crd = coords[0]
        ref = crd.copy()
        res = distances.transform_RtoS(crd, box, backend=backend)
        assert_equal(crd, ref)
        crd = res
        ref = crd.copy()
        res = distances.transform_StoR(crd, box, backend=backend)
        assert_equal(crd, ref)

    @pytest.mark.parametrize("box", boxes)
    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_input_unchanged_calc_bonds(self, coords, box, backend):
        crds = coords[:2]
        refs = [crd.copy() for crd in crds]
        res = distances.calc_bonds(crds[0], crds[1], box=box, backend=backend)
        assert_equal(crds, refs)

    @pytest.mark.parametrize("box", boxes)
    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_input_unchanged_calc_bonds_atomgroup(
        self, coords_atomgroups, box, backend
    ):
        crds = coords_atomgroups[:2]
        refs = [crd.positions.copy() for crd in crds]
        res = distances.calc_bonds(crds[0], crds[1], box=box, backend=backend)
        assert_equal([crd.positions for crd in crds], refs)

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('backend', distopia_conditional_backend())
    def test_input_unchanged_calc_angles(self, coords, box, backend):
        crds = coords[:3]
        refs = [crd.copy() for crd in crds]
        res = distances.calc_angles(crds[0], crds[1], crds[2], box=box,
                                    backend=backend)
        assert_equal(crds, refs)

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('backend', distopia_conditional_backend())
    def test_input_unchanged_calc_angles_atomgroup(self, coords_atomgroups,
                                                   box, backend):
        crds = coords_atomgroups[:3]
        refs = [crd.positions.copy() for crd in crds]
        res = distances.calc_angles(crds[0], crds[1], crds[2], box=box,
                                    backend=backend)
        assert_equal([crd.positions for crd in crds], refs)

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('backend', distopia_conditional_backend())
    def test_input_unchanged_calc_dihedrals(self, coords, box, backend):
        crds = coords
        refs = [crd.copy() for crd in crds]
        res = distances.calc_dihedrals(crds[0], crds[1], crds[2], crds[3],
                                       box=box, backend=backend)
        assert_equal(crds, refs)

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('backend', distopia_conditional_backend())
    def test_input_unchanged_calc_dihedrals_atomgroup(self, coords_atomgroups,
                                                      box, backend):
        crds = coords_atomgroups
        refs = [crd.positions.copy() for crd in crds]
        res = distances.calc_dihedrals(crds[0], crds[1], crds[2], crds[3],
                                       box=box, backend=backend)
        assert_equal([crd.positions for crd in crds], refs)

    @pytest.mark.parametrize('box', boxes[:2])
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_input_unchanged_apply_PBC(self, coords, box, backend):
        crd = coords[0]
        ref = crd.copy()
        res = distances.apply_PBC(crd, box, backend=backend)
        assert_equal(crd, ref)

    @pytest.mark.parametrize('box', boxes[:2])
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_input_unchanged_apply_PBC_atomgroup(self, coords_atomgroups, box,
                                                 backend):
        crd = coords_atomgroups[0]
        ref = crd.positions.copy()
        res = distances.apply_PBC(crd, box, backend=backend)
        assert_equal(crd.positions, ref)

class TestEmptyInputCoordinates(object):
    """Tests ensuring that the following functions in MDAnalysis.lib.distances
    do not choke on empty input coordinate arrays:
      * distance_array
      * self_distance_array
      * capped_distance
      * self_capped_distance
      * transform_RtoS
      * transform_StoR
      * calc_bonds
      * calc_angles
      * calc_dihedrals
      * apply_PBC
    """

    max_cut = 0.25  # max_cutoff parameter for *capped_distance()
    min_cut = 0.0  # optional min_cutoff parameter for *capped_distance()

    boxes = ([1.0, 1.0, 1.0, 90.0, 90.0, 90.0],  # orthorhombic
             [1.0, 1.0, 1.0, 80.0, 80.0, 80.0],  # triclinic
             None)  # no PBC

    @staticmethod
    @pytest.fixture()
    def empty_coord():
        # empty coordinate array:
        return np.empty((0, 3), dtype=np.float32)

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('backend', distopia_conditional_backend())
    def test_empty_input_distance_array(self, empty_coord, box, backend):
        res = distances.distance_array(empty_coord, empty_coord, box=box,
                                       backend=backend)
        assert_equal(res, np.empty((0, 0), dtype=np.float64))

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('backend', distopia_conditional_backend())
    def test_empty_input_self_distance_array(self, empty_coord, box, backend):
        res = distances.self_distance_array(empty_coord, box=box,
                                            backend=backend)
        assert_equal(res, np.empty((0,), dtype=np.float64))

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('min_cut', [min_cut, None])
    @pytest.mark.parametrize('ret_dist', [False, True])
    @pytest.mark.parametrize('met', ["bruteforce", "pkdtree", "nsgrid", None])
    def test_empty_input_capped_distance(self, empty_coord, min_cut, box, met,
                                         ret_dist):
        res = distances.capped_distance(empty_coord, empty_coord,
                                        max_cutoff=self.max_cut,
                                        min_cutoff=min_cut, box=box, method=met,
                                        return_distances=ret_dist)
        if ret_dist:
            assert_equal(res[0], np.empty((0, 2), dtype=np.int64))
            assert_equal(res[1], np.empty((0,), dtype=np.float64))
        else:
            assert_equal(res, np.empty((0, 2), dtype=np.int64))

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('min_cut', [min_cut, None])
    @pytest.mark.parametrize('ret_dist', [False, True])
    @pytest.mark.parametrize('met', ["bruteforce", "pkdtree", "nsgrid", None])
    def test_empty_input_self_capped_distance(self, empty_coord, min_cut, box,
                                              met, ret_dist):
        res = distances.self_capped_distance(empty_coord,
                                             max_cutoff=self.max_cut,
                                             min_cutoff=min_cut, box=box,
                                             method=met, return_distances=ret_dist)
        if ret_dist:
            assert_equal(res[0], np.empty((0, 2), dtype=np.int64))
            assert_equal(res[1], np.empty((0,), dtype=np.float64))
        else:
            assert_equal(res, np.empty((0, 2), dtype=np.int64))

    @pytest.mark.parametrize('box', boxes[:2])
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_empty_input_transform_RtoS(self, empty_coord, box, backend):
        res = distances.transform_RtoS(empty_coord, box, backend=backend)
        assert_equal(res, empty_coord)

    @pytest.mark.parametrize('box', boxes[:2])
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_empty_input_transform_StoR(self, empty_coord, box, backend):
        res = distances.transform_StoR(empty_coord, box, backend=backend)
        assert_equal(res, empty_coord)

    @pytest.mark.parametrize("box", boxes)
    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_empty_input_calc_bonds(self, empty_coord, box, backend):
        res = distances.calc_bonds(empty_coord, empty_coord, box=box,
                                   backend=backend)
        assert_equal(res, np.empty((0,), dtype=np.float64))

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('backend', distopia_conditional_backend())
    def test_empty_input_calc_angles(self, empty_coord, box, backend):
        res = distances.calc_angles(empty_coord, empty_coord, empty_coord,
                                    box=box, backend=backend)
        assert_equal(res, np.empty((0,), dtype=np.float64))

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('backend', distopia_conditional_backend())
    def test_empty_input_calc_dihedrals(self, empty_coord, box, backend):
        res = distances.calc_dihedrals(empty_coord, empty_coord, empty_coord,
                                       empty_coord, box=box, backend=backend)
        assert_equal(res, np.empty((0,), dtype=np.float64))

    @pytest.mark.parametrize('box', boxes[:2])
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_empty_input_apply_PBC(self, empty_coord, box, backend):
        res = distances.apply_PBC(empty_coord, box, backend=backend)
        assert_equal(res, empty_coord)


class TestOutputTypes(object):
    """Tests ensuring that the following functions in MDAnalysis.lib.distances
    return results of the types stated in the docs:
      * distance_array:
        - numpy.ndarray (shape=(n, m), dtype=numpy.float64)
      * self_distance_array:
        - numpy.ndarray (shape=(n*(n-1)//2,), dtype=numpy.float64)
      * capped_distance:
        - numpy.ndarray (shape=(n, 2), dtype=numpy.int64)
        - numpy.ndarray (shape=(n,), dtype=numpy.float64) (optional)
      * self_capped_distance:
        - numpy.ndarray (shape=(n, 2), dtype=numpy.int64)
        - numpy.ndarray (shape=(n,), dtype=numpy.float64)
      * transform_RtoS:
        - numpy.ndarray (shape=input.shape, dtype=numpy.float32)
      * transform_StoR:
        - numpy.ndarray (shape=input.shape, dtype=numpy.float32)
      * calc_bonds:
        - numpy.ndarray (shape=(n,), dtype=numpy.float64) for at least one
          shape (n,3) input, or numpy.float64 if all inputs are of shape (3,)
      * calc_angles:
        - numpy.ndarray (shape=(n,), dtype=numpy.float64) for at least one
          shape (n,3) input, or numpy.float64 if all inputs are of shape (3,)
      * calc_dihedrals:
        - numpy.ndarray (shape=(n,), dtype=numpy.float64) for at least one
          shape (n,3) input, or numpy.float64 for if all inputs are of
          shape (3,)
      * apply_PBC:
        - numpy.ndarray (shape=input.shape, dtype=numpy.float32)
    """
    max_cut = 0.25  # max_cutoff parameter for *capped_distance()
    min_cut = 0.0  # optional min_cutoff parameter for *capped_distance()

    boxes = ([1.0, 1.0, 1.0, 90.0, 90.0, 90.0],  # orthorhombic
             [1.0, 1.0, 1.0, 80.0, 80.0, 80.0],  # triclinic
             None)  # no PBC

    coords = [np.empty((0, 3), dtype=np.float32),  # empty coord array
              np.array([[0.1, 0.1, 0.1]], dtype=np.float32),  # coord array
              np.array([0.1, 0.1, 0.1], dtype=np.float32),  # single coord
              np.array([[-1.1, -1.1, -1.1]], dtype=np.float32)]  # outside box

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('incoords', list(comb(coords, 2)))
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_output_type_distance_array(self, incoords, box, backend):
        res = distances.distance_array(*incoords, box=box, backend=backend)
        assert type(res) == np.ndarray
        assert res.shape == (incoords[0].shape[0] % 2, incoords[1].shape[0] % 2)
        assert res.dtype.type == np.float64

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('incoords', coords)
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_output_type_self_distance_array(self, incoords, box, backend):
        res = distances.self_distance_array(incoords, box=box, backend=backend)
        assert type(res) == np.ndarray
        assert res.shape == (0,)
        assert res.dtype.type == np.float64

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('min_cut', [min_cut, None])
    @pytest.mark.parametrize('ret_dist', [False, True])
    @pytest.mark.parametrize('incoords', list(comb(coords, 2)))
    @pytest.mark.parametrize('met', ["bruteforce", "pkdtree", "nsgrid", None])
    def test_output_type_capped_distance(self, incoords, min_cut, box, met,
                                         ret_dist):
        res = distances.capped_distance(*incoords, max_cutoff=self.max_cut,
                                        min_cutoff=min_cut, box=box, method=met,
                                        return_distances=ret_dist)
        if ret_dist:
            pairs, dist = res
        else:
            pairs = res
        assert type(pairs) == np.ndarray
        assert pairs.dtype.type == np.intp
        assert pairs.ndim == 2
        assert pairs.shape[1] == 2
        if ret_dist:
            assert type(dist) == np.ndarray
            assert dist.dtype.type == np.float64
            assert dist.shape == (pairs.shape[0],)

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('min_cut', [min_cut, None])
    @pytest.mark.parametrize('ret_dist', [False, True])
    @pytest.mark.parametrize('incoords', coords)
    @pytest.mark.parametrize('met', ["bruteforce", "pkdtree", "nsgrid", None])
    def test_output_type_self_capped_distance(self, incoords, min_cut, box,
                                              met, ret_dist):
        res = distances.self_capped_distance(incoords,
                                                     max_cutoff=self.max_cut,
                                                     min_cutoff=min_cut,
                                                     box=box, method=met,
                                                     return_distances=ret_dist)
        if ret_dist:
            pairs, dist = res
        else:
            pairs = res
        assert type(pairs) == np.ndarray
        assert pairs.dtype.type == np.intp
        assert pairs.ndim == 2
        assert pairs.shape[1] == 2
        if ret_dist:
            assert type(dist) == np.ndarray
            assert dist.dtype.type == np.float64
            assert dist.shape == (pairs.shape[0],)

    @pytest.mark.parametrize('box', boxes[:2])
    @pytest.mark.parametrize('incoords', coords)
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_output_dtype_transform_RtoS(self, incoords, box, backend):
        res = distances.transform_RtoS(incoords, box, backend=backend)
        assert type(res) == np.ndarray
        assert res.dtype.type == np.float32
        assert res.shape == incoords.shape

    @pytest.mark.parametrize('box', boxes[:2])
    @pytest.mark.parametrize('incoords', coords)
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_output_dtype_transform_RtoS(self, incoords, box, backend):
        res = distances.transform_RtoS(incoords, box, backend=backend)
        assert type(res) == np.ndarray
        assert res.dtype.type == np.float32
        assert res.shape == incoords.shape

    @pytest.mark.parametrize("box", boxes)
    @pytest.mark.parametrize("incoords", [2 * [coords[0]]] + list(comb(coords[1:], 2)))
    @pytest.mark.parametrize("backend", distopia_conditional_backend())
    def test_output_type_calc_bonds(self, incoords, box, backend):
        res = distances.calc_bonds(*incoords, box=box, backend=backend)
        maxdim = max([crd.ndim for crd in incoords])
        if maxdim == 1:
            assert type(res) == np.float64
        else:
            assert type(res) == np.ndarray
            assert res.dtype.type == np.float64
            coord = [crd for crd in incoords if crd.ndim == maxdim][0]
            assert res.shape == (coord.shape[0],)

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('incoords',
                             [3 * [coords[0]]] + list(comb(coords[1:], 3)))
    @pytest.mark.parametrize('backend', distopia_conditional_backend())
    def test_output_type_calc_angles(self, incoords, box, backend):
        res = distances.calc_angles(*incoords, box=box, backend=backend)
        maxdim = max([crd.ndim for crd in incoords])
        if maxdim == 1:
            assert type(res) == np.float64
        else:
            assert type(res) == np.ndarray
            assert res.dtype.type == np.float64
            coord = [crd for crd in incoords if crd.ndim == maxdim][0]
            assert res.shape == (coord.shape[0],)

    @pytest.mark.parametrize('box', boxes)
    @pytest.mark.parametrize('incoords',
                             [4 * [coords[0]]] + list(comb(coords[1:], 4)))
    @pytest.mark.parametrize('backend', distopia_conditional_backend())
    def test_output_type_calc_dihedrals(self, incoords, box, backend):
        res = distances.calc_dihedrals(*incoords, box=box, backend=backend)
        maxdim = max([crd.ndim for crd in incoords])
        if maxdim == 1:
            assert type(res) == np.float64
        else:
            assert type(res) == np.ndarray
            assert res.dtype.type == np.float64
            coord = [crd for crd in incoords if crd.ndim == maxdim][0]
            assert res.shape == (coord.shape[0],)

    @pytest.mark.parametrize('box', boxes[:2])
    @pytest.mark.parametrize('incoords', coords)
    @pytest.mark.parametrize('backend', ['serial', 'openmp'])
    def test_output_type_apply_PBC(self, incoords, box, backend):
        res = distances.apply_PBC(incoords, box, backend=backend)
        assert type(res) == np.ndarray
        assert res.dtype.type == np.float32
        assert res.shape == incoords.shape


class TestDistanceBackendSelection(object):
    @staticmethod
    @pytest.fixture()
    def backend_selection_pos():
        positions = np.random.rand(10, 3)
        N = positions.shape[0]
        result = np.empty(N * (N - 1) // 2, dtype=np.float64)

        return positions, result

    @pytest.mark.parametrize('backend', [
        "serial", "Serial", "SeRiAL", "SERIAL",
        "openmp", "OpenMP", "oPENmP", "OPENMP",
    ])
    def test_case_insensitivity(self, backend, backend_selection_pos):
        positions, result = backend_selection_pos
        try:
            distances._run("calc_self_distance_array", args=(positions, result),
                           backend=backend)
        except RuntimeError:
            pytest.fail("Failed to understand backend {0}".format(backend))

    def test_wront_backend(self, backend_selection_pos):
        positions, result = backend_selection_pos
        with pytest.raises(ValueError):
            distances._run("calc_self_distance_array", args=(positions, result),
                           backend="not implemented stuff")

def test_used_openmpflag():
    assert isinstance(distances.USED_OPENMP, bool)


# test both orthognal and triclinic boxes
@pytest.mark.parametrize('box', (np.eye(3) * 10, np.array([[10, 0, 0], [2, 10, 0], [2, 2, 10]])))
# try shifts of -2 to +2 in each dimension, and all combinations of shifts
@pytest.mark.parametrize('shift', itertools.product(range(-2, 3), range(-2, 3), range(-2, 3)))
@pytest.mark.parametrize('dtype', (np.float32, np.float64))
def test_minimize_vectors(box, shift, dtype):
    # test vectors pointing in all directions
    # these currently all obey minimum convention as they're much smaller than the box
    vec = np.array(list(itertools.product(range(-1, 2), range(-1, 2), range(-1, 2))), dtype=dtype)
    box = box.astype(dtype)

    # box is 3x3 representation
    # multiply by shift, then sum xyz components then add these to the vector
    # this technically doesn't alter the vector because of periodic boundaries
    shifted_vec = (vec + (box.T * shift).sum(axis=1)).astype(dtype)

    box2 = mdamath.triclinic_box(*box).astype(dtype)

    res = distances.minimize_vectors(shifted_vec, box2)

    assert_allclose(res, vec, atol=0.00001)
    assert res.dtype == dtype
