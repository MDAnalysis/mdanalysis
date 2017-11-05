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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import division, absolute_import
import MDAnalysis
import MDAnalysis.lib.distances

import numpy as np
import pytest
from numpy.testing import assert_almost_equal, assert_equal

from MDAnalysis.tests.datafiles import PSF, DCD, TRIC
from MDAnalysis.lib import mdamath


@pytest.fixture()
def ref_system():
    box = np.array([1., 1., 2.], dtype=np.float32)
    points = np.array(
        [
            [0, 0, 0], [1, 1, 2], [1, 0, 2],  # identical under PBC
            [0.5, 0.5, 1.5],
        ], dtype=np.float32)
    ref = points[0:1]
    conf = points[1:]

    return box, points, ref, conf


@pytest.mark.parametrize('backend', ['serial', 'openmp'])
class TestDistanceArray(object):
    @staticmethod
    def _dist(x, ref):
        ref = np.asarray(ref, dtype=np.float32)
        r = x - ref
        return np.sqrt(np.dot(r, r))

    def test_noPBC(self, backend, ref_system):
        box, points, ref, conf = ref_system

        d = MDAnalysis.lib.distances.distance_array(ref, points, backend=backend)

        assert_almost_equal(d, np.array([[
            self._dist(points[0], ref[0]),
            self._dist(points[1], ref[0]),
            self._dist(points[2], ref[0]),
            self._dist(points[3], ref[0])]
        ]))

    def test_PBC(self, backend, ref_system):
        box, points, ref, conf = ref_system

        d = MDAnalysis.lib.distances.distance_array(ref, points, box=box, backend=backend)

        assert_almost_equal(d, np.array([[0., 0., 0., self._dist(points[3], ref=[1, 1, 2])]]))

    def test_PBC2(self, backend):
        a = np.array([7.90146923, -13.72858524, 3.75326586], dtype=np.float32)
        b = np.array([-1.36250901, 13.45423985, -0.36317623], dtype=np.float32)
        box = np.array([5.5457325, 5.5457325, 5.5457325], dtype=np.float32)

        def mindist(a, b, box):
            x = a - b
            return np.linalg.norm(x - np.rint(x / box) * box)

        ref = mindist(a, b, box)
        val = MDAnalysis.lib.distances.distance_array(np.array([a]), np.array([b]),
                                                      box=box, backend=backend)[0, 0]

        assert_almost_equal(val, ref, decimal=6,
                            err_msg="Issue 151 not correct (PBC in distance array)")

@pytest.fixture()
def DCD_Universe():
    universe = MDAnalysis.Universe(PSF, DCD)
    trajectory = universe.trajectory

    return universe, trajectory

@pytest.mark.parametrize('backend', ['serial', 'openmp'])
class TestDistanceArrayDCD(object):
    # reasonable precision so that tests succeed on 32 and 64 bit machines
    # (the reference values were obtained on 64 bit)
    # Example:
    #   Items are not equal: wrong maximum distance value
    #   ACTUAL: 52.470254967456412
    #   DESIRED: 52.470257062419059
    prec = 5

    def test_simple(self, DCD_Universe, backend):
        U, trajectory = DCD_Universe
        trajectory.rewind()
        x0 = U.atoms.positions
        trajectory[10]
        x1 = U.atoms.positions
        d = MDAnalysis.lib.distances.distance_array(x0, x1, backend=backend)
        assert_equal(d.shape, (3341, 3341), "wrong shape (should be (Natoms,Natoms))")
        assert_almost_equal(d.min(), 0.11981228170520701, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 53.572192429459619, self.prec,
                            err_msg="wrong maximum distance value")

    def test_outarray(self, DCD_Universe, backend):
        U, trajectory = DCD_Universe
        trajectory.rewind()
        x0 = U.atoms.positions
        trajectory[10]
        x1 = U.atoms.positions
        natoms = len(U.atoms)
        d = np.zeros((natoms, natoms), np.float64)
        MDAnalysis.lib.distances.distance_array(x0, x1, result=d, backend=backend)
        assert_equal(d.shape, (natoms, natoms), "wrong shape, shoud be  (Natoms,Natoms) entries")
        assert_almost_equal(d.min(), 0.11981228170520701, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 53.572192429459619, self.prec,
                            err_msg="wrong maximum distance value")

    def test_periodic(self, DCD_Universe, backend):
        # boring with the current dcd as that has no PBC
        U, trajectory = DCD_Universe
        trajectory.rewind()
        x0 = U.atoms.positions
        trajectory[10]
        x1 = U.atoms.positions
        d = MDAnalysis.lib.distances.distance_array(x0, x1, box=U.coord.dimensions,
                                                    backend=backend)
        assert_equal(d.shape, (3341, 3341), "should be square matrix with Natoms entries")
        assert_almost_equal(d.min(), 0.11981228170520701, self.prec,
                            err_msg="wrong minimum distance value with PBC")
        assert_almost_equal(d.max(), 53.572192429459619, self.prec,
                            err_msg="wrong maximum distance value with PBC")



@pytest.mark.parametrize('backend', ['serial', 'openmp'])
class TestSelfDistanceArrayDCD(object):
    prec = 5

    def test_simple(self, DCD_Universe, backend):
        U, trajectory = DCD_Universe
        trajectory.rewind()
        x0 = U.atoms.positions
        d = MDAnalysis.lib.distances.self_distance_array(x0, backend=backend)
        N = 3341 * (3341 - 1) / 2
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 52.4702570624190590, self.prec,
                            err_msg="wrong maximum distance value")

    def test_outarray(self, DCD_Universe, backend):
        U, trajectory = DCD_Universe
        trajectory.rewind()
        x0 = U.atoms.positions
        natoms = len(U.atoms)
        N = natoms * (natoms - 1) // 2
        d = np.zeros((N,), np.float64)
        MDAnalysis.lib.distances.self_distance_array(x0, result=d, backend=backend)
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 52.4702570624190590, self.prec,
                            err_msg="wrong maximum distance value")

    def test_periodic(self, DCD_Universe, backend):
        # boring with the current dcd as that has no PBC
        U, trajectory = DCD_Universe
        trajectory.rewind()
        x0 = U.atoms.positions
        natoms = len(U.atoms)
        N = natoms * (natoms - 1) / 2
        d = MDAnalysis.lib.distances.self_distance_array(x0, box=U.coord.dimensions,
                                                         backend=backend)
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, self.prec,
                            err_msg="wrong minimum distance value with PBC")
        assert_almost_equal(d.max(), 52.4702570624190590, self.prec,
                            err_msg="wrong maximum distance value with PBC")


@pytest.mark.parametrize('backend', ['serial', 'openmp'])
class TestTriclinicDistances(object):
    """Unit tests for the Triclinic PBC functions.
    Tests:
      # transforming to and from S space (fractional coords)
      mda.lib.distances.transform_StoR
      mda.lib.distances.transform_RtoS
      # distance calculations with PBC
      mda.lib.distances.self_distance_array
      mda.lib.distances.distance_array
    """

    prec = 2

    @staticmethod
    @pytest.fixture()
    def TRIC():
        return MDAnalysis.Universe(TRIC)

    @staticmethod
    @pytest.fixture()
    def box(TRIC):
        return MDAnalysis.coordinates.core.triclinic_vectors(TRIC.dimensions)

    @staticmethod
    @pytest.fixture()
    def boxV(box):
        return MDAnalysis.coordinates.core.triclinic_box(box[0], box[1], box[2])

    @staticmethod
    @pytest.fixture()
    def S_mol(TRIC):
        S_mol1 = np.array([TRIC.atoms[383].position])
        S_mol2 = np.array([TRIC.atoms[390].position])

        return S_mol1, S_mol2

    def test_transforms(self, S_mol, box, boxV, backend):
        from MDAnalysis.lib.distances import transform_StoR, transform_RtoS
        # To check the cython coordinate transform, the same operation is done in numpy
        # Is a matrix multiplication of Coords x Box = NewCoords, so can use np.dot
        S_mol1, S_mol2 = S_mol
        # Test transformation
        R_mol1 = transform_StoR(S_mol1, box, backend=backend)
        R_np1 = np.dot(S_mol1, box)

        # Test transformation when given box in different form
        R_mol2 = transform_StoR(S_mol2, boxV, backend=backend)
        R_np2 = np.dot(S_mol2, box)

        assert_almost_equal(R_mol1, R_np1, self.prec, err_msg="StoR transform failed with box")
        assert_almost_equal(R_mol2, R_np2, self.prec, err_msg="StoR transform failed with boxV")

        # Round trip test
        # boxV here althought initial transform with box
        S_test1 = transform_RtoS(R_mol1, boxV, backend=backend)
        # and vice versa, should still work
        S_test2 = transform_RtoS(R_mol2, box, backend=backend)

        assert_almost_equal(S_test1, S_mol1, self.prec, err_msg="Round trip failed in transform")
        assert_almost_equal(S_test2, S_mol2, self.prec, err_msg="Round trip failed in transform")

    def test_selfdist(self, S_mol, box, boxV, backend):
        from MDAnalysis.lib.distances import self_distance_array
        from MDAnalysis.lib.distances import transform_StoR

        S_mol1, S_mol2 = S_mol

        R_coords = transform_StoR(S_mol1, box, backend=backend)
        # Transform functions are tested elsewhere so taken as working here
        dists = self_distance_array(R_coords, box=box, backend=backend)
        # Manually calculate self_distance_array
        manual = np.zeros(len(dists), dtype=np.float64)
        distpos = 0
        for i, Ri in enumerate(R_coords):
            for Rj in R_coords[i + 1:]:
                Rij = Rj - Ri
                Rij -= round(Rij[2] / box[2][2]) * box[2]
                Rij -= round(Rij[1] / box[1][1]) * box[1]
                Rij -= round(Rij[0] / box[0][0]) * box[0]
                Rij = np.linalg.norm(Rij)  # find norm of Rij vector
                manual[distpos] = Rij  # and done, phew
                distpos += 1

        assert_almost_equal(dists, manual, self.prec,
                            err_msg="self_distance_array failed with input 1")

        # Do it again for input 2 (has wider separation in points)
        # Also use boxV here in self_dist calculation
        R_coords = transform_StoR(S_mol2, box, backend=backend)
        # Transform functions are tested elsewhere so taken as working here
        dists = self_distance_array(R_coords, box=boxV, backend=backend)
        # Manually calculate self_distance_array
        manual = np.zeros(len(dists), dtype=np.float64)
        distpos = 0
        for i, Ri in enumerate(R_coords):
            for Rj in R_coords[i + 1:]:
                Rij = Rj - Ri
                Rij -= round(Rij[2] / box[2][2]) * box[2]
                Rij -= round(Rij[1] / box[1][1]) * box[1]
                Rij -= round(Rij[0] / box[0][0]) * box[0]
                Rij = np.linalg.norm(Rij)  # find norm of Rij vector
                manual[distpos] = Rij  # and done, phew
                distpos += 1

        assert_almost_equal(dists, manual, self.prec,
                            err_msg="self_distance_array failed with input 2")

    def test_distarray(self, S_mol, box, boxV, backend):
        from MDAnalysis.lib.distances import distance_array
        from MDAnalysis.lib.distances import transform_StoR

        S_mol1, S_mol2 = S_mol

        R_mol1 = transform_StoR(S_mol1, box, backend=backend)
        R_mol2 = transform_StoR(S_mol2, box, backend=backend)

        # Try with box
        dists = distance_array(R_mol1, R_mol2, box=box, backend=backend)
        # Manually calculate distance_array
        manual = np.zeros((len(R_mol1), len(R_mol2)))
        for i, Ri in enumerate(R_mol1):
            for j, Rj in enumerate(R_mol2):
                Rij = Rj - Ri
                Rij -= round(Rij[2] / box[2][2]) * box[2]
                Rij -= round(Rij[1] / box[1][1]) * box[1]
                Rij -= round(Rij[0] / box[0][0]) * box[0]
                Rij = np.linalg.norm(Rij)  # find norm of Rij vector
                manual[i][j] = Rij

        assert_almost_equal(dists, manual, self.prec,
                            err_msg="distance_array failed with box")

        # Now check using boxV
        dists = distance_array(R_mol1, R_mol2, box=boxV, backend=backend)
        assert_almost_equal(dists, manual, self.prec,
                            err_msg="distance_array failed with boxV")

    def test_pbc_dist(self, S_mol, boxV, backend):
        from MDAnalysis.lib.distances import distance_array
        S_mol1, S_mol2 = S_mol

        results = np.array([[37.629944]])
        dists = distance_array(S_mol1, S_mol2, box=boxV,
                               backend=backend)

        assert_almost_equal(dists, results, self.prec,
                            err_msg="distance_array failed to retrieve PBC distance")

    def test_pbc_wrong_wassenaar_distance(self, backend):
        from MDAnalysis.lib.distances import distance_array
        box = MDAnalysis.lib.mdamath.triclinic_vectors([2, 2, 2, 60, 60, 60])
        a, b, c = box
        point_a = a + b
        point_b = .5 * point_a
        dist = distance_array(point_a[np.newaxis, :], point_b[np.newaxis, :],
                              box=box, backend=backend)
        assert_almost_equal(dist[0, 0], 1)
        # check that our distance is different then the wassenaar distance as
        # expected.
        assert np.linalg.norm(point_a - point_b) != dist[0, 0]



@pytest.mark.parametrize('backend', ['serial', 'openmp'])
class TestCythonFunctions(object):
    # Unit tests for calc_bonds calc_angles and calc_dihedrals in lib.distances
    # Tests both numerical results as well as input types as Cython will silently
    # produce nonsensical results if given wrong data types otherwise.
    prec = 5

    @staticmethod
    @pytest.fixture()
    def box():
        return np.array([10., 10., 10.], dtype=np.float32)

    @staticmethod
    @pytest.fixture()
    def triclinic_box():
        return np.array([[10., 0., 0.], [1., 10., 0., ], [1., 0., 10.]], dtype=np.float32)

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
    def wrongtype():
        # declared as float64 and should raise TypeError
        return np.array([[0., 0., 0.], [3., 3., 3.], [3., 3., 3.], [3., 3., 3.]],
                        dtype=np.float64)

    @staticmethod
    @pytest.fixture()
    def wronglength():
        # has a different length to other inputs and should raise ValueError
        return np.array([[0., 0., 0.], [3., 3., 3.]],
                        dtype=np.float32)

    def test_bonds(self, positions, box, backend):
        a, b, c, d = positions
        dists = MDAnalysis.lib.distances.calc_bonds(a, b,
                                                    backend=backend)
        assert_equal(len(dists), 4, err_msg="calc_bonds results have wrong length")
        dists_pbc = MDAnalysis.lib.distances.calc_bonds(a, b, box=box,
                                                        backend=backend)
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

    #Bad input checking
    def test_bonds_wrongtype(self, positions, wrongtype, wronglength, backend):
        a, b, c, d = positions

        with pytest.raises(TypeError):
            MDAnalysis.lib.distances.calc_bonds(a, wrongtype, backend=backend)
        with pytest.raises(TypeError):
            MDAnalysis.lib.distances.calc_bonds(wrongtype, b, backend=backend)

        with pytest.raises(ValueError):
            MDAnalysis.lib.distances.calc_bonds(a, wronglength, backend=backend)

        with pytest.raises(ValueError):
            MDAnalysis.lib.distances.calc_bonds(wronglength, b, backend=
            backend)


    def test_bonds_badbox(self, positions, backend):
        a, b, c, d = positions
        badboxtype = np.array([10., 10., 10.], dtype=np.float64)
        badboxsize = np.array([[10., 10.], [10., 10., ]], dtype=np.float32)

        with pytest.raises(ValueError):
            MDAnalysis.lib.distances.calc_bonds(a,
                                                b, box=badboxsize,
                                                backend=backend)  # Bad box data

        with pytest.raises(TypeError):
            MDAnalysis.lib.distances.calc_bonds(a,
                                                b, box=badboxtype,
                                                backend=backend)  # Bad box type

    def test_bonds_badresult(self, positions, backend):
        a, b, c, d = positions
        badresult = np.zeros(len(a) - 1)
        with pytest.raises(ValueError):
            MDAnalysis.lib.distances.calc_bonds(
                a, b, result=badresult, backend=backend)  # Bad result array

    def test_bonds_triclinic(self, positions, triclinic_box, backend):
        a, b, c, d = positions
        dists = MDAnalysis.lib.distances.calc_bonds(a, b,
                                                    box=triclinic_box, backend=backend)
        reference = np.array([0.0, 1.7320508, 1.4142136, 2.82842712])
        assert_almost_equal(dists, reference, self.prec, err_msg="calc_bonds with triclinic box failed")

    def test_angles(self, positions, backend):
        a, b, c, d = positions
        angles = MDAnalysis.lib.distances.calc_angles(a, b, c,
                                                      backend=backend)
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

    # Check data type checks
    def test_angles_wrongtype(self, positions, wrongtype, wronglength, backend):
        a, b, c, d = positions
        with pytest.raises(TypeError):
            MDAnalysis.lib.distances.calc_angles(
                a, wrongtype, c, backend=backend)

        with pytest.raises(TypeError):
            MDAnalysis.lib.distances.calc_angles()
            wrongtype, b, c, backend = backend

        with pytest.raises(TypeError):
            MDAnalysis.lib.distances.calc_angles(
                a, b, wrongtype, backend=backend)

        with pytest.raises(ValueError):
            MDAnalysis.lib.distances.calc_angles(
                a, wronglength, c,
                backend=backend)  # try inputting arrays of different length

        with pytest.raises(ValueError):
            MDAnalysis.lib.distances.calc_angles(
                wronglength, b, c,
                backend=backend)

        with pytest.raises(ValueError):
            MDAnalysis.lib.distances.calc_angles(
                a, b, wronglength,
                backend=backend)

    def test_angles_bad_result(self, positions, backend):
        a, b, c, d = positions
        badresult = np.zeros(len(a) - 1)
        with pytest.raises(ValueError):
            MDAnalysis.lib.distances.calc_angles(
                a, b, c, result=badresult, backend=backend)  # Bad result array

    def test_dihedrals(self, positions, backend):
        a, b, c, d = positions
        dihedrals = MDAnalysis.lib.distances.calc_dihedrals(a, b,
                                                            c, d,
                                                            backend=backend)
        # Check calculated values
        assert_equal(len(dihedrals), 4, err_msg="calc_dihedrals results have wrong length")
        assert np.isnan(dihedrals[0]), "Zero length dihedral failed"
        assert np.isnan(dihedrals[1]), "Straight line dihedral failed"
        assert_almost_equal(dihedrals[2], np.pi, self.prec, err_msg="180 degree dihedral failed")
        assert_almost_equal(dihedrals[3], -0.50714064, self.prec,
                            err_msg="arbitrary dihedral angle failed")

    # Check data type checks
    def test_dihedrals_wrongtype(self, positions, wrongtype, backend):
        a, b, c, d = positions
        with pytest.raises(TypeError):
            MDAnalysis.lib.distances.calc_dihedrals(
                a, wrongtype, c, d,
                backend=backend)  # try inputting float64 values

        with pytest.raises(TypeError):
            MDAnalysis.lib.distances.calc_dihedrals(
                wrongtype, b, c, d,
                backend=backend)

        with pytest.raises(TypeError):
            MDAnalysis.lib.distances.calc_dihedrals(
                a, b, wrongtype, d,
                backend=backend)

        with pytest.raises(TypeError):
            MDAnalysis.lib.distances.calc_dihedrals(
                a, b, c, wrongtype,
                backend=backend)

    def test_dihedrals_wronglength(self, positions, wronglength, backend):
        a, b, c, d = positions
        with pytest.raises(ValueError):
            MDAnalysis.lib.distances.calc_dihedrals(
                a, wronglength, c, d,
                backend=backend)

        with pytest.raises(ValueError):
            MDAnalysis.lib.distances.calc_dihedrals(
                wronglength, b, c, d,
                backend=backend)

        with pytest.raises(ValueError):
            MDAnalysis.lib.distances.calc_dihedrals(
                a, b, wronglength, d,
                backend=backend)

        with pytest.raises(ValueError):
            MDAnalysis.lib.distances.calc_dihedrals(
                a, b, c, wronglength,
                backend=backend)

    def test_dihedrals_bad_result(self, positions, backend):
        a, b, c, d = positions
        badresult = np.zeros(len(a) - 1)

        with pytest.raises(ValueError):
            MDAnalysis.lib.distances.calc_dihedrals(
                a, b, c, d, result=badresult,
                backend=backend)  # Bad result array



    def test_numpy_compliance(self, positions, backend):
        a, b, c, d = positions
        # Checks that the cython functions give identical results to the numpy versions
        bonds = MDAnalysis.lib.distances.calc_bonds(a, b,
                                                    backend=backend)
        angles = MDAnalysis.lib.distances.calc_angles(a, b, c,
                                                      backend=backend)
        dihedrals = MDAnalysis.lib.distances.calc_dihedrals(a, b,
                                                            c, d,
                                                            backend=backend)

        bonds_numpy = np.array([mdamath.norm(y - x) for x, y in zip(a, b)])
        vec1 = a - b
        vec2 = c - b
        angles_numpy = np.array([mdamath.angle(x, y) for x, y in zip(vec1, vec2)])
        ab = a - b
        bc = b - c
        cd = c - d
        dihedrals_numpy = np.array([mdamath.dihedral(x, y, z) for x, y, z in zip(ab, bc, cd)])

        assert_almost_equal(bonds, bonds_numpy, self.prec,
                            err_msg="Cython bonds didn't match numpy calculations")
        # numpy 0 angle returns NaN rather than 0
        assert_almost_equal(angles[1:], angles_numpy[1:], self.prec,
                            err_msg="Cython angles didn't match numpy calcuations")
        assert_almost_equal(dihedrals, dihedrals_numpy, self.prec,
                            err_msg="Cython dihedrals didn't match numpy calculations")


@pytest.mark.parametrize('backend', ['serial', 'openmp'])
class Test_apply_PBC(object):
    prec = 6

    def test_ortho_PBC(self, backend):
        from MDAnalysis.lib.distances import apply_PBC

        U = MDAnalysis.Universe(PSF, DCD)
        atoms = U.atoms.positions
        box1 = np.array([2.5, 2.5, 3.5], dtype=np.float32)
        box2 = np.array([2.5, 2.5, 3.5, 90., 90., 90.], dtype=np.float32)

        cyth1 = apply_PBC(atoms, box1, backend=backend)
        cyth2 = apply_PBC(atoms, box2, backend=backend)
        reference = atoms - np.floor(atoms / box1) * box1

        assert_almost_equal(cyth1, reference, self.prec,
                            err_msg="Ortho apply_PBC #1 failed comparison with np")
        assert_almost_equal(cyth2, reference, self.prec,
                            err_msg="Ortho apply_PBC #2 failed comparison with np")

    def test_tric_PBC(self, backend):
        from MDAnalysis.lib.distances import apply_PBC

        U = MDAnalysis.Universe(TRIC)
        atoms = U.atoms.positions
        box1 = U.dimensions
        box2 = MDAnalysis.coordinates.core.triclinic_vectors(box1)

        def numpy_PBC(coords, box):
            # move to fractional coordinates
            fractional = MDAnalysis.lib.distances.transform_RtoS(coords, box)
            # move fractional coordinates to central cell
            fractional -= np.floor(fractional)
            # move back to real coordinates
            return MDAnalysis.lib.distances.transform_StoR(fractional, box)

        cyth1 = apply_PBC(atoms, box1, backend=backend)
        cyth2 = apply_PBC(atoms, box2, backend=backend)
        reference = numpy_PBC(atoms, box2)

        assert_almost_equal(cyth1, reference, decimal=4,
                            err_msg="Triclinic apply_PBC failed comparison with np")
        assert_almost_equal(cyth2, reference, decimal=4,
                            err_msg="Triclinic apply_PBC failed comparison with np")

        box = np.array([10, 7, 3, 45, 60, 90], dtype=np.float32)
        r = np.array([[5.75, 0.36066014, 0.75]], dtype=np.float32)
        r_in_cell = MDAnalysis.lib.distances.apply_PBC(r, box)[0]
        assert_almost_equal([5.75, 7.3606596, 0.75],
                            r_in_cell, self.prec)

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
        from MDAnalysis.lib.distances import calc_angles
        # Shift atom coordinates a few box lengths in random directions and see if we still get same results
        a, b, c, d, box = positions
        a2 = (a + box * (-1, 0, 0)).astype(np.float32)  # seem to get converted to float64 otherwise
        b2 = (b + box * (1, 0, 1)).astype(np.float32)
        c2 = (c + box * (-2, 5, -7)).astype(np.float32)

        ref = calc_angles(a, b, c, backend=backend)

        test1 = calc_angles(a2, b, c, box=box, backend=backend)
        test2 = calc_angles(a, b2, c, box=box, backend=backend)
        test3 = calc_angles(a, b, c2, box=box, backend=backend)
        test4 = calc_angles(a2, b2, c2, box=box, backend=backend)

        for val in [test1, test2, test3, test4]:
            assert_almost_equal(ref, val, self.prec, err_msg="Min image in angle calculation failed")

    def test_dihedrals(self, positions, backend):
        from MDAnalysis.lib.distances import calc_dihedrals
        a, b, c, d, box = positions
        a2 = (a + box * (-1, 0, 0)).astype(np.float32)
        b2 = (b + box * (1, 0, 1)).astype(np.float32)
        c2 = (c + box * (-2, 5, -7)).astype(np.float32)
        d2 = (d + box * (0, -5, 0)).astype(np.float32)

        ref = calc_dihedrals(a, b, c, d, backend=backend)

        test1 = calc_dihedrals(a2, b, c, d, box=box,
                               backend=backend)
        test2 = calc_dihedrals(a, b2, c, d, box=box,
                               backend=backend)
        test3 = calc_dihedrals(a, b, c2, d, box=box,
                               backend=backend)
        test4 = calc_dihedrals(a, b, c, d2, box=box,
                               backend=backend)
        test5 = calc_dihedrals(a2, b2, c2, d2, box=box,
                               backend=backend)

        for val in [test1, test2, test3, test4, test5]:
            assert_almost_equal(ref, val, self.prec, err_msg="Min image in dihedral calculation failed")


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
            MDAnalysis.lib.distances._run("calc_self_distance_array",
                                          args=(positions, result),
                                          backend=backend)
        except RuntimeError:
            pytest.fail("Failed to understand backend {0}".format(backend))

    def test_wront_backend(self, backend_selection_pos):
        positions, result = backend_selection_pos
        with pytest.raises(ValueError):
            MDAnalysis.lib.distances._run("calc_self_distance_array",
                                          args=(positions, result),
                                          backend="not implemented stuff")

def test_used_openmpflag():
    assert isinstance(MDAnalysis.lib.distances.USED_OPENMP, bool)
