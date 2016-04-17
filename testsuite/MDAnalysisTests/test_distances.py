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
import MDAnalysis
import MDAnalysis.lib.distances

import numpy as np
from numpy.testing import (TestCase, dec, raises, assert_,
                           assert_almost_equal, assert_equal, assert_raises,)

from nose.plugins.attrib import attr

from MDAnalysis.tests.datafiles import PSF, DCD, TRIC
from MDAnalysis.lib import mdamath
from MDAnalysisTests import parser_not_found

class _TestDistanceArray(TestCase):
    # override backend in test classes
    backend = None
    def setUp(self):
        self.box = np.array([1., 1., 2.], dtype=np.float32)
        self.points = np.array(
            [
                [0, 0, 0], [1, 1, 2], [1, 0, 2],  # identical under PBC
                [0.5, 0.5, 1.5],
            ], dtype=np.float32)
        self.ref = self.points[0:1]
        self.conf = self.points[1:]

    def _dist(self, n, ref=None):
        if ref is None:
            ref = self.ref[0]
        else:
            ref = np.asarray(ref, dtype=np.float32)
        x = self.points[n]
        r = x - ref
        return np.sqrt(np.dot(r, r))

    def test_noPBC(self):
        d = MDAnalysis.lib.distances.distance_array(self.ref, self.points,
                                                    backend=self.backend)
        assert_almost_equal(d, np.array([[self._dist(0), self._dist(1), self._dist(2), self._dist(3)]]))

    def test_PBC(self):
        d = MDAnalysis.lib.distances.distance_array(self.ref, self.points,
                                                    box=self.box, backend=self.backend)
        assert_almost_equal(d, np.array([[0., 0., 0., self._dist(3, ref=[1, 1, 2])]]))

    def test_PBC2(self):
        a = np.array([7.90146923, -13.72858524, 3.75326586], dtype=np.float32)
        b = np.array([-1.36250901, 13.45423985, -0.36317623], dtype=np.float32)
        box = np.array([5.5457325, 5.5457325, 5.5457325], dtype=np.float32)

        def mindist(a, b, box):
            x = a - b
            return np.linalg.norm(x - np.rint(x / box) * box)

        ref = mindist(a, b, box)
        val = MDAnalysis.lib.distances.distance_array(np.array([a]), np.array([b]),
                                                      box=box, backend=self.backend)[0, 0]

        assert_almost_equal(val, ref, decimal=6,
                            err_msg="Issue 151 not correct (PBC in distance array)")

class TestDistanceArray_Serial(_TestDistanceArray):
    backend = "serial"

class TestDistanceArray_OpenMP(_TestDistanceArray):
    backend = "OpenMP"

class _TestDistanceArrayDCD(TestCase):
    backend = None

    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.trajectory = self.universe.trajectory
        self.ca = self.universe.select_atoms('name CA')
        # reasonable precision so that tests succeed on 32 and 64 bit machines
        # (the reference values were obtained on 64 bit)
        # Example:
        #   Items are not equal: wrong maximum distance value
        #   ACTUAL: 52.470254967456412
        #   DESIRED: 52.470257062419059
        self.prec = 5

    def tearDown(self):
        del self.universe
        del self.trajectory
        del self.ca

    @attr('issue')
    def test_simple(self):
        U = self.universe
        self.trajectory.rewind()
        x0 = U.atoms.positions
        self.trajectory[10]
        x1 = U.atoms.positions
        d = MDAnalysis.lib.distances.distance_array(x0, x1, backend=self.backend)
        assert_equal(d.shape, (3341, 3341), "wrong shape (should be (Natoms,Natoms))")
        assert_almost_equal(d.min(), 0.11981228170520701, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 53.572192429459619, self.prec,
                            err_msg="wrong maximum distance value")

    def test_outarray(self):
        U = self.universe
        self.trajectory.rewind()
        x0 = U.atoms.positions
        self.trajectory[10]
        x1 = U.atoms.positions
        natoms = len(U.atoms)
        d = np.zeros((natoms, natoms), np.float64)
        MDAnalysis.lib.distances.distance_array(x0, x1, result=d, backend=self.backend)
        assert_equal(d.shape, (natoms, natoms), "wrong shape, shoud be  (Natoms,Natoms) entries")
        assert_almost_equal(d.min(), 0.11981228170520701, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 53.572192429459619, self.prec,
                            err_msg="wrong maximum distance value")

    def test_periodic(self):
        # boring with the current dcd as that has no PBC
        U = self.universe
        self.trajectory.rewind()
        x0 = U.atoms.positions
        self.trajectory[10]
        x1 = U.atoms.positions
        d = MDAnalysis.lib.distances.distance_array(x0, x1, box=U.coord.dimensions,
                                                    backend=self.backend)
        assert_equal(d.shape, (3341, 3341), "should be square matrix with Natoms entries")
        assert_almost_equal(d.min(), 0.11981228170520701, self.prec,
                            err_msg="wrong minimum distance value with PBC")
        assert_almost_equal(d.max(), 53.572192429459619, self.prec,
                            err_msg="wrong maximum distance value with PBC")


class TestDistanceArrayDCD_Serial(_TestDistanceArrayDCD):
    backend = "serial"

class TestDistanceArrayDCD_OpenMP(_TestDistanceArrayDCD):
    backend = "OpenMP"

class _TestSelfDistanceArrayDCD(TestCase):
    backend = None

    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.trajectory = self.universe.trajectory
        self.ca = self.universe.select_atoms('name CA')
        # see comments above on precision
        self.prec = 5

    def tearDown(self):
        del self.universe
        del self.trajectory
        del self.ca

    def test_simple(self):
        U = self.universe
        self.trajectory.rewind()
        x0 = U.atoms.positions
        d = MDAnalysis.lib.distances.self_distance_array(x0, backend=self.backend)
        N = 3341 * (3341 - 1) / 2
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 52.4702570624190590, self.prec,
                            err_msg="wrong maximum distance value")

    def test_outarray(self):
        U = self.universe
        self.trajectory.rewind()
        x0 = U.atoms.positions
        natoms = len(U.atoms)
        N = natoms * (natoms - 1) / 2
        d = np.zeros((N,), np.float64)
        MDAnalysis.lib.distances.self_distance_array(x0, result=d, backend=self.backend)
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 52.4702570624190590, self.prec,
                            err_msg="wrong maximum distance value")

    def test_periodic(self):
        # boring with the current dcd as that has no PBC
        U = self.universe
        self.trajectory.rewind()
        x0 = U.atoms.positions
        natoms = len(U.atoms)
        N = natoms * (natoms - 1) / 2
        d = MDAnalysis.lib.distances.self_distance_array(x0, box=U.coord.dimensions,
                                                         backend=self.backend)
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, self.prec,
                            err_msg="wrong minimum distance value with PBC")
        assert_almost_equal(d.max(), 52.4702570624190590, self.prec,
                            err_msg="wrong maximum distance value with PBC")

class TestSelfDistanceArrayDCD_Serial(_TestSelfDistanceArrayDCD):
    backend = "serial"

class TestSelfDistanceArrayDCD_OpenMP(_TestSelfDistanceArrayDCD):
    backend = "OpenMP"

class _TestTriclinicDistances(TestCase):
    """Unit tests for the Triclinic PBC functions.
    Tests:
      # transforming to and from S space (fractional coords)
      mda.lib.distances.transform_StoR
      mda.lib.distances.transform_RtoS
      # distance calculations with PBC
      mda.lib.distances.self_distance_array
      mda.lib.distances.distance_array
    """
    backend = None

    def setUp(self):
        self.universe = MDAnalysis.Universe(TRIC)
        self.prec = 2

        self.box = MDAnalysis.coordinates.core.triclinic_vectors(self.universe.dimensions)
        self.boxV = MDAnalysis.coordinates.core.triclinic_box(self.box[0], self.box[1], self.box[2])

        self.S_mol1 = np.array([self.universe.atoms[383].pos])
        self.S_mol2 = np.array([self.universe.atoms[390].pos])

    def tearDown(self):
        del self.universe
        del self.boxV
        del self.box
        del self.S_mol1
        del self.S_mol2
        del self.prec

    def test_transforms(self):
        from MDAnalysis.lib.distances import transform_StoR, transform_RtoS
        # To check the cython coordinate transform, the same operation is done in numpy
        # Is a matrix multiplication of Coords x Box = NewCoords, so can use np.dot

        # Test transformation
        R_mol1 = transform_StoR(self.S_mol1, self.box, backend=self.backend)
        R_np1 = np.dot(self.S_mol1, self.box)

        # Test transformation when given box in different form
        R_mol2 = transform_StoR(self.S_mol2, self.boxV, backend=self.backend)
        R_np2 = np.dot(self.S_mol2, self.box)

        assert_almost_equal(R_mol1, R_np1, self.prec, err_msg="StoR transform failed with box")
        assert_almost_equal(R_mol2, R_np2, self.prec, err_msg="StoR transform failed with boxV")

        # Round trip test
        # boxV here althought initial transform with box
        S_test1 = transform_RtoS(R_mol1, self.boxV, backend=self.backend)
        # and vice versa, should still work
        S_test2 = transform_RtoS(R_mol2, self.box, backend=self.backend)

        assert_almost_equal(S_test1, self.S_mol1, self.prec, err_msg="Round trip failed in transform")
        assert_almost_equal(S_test2, self.S_mol2, self.prec, err_msg="Round trip failed in transform")

    def test_selfdist(self):
        from MDAnalysis.lib.distances import self_distance_array
        from MDAnalysis.lib.distances import transform_RtoS, transform_StoR

        R_coords = transform_StoR(self.S_mol1, self.box, backend=self.backend)
        # Transform functions are tested elsewhere so taken as working here
        dists = self_distance_array(R_coords, box=self.box, backend=self.backend)
        # Manually calculate self_distance_array
        manual = np.zeros(len(dists), dtype=np.float64)
        distpos = 0
        for i, Ri in enumerate(R_coords):
            for Rj in R_coords[i + 1:]:
                Rij = Rj - Ri
                Rij -= round(Rij[2] / self.box[2][2]) * self.box[2]
                Rij -= round(Rij[1] / self.box[1][1]) * self.box[1]
                Rij -= round(Rij[0] / self.box[0][0]) * self.box[0]
                Rij = np.linalg.norm(Rij)  # find norm of Rij vector
                manual[distpos] = Rij  # and done, phew
                distpos += 1

        assert_almost_equal(dists, manual, self.prec,
                            err_msg="self_distance_array failed with input 1")

        # Do it again for input 2 (has wider separation in points)
        # Also use boxV here in self_dist calculation
        R_coords = transform_StoR(self.S_mol2, self.box, backend=self.backend)
        # Transform functions are tested elsewhere so taken as working here
        dists = self_distance_array(R_coords, box=self.boxV, backend=self.backend)
        # Manually calculate self_distance_array
        manual = np.zeros(len(dists), dtype=np.float64)
        distpos = 0
        for i, Ri in enumerate(R_coords):
            for Rj in R_coords[i + 1:]:
                Rij = Rj - Ri
                Rij -= round(Rij[2] / self.box[2][2]) * self.box[2]
                Rij -= round(Rij[1] / self.box[1][1]) * self.box[1]
                Rij -= round(Rij[0] / self.box[0][0]) * self.box[0]
                Rij = np.linalg.norm(Rij)  # find norm of Rij vector
                manual[distpos] = Rij  # and done, phew
                distpos += 1

        assert_almost_equal(dists, manual, self.prec,
                            err_msg="self_distance_array failed with input 2")

    def test_distarray(self):
        from MDAnalysis.lib.distances import distance_array
        from MDAnalysis.lib.distances import transform_StoR, transform_RtoS

        R_mol1 = transform_StoR(self.S_mol1, self.box, backend=self.backend)
        R_mol2 = transform_StoR(self.S_mol2, self.box, backend=self.backend)

        # Try with box
        dists = distance_array(R_mol1, R_mol2, box=self.box, backend=self.backend)
        # Manually calculate distance_array
        manual = np.zeros((len(R_mol1), len(R_mol2)))
        for i, Ri in enumerate(R_mol1):
            for j, Rj in enumerate(R_mol2):
                Rij = Rj - Ri
                Rij -= round(Rij[2] / self.box[2][2]) * self.box[2]
                Rij -= round(Rij[1] / self.box[1][1]) * self.box[1]
                Rij -= round(Rij[0] / self.box[0][0]) * self.box[0]
                Rij = np.linalg.norm(Rij)  # find norm of Rij vector
                manual[i][j] = Rij

        assert_almost_equal(dists, manual, self.prec,
                            err_msg="distance_array failed with box")

        # Now check using boxV
        dists = distance_array(R_mol1, R_mol2, box=self.boxV, backend=self.backend)
        assert_almost_equal(dists, manual, self.prec,
                            err_msg="distance_array failed with boxV")

    def test_pbc_dist(self):
        from MDAnalysis.lib.distances import distance_array

        results = np.array([[37.629944]])

        dists = distance_array(self.S_mol1, self.S_mol2, box=self.boxV,
                               backend=self.backend)

        assert_almost_equal(dists, results, self.prec,
                            err_msg="distance_array failed to retrieve PBC distance")

class TestTriclinicDistances_Serial(_TestTriclinicDistances):
    backend = "serial"

class TestTriclinicDistances_OpenMP(_TestTriclinicDistances):
    backend = "OpenMP"

class _TestCythonFunctions(TestCase):
    # Unit tests for calc_bonds calc_angles and calc_dihedrals in lib.distances
    # Tests both numerical results as well as input types as Cython will silently
    # produce nonsensical results if given wrong data types otherwise.
    backend = None

    def setUp(self):
        self.prec = 5
        self.box = np.array([10., 10., 10.], dtype=np.float32)
        self.box2 = np.array([[10., 0., 0.], [1., 10., 0., ], [1., 0., 10.]], dtype=np.float32)
        # dummy atom data
        self.a = np.array([[0., 0., 0.], [0., 0., 0.], [0., 11., 0.], [1., 1., 1.]], dtype=np.float32)
        self.b = np.array([[0., 0., 0.], [1., 1., 1.], [0., 0., 0.], [29., -21., 99.]], dtype=np.float32)
        self.c = np.array([[0., 0., 0.], [2., 2., 2.], [11., 0., 0.], [1., 9., 9.]], dtype=np.float32)
        self.d = np.array([[0., 0., 0.], [3., 3., 3.], [11., -11., 0.], [65., -65., 65.]], dtype=np.float32)
        self.wrongtype = np.array([[0., 0., 0.], [3., 3., 3.], [3., 3., 3.], [3., 3., 3.]],
                                  dtype=np.float64)  # declared as float64 and should raise TypeError
        self.wronglength = np.array([[0., 0., 0.], [3., 3., 3.]],
                                    dtype=np.float32)  # has a different length to other inputs and should raise
                                    # ValueError

    def tearDown(self):
        del self.box
        del self.box2
        del self.a
        del self.b
        del self.c
        del self.d
        del self.wrongtype
        del self.wronglength

    def test_bonds(self):
        dists = MDAnalysis.lib.distances.calc_bonds(self.a, self.b,
                                                    backend=self.backend)
        assert_equal(len(dists), 4, err_msg="calc_bonds results have wrong length")
        dists_pbc = MDAnalysis.lib.distances.calc_bonds(self.a, self.b, box=self.box,
                                                        backend=self.backend)
        #tests 0 length
        assert_almost_equal(dists[0], 0.0, self.prec, err_msg="Zero length calc_bonds fail")
        assert_almost_equal(dists[1], 1.7320508075688772, self.prec,
                            err_msg="Standard length calc_bonds fail")  # arbitrary length check
        #PBC checks, 2 without, 2 with
        assert_almost_equal(dists[2], 11.0, self.prec,
                            err_msg="PBC check #1 w/o box")  # pbc check 1, subtract single box length
        assert_almost_equal(dists_pbc[2], 1.0, self.prec,
                            err_msg="PBC check #1 with box")
        assert_almost_equal(dists[3], 104.26888318, self.prec,  # pbc check 2, subtract multiple box
                            err_msg="PBC check #2 w/o box")  # lengths in all directions
        assert_almost_equal(dists_pbc[3], 3.46410072, self.prec,
                            err_msg="PBC check #w with box")
        #Bad input checking
    def test_bonds_wrongtype(self):
        assert_raises(TypeError, MDAnalysis.lib.distances.calc_bonds, self.a,
                      self.wrongtype, backend=self.backend)
        assert_raises(TypeError, MDAnalysis.lib.distances.calc_bonds,
                      self.wrongtype, self.b, backend=self.backend)
        assert_raises(ValueError, MDAnalysis.lib.distances.calc_bonds, self.a,
                      self.wronglength, backend=self.backend)
        assert_raises(ValueError, MDAnalysis.lib.distances.calc_bonds,
                      self.wronglength, self.b, backend=self.backend)

    def test_bonds_badbox(self):
        badboxtype = np.array([10., 10., 10.], dtype=np.float64)
        badboxsize = np.array([[10., 10.], [10., 10., ]], dtype=np.float32)

        assert_raises(ValueError, MDAnalysis.lib.distances.calc_bonds, self.a,
                      self.b, box=badboxsize, backend=self.backend)  # Bad box data
        assert_raises(TypeError, MDAnalysis.lib.distances.calc_bonds, self.a,
                      self.b, box=badboxtype, backend=self.backend)  # Bad box type

    def test_bonds_badresult(self):
        badresult = np.zeros(len(self.a) - 1)
        assert_raises(ValueError, MDAnalysis.lib.distances.calc_bonds,
                      self.a, self.b, result=badresult, backend=self.backend)  # Bad result array

    def test_bonds_triclinic(self):
        dists = MDAnalysis.lib.distances.calc_bonds(self.a, self.b,
                                                    box=self.box2, backend=self.backend)
        reference = np.array([0.0, 1.7320508, 1.4142136, 2.82842712])
        assert_almost_equal(dists, reference, self.prec, err_msg="calc_bonds with triclinic box failed")

    def test_angles(self):
        angles = MDAnalysis.lib.distances.calc_angles(self.a, self.b, self.c,
                                                      backend=self.backend)
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
    def test_angles_wrongtype(self):
        assert_raises(TypeError, MDAnalysis.lib.distances.calc_angles,
                      self.a, self.wrongtype, self.c, backend=self.backend)  # try inputting float64 values
        assert_raises(TypeError, MDAnalysis.lib.distances.calc_angles,
                      self.wrongtype, self.b, self.c, backend=self.backend)
        assert_raises(TypeError, MDAnalysis.lib.distances.calc_angles,
                      self.a, self.b, self.wrongtype, backend=self.backend)
        assert_raises(ValueError, MDAnalysis.lib.distances.calc_angles,
                      self.a, self.wronglength, self.c,
                      backend=self.backend)  # try inputting arrays of different length
        assert_raises(ValueError, MDAnalysis.lib.distances.calc_angles,
                      self.wronglength, self.b, self.c,
                      backend=self.backend)
        assert_raises(ValueError, MDAnalysis.lib.distances.calc_angles,
                      self.a, self.b, self.wronglength,
                      backend=self.backend)

    def test_angles_bad_result(self):
        badresult = np.zeros(len(self.a) - 1)
        assert_raises(ValueError, MDAnalysis.lib.distances.calc_angles,
                      self.a, self.b, self.c, result=badresult, backend=self.backend)  # Bad result array

    def test_dihedrals(self):
        dihedrals = MDAnalysis.lib.distances.calc_dihedrals(self.a, self.b,
                                                            self.c, self.d,
                                                            backend=self.backend)
        # Check calculated values
        assert_equal(len(dihedrals), 4, err_msg="calc_dihedrals results have wrong length")
        #        assert_almost_equal(dihedrals[0], 0.0, self.prec, err_msg="Zero length dihedral failed")
        #        assert_almost_equal(dihedrals[1], 0.0, self.prec, err_msg="Straight line dihedral failed")
        assert_almost_equal(dihedrals[2], np.pi, self.prec, err_msg="180 degree dihedral failed")
        assert_almost_equal(dihedrals[3], 0.50714064, self.prec,
                            err_msg="arbitrary dihedral angle failed")
    # Check data type checks
    def test_dihedrals_wrongtype(self):
        assert_raises(TypeError, MDAnalysis.lib.distances.calc_dihedrals,
                      self.a, self.wrongtype, self.c, self.d,
                      backend=self.backend)  # try inputting float64 values
        assert_raises(TypeError, MDAnalysis.lib.distances.calc_dihedrals,
                      self.wrongtype, self.b, self.c, self.d,
                      backend=self.backend)
        assert_raises(TypeError, MDAnalysis.lib.distances.calc_dihedrals,
                      self.a, self.b, self.wrongtype, self.d,
                      backend=self.backend)
        assert_raises(TypeError, MDAnalysis.lib.distances.calc_dihedrals,
                      self.a, self.b, self.c, self.wrongtype,
                      backend=self.backend)

    def test_dihedrals_wronglength(self):
        assert_raises(ValueError, MDAnalysis.lib.distances.calc_dihedrals,
                      self.a, self.wronglength, self.c, self.d,
                      backend=self.backend)
        assert_raises(ValueError, MDAnalysis.lib.distances.calc_dihedrals,
                      self.wronglength, self.b, self.c, self.d,
                      backend=self.backend)
        assert_raises(ValueError, MDAnalysis.lib.distances.calc_dihedrals,
                      self.a, self.b, self.wronglength, self.d,
                      backend=self.backend)
        assert_raises(ValueError, MDAnalysis.lib.distances.calc_dihedrals,
                      self.a, self.b, self.c, self.wronglength,
                      backend=self.backend)

    def test_dihedrals_bad_result(self):
        badresult = np.zeros(len(self.a) - 1)

        assert_raises(ValueError, MDAnalysis.lib.distances.calc_dihedrals,
                      self.a, self.b, self.c, self.d, result=badresult,
                      backend=self.backend)  # Bad result array

    def test_numpy_compliance(self):
        # Checks that the cython functions give identical results to the numpy versions
        bonds = MDAnalysis.lib.distances.calc_bonds(self.a, self.b,
                                                    backend=self.backend)
        angles = MDAnalysis.lib.distances.calc_angles(self.a, self.b, self.c,
                                                      backend=self.backend)
        dihedrals = MDAnalysis.lib.distances.calc_dihedrals(self.a, self.b,
                                                            self.c, self.d,
                                                            backend=self.backend)

        bonds_numpy = np.array([mdamath.norm(y - x) for x, y in zip(self.a, self.b)])
        vec1 = self.a - self.b
        vec2 = self.c - self.b
        angles_numpy = np.array([mdamath.angle(x, y) for x, y in zip(vec1, vec2)])
        ab = self.b - self.a
        bc = self.c - self.b
        cd = self.d - self.c
        dihedrals_numpy = np.array([mdamath.dihedral(x, y, z) for x, y, z in zip(ab, bc, cd)])

        assert_almost_equal(bonds, bonds_numpy, self.prec,
                            err_msg="Cython bonds didn't match numpy calculations")
        # numpy 0 angle returns NaN rather than 0
        assert_almost_equal(angles[1:], angles_numpy[1:], self.prec,
                            err_msg="Cython angles didn't match numpy calcuations")
        # same issue with first two dihedrals
        assert_almost_equal(dihedrals[2:], dihedrals_numpy[2:], self.prec,
                            err_msg="Cython dihedrals didn't match numpy calculations")

class TestCythonFunctions_Serial(_TestCythonFunctions):
    backend = "serial"

class TestCythonFunctions_OpenMP(_TestCythonFunctions):
    backend = "OpenMP"

class _Test_apply_PBC(TestCase):
    backend = None

    def setUp(self):
        self.prec = 6

    def tearDown(self):
        del self.prec

    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_ortho_PBC(self):
        from MDAnalysis.lib.distances import apply_PBC

        U = MDAnalysis.Universe(PSF, DCD)
        atoms = U.atoms.positions
        box1 = np.array([2.5, 2.5, 3.5], dtype=np.float32)
        box2 = np.array([2.5, 2.5, 3.5, 90., 90., 90.], dtype=np.float32)

        cyth1 = apply_PBC(atoms, box1, backend=self.backend)
        cyth2 = apply_PBC(atoms, box2, backend=self.backend)
        reference = atoms - np.floor(atoms / box1) * box1

        assert_almost_equal(cyth1, reference, self.prec,
                            err_msg="Ortho apply_PBC #1 failed comparison with np")
        assert_almost_equal(cyth2, reference, self.prec,
                            err_msg="Ortho apply_PBC #2 failed comparison with np")

    def test_tric_PBC(self):
        from MDAnalysis.lib.distances import apply_PBC

        U = MDAnalysis.Universe(TRIC)
        atoms = U.atoms.positions
        box1 = U.dimensions
        box2 = MDAnalysis.coordinates.core.triclinic_vectors(box1)

        #print box2
        #print box2.shape

        def numpy_PBC(coords, box):
            coords -= np.array([box[2] * val for val in np.floor(coords[:, 2] / box[2][2])])
            coords -= np.array([box[1] * val for val in np.floor(coords[:, 1] / box[1][1])])
            coords -= np.array([box[0] * val for val in np.floor(coords[:, 0] / box[0][0])])

            return coords

        cyth1 = apply_PBC(atoms, box1, backend=self.backend)
        cyth2 = apply_PBC(atoms, box2, backend=self.backend)
        reference = numpy_PBC(atoms, box2)

        assert_almost_equal(cyth1, reference, self.prec,
                            err_msg="Triclinic apply_PBC failed comparison with np")
        assert_almost_equal(cyth2, reference, self.prec,
                            err_msg="Trlclinic apply_PBC failed comparison with np")

class _Test_apply_PBC_Serial(_Test_apply_PBC):
    backend = "serial"

class _Test_apply_PBC_OpenMP(_Test_apply_PBC):
    backend = "OpenMP"

class _TestPeriodicAngles(TestCase):
    """Test case for properly considering minimum image convention when calculating angles and dihedrals
    (Issue 172)
    """

    def setUp(self):
        self.prec = 5
        self.a = np.array([[0.0, 1.0, 0.0]], dtype=np.float32)
        self.b = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        self.c = np.array([[1.0, 0.0, 0.0]], dtype=np.float32)
        self.d = np.array([[1.0, 0.0, 1.0]], dtype=np.float32)
        self.box = np.array([10.0, 10.0, 10.0], dtype=np.float32)

    def tearDown(self):
        del self.prec
        del self.a
        del self.b
        del self.c
        del self.d
        del self.box

    def test_angles(self):
        from MDAnalysis.lib.distances import calc_angles
        # Shift atom coordinates a few box lengths in random directions and see if we still get same results
        a2 = (self.a + self.box * (-1, 0, 0)).astype(np.float32)  # seem to get converted to float64 otherwise
        b2 = (self.b + self.box * (1, 0, 1)).astype(np.float32)
        c2 = (self.c + self.box * (-2, 5, -7)).astype(np.float32)

        ref = calc_angles(self.a, self.b, self.c, backend=self.backend)

        test1 = calc_angles(a2, self.b, self.c, box=self.box, backend=self.backend)
        test2 = calc_angles(self.a, b2, self.c, box=self.box, backend=self.backend)
        test3 = calc_angles(self.a, self.b, c2, box=self.box, backend=self.backend)
        test4 = calc_angles(a2, b2, c2, box=self.box, backend=self.backend)

        for val in [test1, test2, test3, test4]:
            assert_almost_equal(ref, val, self.prec, err_msg="Min image in angle calculation failed")

    def test_dihedrals(self):
        from MDAnalysis.lib.distances import calc_dihedrals

        a2 = (self.a + self.box * (-1, 0, 0)).astype(np.float32)
        b2 = (self.b + self.box * (1, 0, 1)).astype(np.float32)
        c2 = (self.c + self.box * (-2, 5, -7)).astype(np.float32)
        d2 = (self.d + self.box * (0, -5, 0)).astype(np.float32)

        ref = calc_dihedrals(self.a, self.b, self.c, self.d, backend=self.backend)

        test1 = calc_dihedrals(a2, self.b, self.c, self.d, box=self.box,
                               backend=self.backend)
        test2 = calc_dihedrals(self.a, b2, self.c, self.d, box=self.box,
                               backend=self.backend)
        test3 = calc_dihedrals(self.a, self.b, c2, self.d, box=self.box,
                               backend=self.backend)
        test4 = calc_dihedrals(self.a, self.b, self.c, d2, box=self.box,
                               backend=self.backend)
        test5 = calc_dihedrals(a2, b2, c2, d2, box=self.box,
                               backend=self.backend)

        for val in [test1, test2, test3, test4, test5]:
            assert_almost_equal(ref, val, self.prec, err_msg="Min image in dihedral calculation failed")

class TestPeriodicAngles_Serial(_TestPeriodicAngles):
    backend = "serial"

class TestPeriodicAngles_OpenMP(_TestPeriodicAngles):
    backend = "OpenMP"

class TestDistanceBackendSelection(object):
    def __init__(self):
        self.positions = np.random.rand(10, 3)
        N = self.positions.shape[0]
        self.result = np.empty(N*(N-1)/2, dtype=np.float64)

    def _case_insensitivity_test(self, backend):
        try:
            MDAnalysis.lib.distances._run("calc_self_distance_array",
                                          args=(self.positions, self.result),
                                          backend=backend)
        except RuntimeError:
            raise AssertionError("Failed to understand backend {0}".format(backend))

    def test_case_insensitivity(self):
        for backend in ("serial", "Serial", "SeRiAL", "SERIAL",
                        "openmp", "OpenMP", "oPENmP", "OPENMP"):
            yield self._case_insensitivity_test, backend

    @raises(ValueError)
    def test_missing_backend_raises_ValueError(self):
        MDAnalysis.lib.distances._run("calc_self_distance_array",
                                      args=(self.positions, self.result),
                                      backend="not_implemented_stuff")


