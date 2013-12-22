# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

import MDAnalysis
import MDAnalysis.core.distances

import numpy as np
from numpy.testing import *
del test
from nose.plugins.attrib import attr

from MDAnalysis.tests.datafiles import PSF, DCD, TRIC



class TestDistanceArray(TestCase):
    def setUp(self):
        self.box = np.array([1.,1.,2.], dtype=np.float32)
        self.points = np.array([[0,0,0], [1,1,2], [1,0,2],  # identical under PBC
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
        return np.sqrt(np.dot(r,r))

    def test_noPBC(self):
        d = MDAnalysis.core.distances.distance_array(self.ref, self.points)
        assert_almost_equal(d, np.array([[self._dist(0), self._dist(1), self._dist(2),
                                          self._dist(3),
                                          ]]))

    def test_PBC(self):
        d = MDAnalysis.core.distances.distance_array(self.ref, self.points, box=self.box)
        assert_almost_equal(d, np.array([[ 0., 0., 0.,
                                           self._dist(3, ref=[1,1,2]),
                                           ]]))

    def test_PBC2(self):
        a = np.array([  7.90146923, -13.72858524,   3.75326586], dtype=np.float32)
        b = np.array([ -1.36250901,  13.45423985,  -0.36317623], dtype=np.float32)
        box = np.array([5.5457325, 5.5457325, 5.5457325], dtype=np.float32)

        def mindist(a, b, box):
            x = a - b
            return np.linalg.norm(x - np.rint(x/box) * box)

        ref = mindist(a, b, box)
        val = MDAnalysis.core.distances.distance_array(np.array([a]), np.array([b]), box)[0,0]

        assert_almost_equal(val, ref, decimal=6, err_msg="Issue 151 not correct (PBC in distance array)")


class TestDistanceArrayDCD(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.trajectory = self.universe.trajectory
        self.ca = self.universe.selectAtoms('name CA')
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
        x0 = U.atoms.coordinates(copy=True)
        self.trajectory[10]
        x1 = U.atoms.coordinates(copy=True)
        d = MDAnalysis.core.distances.distance_array(x0, x1)
        assert_equal(d.shape, (3341, 3341), "wrong shape (should be (Natoms,Natoms))")
        assert_almost_equal(d.min(), 0.11981228170520701, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 53.572192429459619, self.prec,
                            err_msg="wrong maximum distance value")

    def test_outarray(self):
        U = self.universe
        self.trajectory.rewind()
        x0 = U.atoms.coordinates(copy=True)
        self.trajectory[10]
        x1 = U.atoms.coordinates(copy=True)
        natoms = len(U.atoms)
        d = np.zeros((natoms, natoms), np.float64)
        MDAnalysis.core.distances.distance_array(x0, x1, result=d)
        assert_equal(d.shape, (natoms, natoms), "wrong shape, shoud be  (Natoms,Natoms) entries")
        assert_almost_equal(d.min(), 0.11981228170520701, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 53.572192429459619, self.prec,
                            err_msg="wrong maximum distance value")

    def test_periodic(self):
        # boring with the current dcd as that has no PBC
        U = self.universe
        self.trajectory.rewind()
        x0 = U.atoms.coordinates(copy=True)
        self.trajectory[10]
        x1 = U.atoms.coordinates(copy=True)
        d = MDAnalysis.core.distances.distance_array(x0, x1, box=U.coord.dimensions)
        assert_equal(d.shape, (3341, 3341), "should be square matrix with Natoms entries")
        assert_almost_equal(d.min(), 0.11981228170520701, self.prec,
                            err_msg="wrong minimum distance value with PBC")
        assert_almost_equal(d.max(), 53.572192429459619, self.prec,
                            err_msg="wrong maximum distance value with PBC")


class TestSelfDistanceArrayDCD(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.trajectory = self.universe.trajectory
        self.ca = self.universe.selectAtoms('name CA')
        # see comments above on precision
        self.prec = 5

    def tearDown(self):
        del self.universe
        del self.trajectory
        del self.ca

    def test_simple(self):
        U = self.universe
        self.trajectory.rewind()
        x0 = U.atoms.coordinates(copy=True)
        d = MDAnalysis.core.distances.self_distance_array(x0)
        N = 3341 * (3341 - 1) / 2
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 52.4702570624190590, self.prec,
                            err_msg="wrong maximum distance value")

    def test_outarray(self):
        U = self.universe
        self.trajectory.rewind()
        x0 = U.atoms.coordinates(copy=True)
        natoms = len(U.atoms)
        N = natoms*(natoms-1) / 2
        d = np.zeros((N,), np.float64)
        MDAnalysis.core.distances.self_distance_array(x0, result=d)
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, self.prec,
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 52.4702570624190590, self.prec,
                            err_msg="wrong maximum distance value")

    def test_periodic(self):
        # boring with the current dcd as that has no PBC
        U = self.universe
        self.trajectory.rewind()
        x0 = U.atoms.coordinates(copy=True)
        natoms = len(U.atoms)
        N = natoms*(natoms-1) / 2
        d = MDAnalysis.core.distances.self_distance_array(x0, box=U.coord.dimensions)
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, self.prec,
                            err_msg="wrong minimum distance value with PBC")
        assert_almost_equal(d.max(), 52.4702570624190590, self.prec,
                            err_msg="wrong maximum distance value with PBC")



class TestTriclinicDistances(TestCase):
    """Unit tests for the Triclinic PBC functions.
    Tests:
      # transforming to and form S space (fractional coords)
      mda.core.distances.transform_StoR
      mda.core.distances.transform_RtoS
      # distance calculations with PBC 
      mda.core.distances.self_distance_array
      mda.core.distances.distance_array
    """

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
        from MDAnalysis.core.distances import transform_StoR, transform_RtoS
        # To check the cython coordinate transform, the same operation is done in numpy
        # Is a matrix multiplication of Coords x Box = NewCoords, so can use np.dot

        # Test transformation
        R_mol1 = transform_StoR(self.S_mol1, self.box)
        R_np1 = np.dot(self.S_mol1, self.box)
        
        # Test transformation when given box in different form
        R_mol2 = transform_StoR(self.S_mol2, self.boxV)
        R_np2 = np.dot(self.S_mol2, self.box)

        assert_almost_equal(R_mol1, R_np1, self.prec, err_msg="StoR transform failed with box")
        assert_almost_equal(R_mol2, R_np2, self.prec, err_msg="StoR transform failed with boxV")

        # Round trip test
        S_test1 = transform_RtoS(R_mol1, self.boxV) # boxV here althought initial transform with box
        S_test2 = transform_RtoS(R_mol2, self.box) # and vice versa, should still work
        
        assert_almost_equal(S_test1, self.S_mol1, self.prec, err_msg="Round trip failed in transform")
        assert_almost_equal(S_test2, self.S_mol2, self.prec, err_msg="Round trip failed in transform")
                            
    def test_selfdist(self):
        from MDAnalysis.core.distances import self_distance_array
        from MDAnalysis.core.distances import transform_RtoS, transform_StoR

        R_coords = transform_StoR(self.S_mol1, self.box)
        # Transform functions are tested elsewhere so taken as working here
        dists = self_distance_array(R_coords, box = self.box)
        # Manually calculate self_distance_array
        manual = np.zeros(len(dists), dtype = np.float64)
        distpos = 0
        for i, Ri in enumerate(R_coords):
            for Rj in R_coords[i+1:]:
                Rij = Rj - Ri
                Rij -= round(Rij[2]/self.box[2][2])*self.box[2]
                Rij -= round(Rij[1]/self.box[1][1])*self.box[1]
                Rij -= round(Rij[0]/self.box[0][0])*self.box[0]
                Rij = np.linalg.norm(Rij) # find norm of Rij vector
                manual[distpos] = Rij # and done, phew
                distpos += 1

        assert_almost_equal(dists, manual, self.prec, 
                            err_msg="self_distance_array failed with input 1")

        # Do it again for input 2 (has wider separation in points)
        # Also use boxV here in self_dist calculation
        R_coords = transform_StoR(self.S_mol2, self.box)
        # Transform functions are tested elsewhere so taken as working here
        dists = self_distance_array(R_coords, box = self.boxV)
        # Manually calculate self_distance_array
        manual = np.zeros(len(dists), dtype = np.float64)
        distpos = 0
        for i, Ri in enumerate(R_coords):
            for Rj in R_coords[i+1:]:
                Rij = Rj - Ri
                Rij -= round(Rij[2]/self.box[2][2])*self.box[2]
                Rij -= round(Rij[1]/self.box[1][1])*self.box[1]
                Rij -= round(Rij[0]/self.box[0][0])*self.box[0]
                Rij = np.linalg.norm(Rij) # find norm of Rij vector
                manual[distpos] = Rij # and done, phew
                distpos += 1
                
        assert_almost_equal(dists, manual, self.prec, 
                            err_msg="self_distance_array failed with input 2")
                
    def test_distarray(self):
        from MDAnalysis.core.distances import distance_array
        from MDAnalysis.core.distances import transform_StoR, transform_RtoS

        R_mol1 = transform_StoR(self.S_mol1, self.box)
        R_mol2 = transform_StoR(self.S_mol2, self.box)

        # Try with box
        dists = distance_array(R_mol1, R_mol2, box=self.box)
        # Manually calculate distance_array
        manual = np.zeros((len(R_mol1),len(R_mol2)))
        for i, Ri in enumerate(R_mol1):
            for j, Rj in enumerate(R_mol2):
                Rij = Rj - Ri
                Rij -= round(Rij[2]/self.box[2][2])*self.box[2]
                Rij -= round(Rij[1]/self.box[1][1])*self.box[1]
                Rij -= round(Rij[0]/self.box[0][0])*self.box[0]
                Rij = np.linalg.norm(Rij) # find norm of Rij vector
                manual[i][j] = Rij

        assert_almost_equal(dists, manual, self.prec,
                            err_msg="distance_array failed with box")
    
        # Now check using boxV
        dists = distance_array(R_mol1, R_mol2, box=self.boxV)
        assert_almost_equal(dists, manual, self.prec,
                            err_msg="distance_array failed with boxV")

    def test_pbc_dist(self):
        from MDAnalysis.core.distances import distance_array
        results = np.array([[37.629944]])

        dists = distance_array(self.S_mol1, self.S_mol2, box=self.boxV)

        assert_almost_equal(dists, results, self.prec,
                            err_msg="distance_array failed to retrieve PBC distance")
