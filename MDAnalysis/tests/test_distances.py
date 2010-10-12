import MDAnalysis
import MDAnalysis.core.distances

import numpy as np
from numpy.testing import *
from nose.plugins.attrib import attr

from MDAnalysis.tests.datafiles import PSF,DCD

class TestDistanceArrayDCD(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.dcd = self.universe.trajectory
        self.ca = self.universe.selectAtoms('name CA')
        # reasonable precision so that tests succeed on 32 and 64 bit machines
        # (the reference values were obtained on 64 bit)
        # Example:
        #   Items are not equal: wrong maximum distance value
        #   ACTUAL: 52.470254967456412
        #   DESIRED: 52.470257062419059
        self.prec = 5

    @attr('issue')    
    def test_simple(self):
        U = self.universe
        self.dcd.rewind()
        x0 = U.atoms.coordinates(copy=True)
        self.dcd[10]
        x1 = U.atoms.coordinates(copy=True)
        d = MDAnalysis.core.distances.distance_array(x0, x1)
        assert_equal(d.shape, (3341, 3341), "wrong shape (should be (Natoms,Natoms))")
        assert_almost_equal(d.min(), 0.11981228170520701, self.prec, 
                            err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 53.572192429459619, self.prec,
                            err_msg="wrong maximum distance value")

    def test_outarray(self):
        U = self.universe
        self.dcd.rewind()
        x0 = U.atoms.coordinates(copy=True)
        self.dcd[10]
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
        self.dcd.rewind()
        x0 = U.atoms.coordinates(copy=True)
        self.dcd[10]
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
        self.dcd = self.universe.trajectory
        self.ca = self.universe.selectAtoms('name CA')
        # see comments above on precision
        self.prec = 5
    
    def test_simple(self):
        U = self.universe
        self.dcd.rewind()
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
        self.dcd.rewind()
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
        self.dcd.rewind()
        x0 = U.atoms.coordinates(copy=True)
        natoms = len(U.atoms)
        N = natoms*(natoms-1) / 2
        d = MDAnalysis.core.distances.self_distance_array(x0, box=U.coord.dimensions)
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, self.prec,
                            err_msg="wrong minimum distance value with PBC")
        assert_almost_equal(d.max(), 52.4702570624190590, self.prec,
                            err_msg="wrong maximum distance value with PBC")

        
