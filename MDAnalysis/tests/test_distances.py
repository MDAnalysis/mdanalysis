import MDAnalysis

import numpy as np
from numpy.testing import *
from nose.plugins.attrib import attr

from MDAnalysis.tests.datafiles import PSF,DCD,PDB_small,PDB,XTC,TRR

class TestDistanceArrayDCD(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.dcd = self.universe.trajectory
        self.ca = self.universe.selectAtoms('name CA')

    @attr('issue')    
    def test_simple(self):
        U = self.universe
        self.dcd.rewind()
        x0 = U.atoms.coordinates(copy=True)
        self.dcd[10]
        x1 = U.atoms.coordinates(copy=True)
        d = MDAnalysis.core.distances.distance_array(x0, x1)
        assert_equal(d.shape, (3341, 3341), "wrong shape (should be (Natoms,Natoms))")
        assert_almost_equal(d.min(), 0.11981228170520701, err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 53.572192429459619,  err_msg="wrong maximum distance value")

    @attr('issue')
    def test_copyTrue(self):
        U = self.universe
        self.dcd.rewind()
        x0 = U.atoms.coordinates(copy=True)
        self.dcd[10]
        x1 = U.atoms.coordinates(copy=True)
        d = MDAnalysis.core.distances.distance_array(x0, x1, copy=True)
        assert_equal(d.shape, (3341, 3341), "wrong shape, should be (Natoms,Natoms)")
        assert_almost_equal(d.min(), 0.11981228170520701, err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 53.572192429459619,  err_msg="wrong maximum distance value")

    @dec.knownfailureif(True, "Using copy=False in a distance array is a known failure.")
    @attr('issue')
    def test_copyFalse(self):
        # This will fail at the moment until the underlying C-code is improved.
        # Should give the same answer as copy=True but does not, see Issue 4
        U = self.universe
        self.dcd.rewind()
        x0 = U.atoms.coordinates(copy=True)
        self.dcd[10]
        x1 = U.atoms.coordinates(copy=True)
        d = MDAnalysis.core.distances.distance_array(x0, x1, copy=False)
        assert_equal(d.shape, (3341, 3341), "wrong shape, should be (Natoms,Natoms) entries")
        # should give the same answer as copy=True but does not, see Issue 4
        assert_almost_equal(d.min(), 0.11981228170520701, 
                            err_msg="wrong minimum distance value with copy=False, see Issue 4")
        assert_almost_equal(d.max(), 53.572192429459619,  
                            err_msg="wrong maximum distance value with copy=False, see Issue 4")

    @dec.knownfailureif(True, "Using copy=False in a distance array is a known failure.")
    @attr('issue')
    def test_Issue4(self):
        U = self.universe
        self.dcd.rewind()
        x0 = self.ca.coordinates(copy=True)
        self.dcd[10]
        x1 = self.ca.coordinates(copy=True)
        d_ok = MDAnalysis.core.distances.distance_array(x0, x1, copy=True)
        d_wrong = MDAnalysis.core.distances.distance_array(x0, x1, copy=False)
        assert_equal(d_ok.shape, (len(self.ca), len(self.ca)), "copy=True square matrix with N_CA entries")
        assert_equal(d_wrong.shape, (len(self.ca), len(self.ca)), "copy=False square matrix with N_CA entries")
        # should give the same answer as copy=True but does not, see Issue 4
        assert_almost_equal(d_ok.min(), d_wrong.min(), 
                            err_msg="copy=False distances != copy=True distances (tested min), known bug")
        assert_almost_equal(d_ok.max(), d_wrong.max(), 
                            err_msg="copy=False distances != copy=True distances (tested max), known bug")

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
        assert_almost_equal(d.min(), 0.11981228170520701, err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 53.572192429459619,  err_msg="wrong maximum distance value")

    def test_periodic(self):
        # boring with the current dcd as that has no PBC
        U = self.universe
        self.dcd.rewind()
        x0 = U.atoms.coordinates(copy=True)
        self.dcd[10]
        x1 = U.atoms.coordinates(copy=True)
        d = MDAnalysis.core.distances.distance_array(x0, x1, box=U.coord.dimensions)
        assert_equal(d.shape, (3341, 3341), "should be square matrix with Natoms entries")
        assert_almost_equal(d.min(), 0.11981228170520701, err_msg="wrong minimum distance value with PBC")
        assert_almost_equal(d.max(), 53.572192429459619,  err_msg="wrong maximum distance value with PBC")

        
class TestSelfDistanceArrayDCD(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.dcd = self.universe.trajectory
        self.ca = self.universe.selectAtoms('name CA')
    
    def test_simple(self):
        U = self.universe
        self.dcd.rewind()
        x0 = U.atoms.coordinates(copy=True)
        d = MDAnalysis.core.distances.self_distance_array(x0)
        N = 3341 * (3341 - 1) / 2
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 52.4702570624190590, err_msg="wrong maximum distance value")

    def test_copyTrue(self):
        U = self.universe
        self.dcd.rewind()
        x0 = U.atoms.coordinates(copy=True)
        d = MDAnalysis.core.distances.self_distance_array(x0, copy=True)
        N = 3341 * (3341 - 1) / 2
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 52.4702570624190590, err_msg="wrong maximum distance value")

    @dec.knownfailureif(True, "Using copy=False in a distance array is a known failure.")
    @attr('issue')
    def test_copyFalse(self):
        # This will fail at the moment until the underlying C-code is improved.
        # Should give the same answer as copy=True but does not, see Issue 4
        U = self.universe
        self.dcd.rewind()
        x0 = U.atoms.coordinates(copy=True)
        d = MDAnalysis.core.distances.self_distance_array(x0, copy=False)
        N = 3341 * (3341 - 1) / 2
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        # should give the same answer as copy=True but does not, see Issue 4
        assert_almost_equal(d.min(), 0.92905562402529318, 
                            err_msg="wrong minimum distance value with copy=False, see Issue 4")
        assert_almost_equal(d.max(), 52.4702570624190590,  
                            err_msg="wrong maximum distance value with copy=False, see Issue 4")

    @dec.knownfailureif(True, "Using copy=False in a distance array is a known failure.")
    @attr('issue')
    def test_Issue4(self):
        U = self.universe
        self.dcd.rewind()
        x0 = self.ca.coordinates(copy=True)
        d_ok = MDAnalysis.core.distances.self_distance_array(x0, copy=True)
        d_wrong = MDAnalysis.core.distances.self_distance_array(x0, copy=False)
        N = len(self.ca) * (len(self.ca) - 1) / 2
        assert_equal(d_ok.shape, (N,), "copy=True 1D array with N_CA entries")
        assert_equal(d_wrong.shape, (N,), "copy=False 1D array with N_CA entries")
        # should give the same answer as copy=True but does not, see Issue 4
        assert_almost_equal(d_ok.min(), d_wrong.min(), 
                            err_msg="copy=False distances != copy=True distances (tested min), known bug")
        assert_almost_equal(d_ok.max(), d_wrong.max(), 
                            err_msg="copy=False distances != copy=True distances (tested max), known bug")


    def test_outarray(self):
        U = self.universe
        self.dcd.rewind()
        x0 = U.atoms.coordinates(copy=True)
        natoms = len(U.atoms)
        N = natoms*(natoms-1) / 2
        d = np.zeros((N,), np.float64)
        MDAnalysis.core.distances.self_distance_array(x0, result=d)
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, err_msg="wrong minimum distance value")
        assert_almost_equal(d.max(), 52.4702570624190590, err_msg="wrong maximum distance value")

    def test_periodic(self):
        # boring with the current dcd as that has no PBC
        U = self.universe
        self.dcd.rewind()
        x0 = U.atoms.coordinates(copy=True)
        natoms = len(U.atoms)
        N = natoms*(natoms-1) / 2
        d = MDAnalysis.core.distances.self_distance_array(x0, box=U.coord.dimensions)
        assert_equal(d.shape, (N,), "wrong shape (should be (Natoms*(Natoms-1)/2,))")
        assert_almost_equal(d.min(), 0.92905562402529318, err_msg="wrong minimum distance value with PBC")
        assert_almost_equal(d.max(), 52.4702570624190590, err_msg="wrong maximum distance value with PBC")

        
