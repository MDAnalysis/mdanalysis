import MDAnalysis
import MDAnalysis.analysis.distances as distances

import numpy as np
from numpy.testing import *
from nose.plugins.attrib import attr

from MDAnalysis.tests.datafiles import PSF,DCD

class TestContactMatrix(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.dcd = self.universe.trajectory
        # reasonable precision so that tests succeed on 32 and 64 bit machines
        # (the reference values were obtained on 64 bit)
        # Example:
        #   Items are not equal: wrong maximum distance value
        #   ACTUAL: 52.470254967456412
        #   DESIRED: 52.470257062419059
        self.prec = 5

    def test_numpy(self):
        U = self.universe
        self.dcd.rewind()
        self.dcd[10]
        # small cutoff value as the input file is a protein
        contacts = distances.contact_matrix(U.atoms.coordinates() , cutoff=1.5 , returntype="numpy")
        assert_equal(contacts.shape, (3341, 3341), "wrong shape (should be (Natoms,Natoms))")
        assert_equal(contacts[0][0] , True , "first entry should be a contact")
        assert_equal(contacts[0][-1] , False , "last entry for first atom should be a non-contact")

    def test_sparse(self):
        U = self.universe
        self.dcd.rewind()
        self.dcd[10]
        # Just taking first 50 atoms as the sparse method is slow
        selection = U.selectAtoms('bynum 1:50')
        # small cutoff value as the input file is a protein
        # High progress_meter_freq so progress meter is not printed during test
        contacts = distances.contact_matrix(selection.coordinates() , cutoff=1.0 , returntype="sparse" , suppress_progmet=True)
        assert_equal(contacts.shape, (50, 50), "wrong shape (should be (50,50))")
        assert_equal(contacts[0,0] , False , "entry (0,0) should be a non-contact")
        assert_equal(contacts[19,20] , True , "entry (19,20) should be a contact")


