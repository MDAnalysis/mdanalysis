import MDAnalysis
from MDAnalysis.tests.datafiles import PSF,DCD
import MDAnalysis.core.AtomGroup

from numpy.testing import *
from numpy import array, float32
from nose.plugins.attrib import attr

import tempfile

class TestAtomGroup(TestCase):
    """Tests of AtomGroup; selections are tested separately."""
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.ag = self.universe.atoms  # prototypical AtomGroup

    def test_newAtomGroup(self):
        newag = MDAnalysis.core.AtomGroup.AtomGroup(self.ag[1000:2000:200])
        assert_equal(type(newag), type(self.ag), "Failed to make a new AtomGroup: type mismatch")
        assert_equal(newag.numberOfAtoms(), len(self.ag[1000:2000:200]))
        assert_equal(newag.numberOfResidues(), 5)
        assert_almost_equal(newag.totalMass(),  40.044999999999995) # check any special method

    def test_numberOfAtoms(self):
        assert_equal(self.ag.numberOfAtoms(), 3341)

    def test_numberOfResidues(self):
        assert_equal(self.ag.numberOfResidues(), 214)

    def test_len(self):
        """testing that len(atomgroup) == atomgroup.numberOfAtoms()"""
        assert_equal(len(self.ag), self.ag.numberOfAtoms(), "len and numberOfAtoms() disagree")

    def test_centerOfGeometry(self):
        assert_array_almost_equal(self.ag.centerOfGeometry(), 
                                  array([-0.04223963,  0.0141824 , -0.03505163], dtype=float32))
    def test_centerOfMass(self):
        assert_array_almost_equal(self.ag.centerOfMass(), 
                                  array([-0.01094035,  0.05727601, -0.12885778]))

    def test_charges(self):
        assert_array_almost_equal(self.ag.charges()[1000:2000:200], 
                                  array([-0.09,  0.09, -0.47,  0.51,  0.09]))
        
    def test_coordinates(self):
        assert_array_almost_equal(self.ag.coordinates()[1000:2000:200], 
                                  array([[  3.94543672, -12.4060812 ,  -7.26820087],
                                         [ 13.21632767,   5.879035  , -14.67914867],
                                         [ 12.07735443,  -9.00604534,   4.09301519],
                                         [ 11.35541916,   7.0690732 ,  -0.32511973],
                                         [-13.26763439,   4.90658951,  10.6880455 ]], dtype=float32))
    def test_indices(self):
        assert_array_equal(self.ag.indices()[:5], array([0, 1, 2, 3, 4]))

    def test_principalAxes(self):
        assert_array_almost_equal(self.ag.principalAxes(),
                                  array([[ -9.99925632e-01,   1.21546132e-02,   9.98264877e-04],
                                         [  1.20986911e-02,   9.98951474e-01,  -4.41539838e-02],
                                         [  1.53389276e-03,   4.41386224e-02,   9.99024239e-01]]))

    def test_totalCharge(self):
        assert_almost_equal(self.ag.totalCharge(), -4.0)

    # TODO: add all other methods except selectAtoms(), see test_selections.py

    # add new methods here...


class _WriteAtoms(TestCase):
    """Set up the standard AdK system in implicit solvent."""
    ext = None   # override to test various output writers
    precision = 3

    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        suffix = '.' + self.ext
        fd, self.outfile = tempfile.mkstemp(suffix=suffix)
        
    def universe_from_tmp(self):
        return MDAnalysis.Universe(self.outfile)

    def test_write_atoms(self):
        self.universe.atoms.write(self.outfile)
        u2 = self.universe_from_tmp()
        assert_array_almost_equal(self.universe.atoms.coordinates(), u2.atoms.coordinates(), self.precision,
                                  err_msg="atom coordinate mismatch between original and %s file" % self.ext)

    def test_write_selection(self):
        CA = self.universe.selectAtoms('name CA')
        CA.write(self.outfile)
        u2 = self.universe_from_tmp()
        CA2 = u2.selectAtoms('all')   # check EVERYTHING, otherwise we might get false positives!

        assert_equal(len(u2.atoms), len(CA), "written CA selection does not match original selection")
        assert_almost_equal(CA2.coordinates(), CA.coordinates(), self.precision,
                            err_msg="CA coordinates do not agree with original")

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass

class TestWritePDB(_WriteAtoms):
    ext = "pdb"
    precision = 3

import MDAnalysis.coordinates
class TestWriteCRD(_WriteAtoms):
    ext = "crd"
    precision = 2

    # replicated methods here so that I can decorate as known failures;
    # remove once Issue 40 is resolved.
    @dec.knownfailureif(not 'CRD' in MDAnalysis.coordinates._trajectory_readers, 
                        "CRD reader is not implemented yet (see Issue 40)")
    def test_write_atoms(self):
        super(TestWriteCRD, self).test_write_atoms()

    @dec.knownfailureif(not 'CRD' in MDAnalysis.coordinates._trajectory_readers, 
                        "CRD reader is not implemented yet (see Issue 40)")
    def test_write_selection(self):
        super(TestWriteCRD, self).test_write_selection()

class TestWriteGRO(_WriteAtoms):
    ext = "gro"
    precision = 2


