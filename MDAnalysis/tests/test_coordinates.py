import MDAnalysis as mda

import numpy as np
from numpy.testing import *
from pkg_resources import resource_filename

PSF = resource_filename(__name__, 'data/adk.psf')
DCD = resource_filename(__name__, 'data/adk_dims.dcd')

PDB_small = resource_filename(__name__, 'data/adk_open.pdb')

PDB = resource_filename(__name__, 'data/adk_oplsaa.pdb')
XTC = resource_filename(__name__, 'data/adk_oplsaa.xtc')
TRR = resource_filename(__name__, 'data/adk_oplsaa.trr')

class TestPDBReader(TestCase):
    def test_load_pdb_small(self):
        U = mda.Universe(PDB_small)
        assert_equal(len(U.atoms), 3341, "load Universe from small PDB")
        assert_equal(U.atoms.selectAtoms('resid 150 and name HA2').atoms[0], 
                     U.atoms[2314], "Atom selections")

    def test_load_pdb_big(self):
        U = mda.Universe(PDB)
        assert_equal(len(U.atoms), 47681, "load Universe from big PDB")
        assert_equal(U.atoms.selectAtoms('resid 150 and name HA2').atoms[0], 
                     U.atoms[2314], "Atom selections")
        na = U.selectAtoms('resname NA+')
        assert_equal(len(na), 4, "Atom selection of last atoms in file")


class _TestDCD(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        self.dcd = self.universe.trajectory
        self.ts = self.universe.coord

class TestDCDReader(_TestDCD):
    def test_rewind_dcd(self):
        self.dcd.rewind()
        assert_equal(self.ts.frame, 1, "rewinding to frame 1")

    def test_next_dcd(self):
        self.dcd.rewind()
        self.dcd.next()
        assert_equal(self.ts.frame, 2, "loading frame 2")

    def test_jump_dcd(self):
        self.dcd[15]  # index is 0-based but frames are 1-based
        assert_equal(self.ts.frame, 16, "jumping to frame 16")

class TestDCDCorrel(_TestDCD):
    def setUp(self):
        super(TestDCDCorrel, self).setUp()
        import MDAnalysis.core.Timeseries as TS
        self.collection = TS.TimeseriesCollection()
        C = self.collection
        all = self.universe.atoms
        ca = self.universe.s4AKE.CA
        ca_termini =  mda.core.AtomGroup.AtomGroup([ca[0], ca[-1]])
        phi151 = self.universe.selectAtoms('resid 151').selectAtoms('name HN', 'name N', 'name CA', 'name CB')
        C.addTimeseries(TS.Atom('v', ca_termini))       # 0
        C.addTimeseries(TS.Bond(ca_termini))            # 1
        C.addTimeseries(TS.Bond([ca[0], ca[-1]]))       # 2
        C.addTimeseries(TS.Angle(phi151[1:4]))          # 3
        C.addTimeseries(TS.Dihedral(phi151))            # 4
        C.addTimeseries(TS.Distance('r', ca_termini))   # 5
        C.addTimeseries(TS.CenterOfMass(ca))            # 6
        C.addTimeseries(TS.CenterOfGeometry(ca))        # 7
        C.addTimeseries(TS.CenterOfMass(all))           # 8
        C.addTimeseries(TS.CenterOfGeometry(all))       # 9
        # cannot test WaterDipole because there's no water in the test dcd

    def test_correl(self):
        C = self.collection
        C.compute(self.dcd)
        avg_com_ca  = np.array([ 0.0043688 , -0.27812258, 0.0284051])
        avg_com_all = np.array([-0.10086529, -0.16357276, 0.12724672])
        avg_cog_all = np.array([-0.13554797, -0.20521885, 0.2118998])
        avg_angle = 1.9111695972912988
        avg_phi151 = 0.0088003870749735619
        avg_dist = 9.7960210987736236

        assert_equal(len(C), 10, "Correl: len(collection)")

        assert_equal(C[0].shape, (2, 3, 98), 
                     "Correl: Atom positions")
        assert_array_equal(C[1].__data__, C[2].__data__, 
                           "Correl: Bonds with lists and AtomGroup")
        assert_almost_equal(C[3].__data__.mean(), avg_angle,
                            err_msg="Correl: average Angle")
        assert_almost_equal(C[4].__data__.mean(), avg_phi151,
                            err_msg="Correl: average Dihedral")
        assert_almost_equal(C[5].__data__.mean(), avg_dist,
                            err_msg="Correl: average scalar Distance")
        assert_array_almost_equal(C[6].__data__.mean(axis=1), avg_com_ca,
                                  err_msg="Correl: average CA CentreOfMass")
        assert_array_almost_equal(C[6].__data__, C[7].__data__,
                                  err_msg="Correl: CA CentreOfMass == CentreOfGeometry")
        assert_almost_equal(C[8].__data__.mean(axis=1), avg_com_all,
                            err_msg="Correl: average all CenterOfMass")
        assert_almost_equal(C[9].__data__.mean(axis=1), avg_cog_all,
                            err_msg="Correl: average all CenterOfGeometry")
        
        C.clear()
        assert_equal(len(C), 0, "Correl: clear()")
   

def compute_correl_references():
    universe = MDAnalysis.Universe(PSF,DCD)

    all = universe.atoms
    ca = universe.s4AKE.CA
    ca_termini =  mda.core.AtomGroup.AtomGroup([ca[0], ca[-1]])
    phi151 = universe.selectAtoms('resid 151').selectAtoms('name HN', 'name N', 'name CA', 'name CB')

    C = MDAnalysis.collection
    C.clear()

    C.addTimeseries(TS.Atom('v', ca_termini))       # 0
    C.addTimeseries(TS.Bond(ca_termini))            # 1
    C.addTimeseries(TS.Bond([ca[0], ca[-1]]))       # 2
    C.addTimeseries(TS.Angle(phi151[1:4]))          # 3
    C.addTimeseries(TS.Dihedral(phi151))            # 4
    C.addTimeseries(TS.Distance('r', ca_termini))   # 5
    C.addTimeseries(TS.CenterOfMass(ca))            # 6
    C.addTimeseries(TS.CenterOfGeometry(ca))        # 7
    C.addTimeseries(TS.CenterOfMass(all))           # 8
    C.addTimeseries(TS.CenterOfGeometry(all))       # 9

    C.compute(universe.dcd)

    results = {"avg_angle": C[3].__data__.mean(),
               "avg_phi151": C[4].__data__.mean(),
               "avg_dist": C[5].__data__.mean(),
               "avg_com_ca": C[6].__data__.mean(axis=1), 
               "avg_com_all": C[8].__data__.mean(axis=1),
               "avg_cog_all": C[9].__data__.mean(axis=1),
               }
    C.clear()
    return results


class _TestGromacsReader(TestCase):
    def test_rewind_xdrtrj(self):
        self.trajectory.rewind()
        assert_equal(self.ts.frame, 1, "rewinding to frame 1")

    def test_next_xdrtrj(self):
        self.trajectory.rewind()
        self.trajectory.next()
        assert_equal(self.ts.frame, 2, "loading frame 2")

class TestXTCReader(_TestGromacsReader):
    def setUp(self):
        self.universe = mda.Universe(PDB, XTC)
        self.trajectory = self.universe.trajectory
        self.ts = self.universe.coord

class TestTRRReader(_TestGromacsReader):
    def setUp(self):
        self.universe = mda.Universe(PDB, TRR)
        self.trajectory = self.universe.trajectory
        self.ts = self.universe.coord

