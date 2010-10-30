import MDAnalysis
import MDAnalysis as mda
import MDAnalysis.coordinates

import numpy as np
from numpy.testing import *
from nose.plugins.attrib import attr

from MDAnalysis.tests.datafiles import PSF,DCD,DCD_empty,PDB_small,PDB,CRD,XTC,TRR,GRO,XYZ,XYZ_bz2,XYZ_psf

import os.path, tempfile
import itertools

def atom_distance(a, b):
    """Calculate the distance between two atoms a and b."""
    r = a.pos - b.pos
    return np.sqrt(np.sum(r**2))

class RefAdKSmall(object):
    """Mixin class to provide comparison numbers.

    Based on small PDB with AdK (:data:`PDB_small`).

    .. Note:: All distances must be in ANGSTROEM as this is the
       MDAnalysis default unit. All readers must return Angstroem by
       default.
    """
    ref_coordinates = {
        'A10CA': np.array([ -1.198, 7.937, 22.654]),   # G11:CA, copied frm adk_open.pdb
        }
    ref_distances = {'endtoend': 11.016959}
    ref_E151HA2_index = 2314
    ref_numatoms = 3341

class RefAdK(object):
    """Mixin class to provide comparison numbers.

    Based on PDB/GRO with AdK in water + Na+ (:data:`PDB`).

    .. Note:: All distances must be in ANGSTROEM as this is the
       MDAnalysis default unit. All readers must return Angstroem by
       default.
    """
    ref_coordinates = {
        'A10CA': np.array([ 62.97600174,  62.08800125,  20.2329998 ]),  # Angstroem as MDAnalysis unit!!
        }
    ref_distances = {'endtoend': 9.3513174}
    ref_E151HA2_index = 2314
    ref_numatoms = 47681
    ref_Na_sel_size = 4
    ref_unitcell = np.array([ 80.017,  80.017,  80.017,  90., 60., 60.], dtype=np.float32)

class Ref2r9r(object):
    """Mixin class to provide comparison numbers.

    Based on S6 helices of chimeric Kv channel 

    .. Note:: All distances must be in ANGSTROEM as this is the
       MDAnalysis default unit. All readers must return Angstroem by
       default.
    """
    ref_numatoms = 1284
    ref_sum_centre_of_geometry = -98.24146
    ref_numframes = 10

class TestXYZReader(TestCase, Ref2r9r):
    def setUp(self):
        self.universe = mda.Universe(XYZ_psf, XYZ)
        self.prec = 3 # 4 decimals in xyz file

    def tearDown(self):
        del self.universe

    def test_load_xyz(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_numatoms, "load Universe from PSF and XYZ")
   
    def test_numatoms(self):
        assert_equal(self.universe.trajectory.numatoms, self.ref_numatoms, "wrong number of atoms")

    def test_numframes(self):
        assert_equal(self.universe.trajectory.numframes, self.ref_numframes, "wrong number of frames in xyz")

    def test_sum_centres_of_geometry(self):
        centreOfGeometry=0
        
        for i in self.universe.trajectory:
            sel = self.universe.selectAtoms("all")      
            centreOfGeometry+=sum(sel.centerOfGeometry())  
        
        assert_almost_equal(centreOfGeometry, self.ref_sum_centre_of_geometry, self.prec,
                            err_msg="sum of centers of geometry over the trajectory do not match")


class TestCompressedXYZReader(TestCase, Ref2r9r):
    def setUp(self):
        self.universe = mda.Universe(XYZ_psf, XYZ_bz2)
        self.prec = 3 # 4 decimals in xyz file

    def tearDown(self):
        del self.universe

    def test_load_xyz(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_numatoms, "load Universe from PSF and XYZ")
   
    def test_numatoms(self):
        assert_equal(self.universe.trajectory.numatoms, self.ref_numatoms, "wrong number of atoms")

    def test_numframes(self):
        assert_equal(self.universe.trajectory.numframes, self.ref_numframes, "wrong number of frames in xyz")

    def test_sum_centres_of_geometry(self):
        centreOfGeometry=0
        
        for i in self.universe.trajectory:
            sel = self.universe.selectAtoms("all")      
            centreOfGeometry+=sum(sel.centerOfGeometry())  
        
        assert_almost_equal(centreOfGeometry, self.ref_sum_centre_of_geometry, self.prec,
                            err_msg="sum of centers of geometry over the trajectory do not match")



class _SingleFrameReader(TestCase, RefAdKSmall):
    # see TestPDBReader how to set up!

    def tearDown(self):
        del self.universe

    def test_load_file(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_numatoms, "load Universe from file %s" % U.trajectory.filename)
        assert_equal(U.atoms.selectAtoms('resid 150 and name HA2').atoms[0], 
                     U.atoms[self.ref_E151HA2_index], "Atom selections")

    def test_numatoms(self):
        assert_equal(self.universe.trajectory.numatoms, self.ref_numatoms, "wrong number of atoms")

    def test_numres(self):
        assert_equal(self.universe.atoms.numberOfResidues(), 214, "wrong number of residues")

    def test_numframes(self):
        assert_equal(self.universe.trajectory.numframes, 1, "wrong number of frames in pdb")

    def test_coordinates(self):
        A10CA = self.universe.s4AKE.CA[10]
        # restrict accuracy to maximum in PDB files (3 decimals)
        assert_almost_equal(A10CA.pos, self.ref_coordinates['A10CA'], 3,
                            err_msg="wrong coordinates for A10:CA")
        
    def test_distances(self):
        NTERM = self.universe.s4AKE.N[0]
        CTERM = self.universe.s4AKE.C[-1]
        d = atom_distance(NTERM, CTERM)
        assert_almost_equal(d, self.ref_distances['endtoend'], self.prec,
                            err_msg="distance between M1:N and G214:C")

class TestPDBReader(_SingleFrameReader):
    def setUp(self):
        self.universe = mda.Universe(PDB_small) 
        self.prec = 3  # 3 decimals in PDB spec http://www.wwpdb.org/documentation/format32/sect9.html#ATOM

class TestPSF_CRDReader(_SingleFrameReader):
    def setUp(self):
        self.universe = mda.Universe(PSF, CRD)
        self.prec = 5  # precision in CRD (at least we are writing %9.5f)

class TestPSF_PDBReader(TestPDBReader):
    def setUp(self):
        self.universe = mda.Universe(PSF, PDB_small)
        self.prec = 3  # 3 decimals in PDB spec http://www.wwpdb.org/documentation/format32/sect9.html#ATOM

class TestPrimitivePDBReader(TestPDBReader):
    def setUp(self):
        self.universe = mda.Universe(PDB_small, permissive=True) 
        self.prec = 3  # 3 decimals in PDB spec http://www.wwpdb.org/documentation/format32/sect9.html#ATOM

class TestPSF_PrimitivePDBReader(TestPDBReader):
    def setUp(self):
        self.universe = mda.Universe(PSF, PDB_small, permissive=True) 
        self.prec = 3  # 3 decimals in PDB spec http://www.wwpdb.org/documentation/format32/sect9.html#ATOM

class TestGROReader(TestCase, RefAdK):
    def setUp(self):
	self.universe = mda.Universe(GRO)
        self.ts = self.universe.trajectory.ts
        self.prec = 2  # lower prec in gro!! (3 decimals nm -> 2 decimals in Angstroem)

    def tearDown(self):
        del self.universe
        del self.ts

    def test_load_gro(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_numatoms, "load Universe from small GRO")
        assert_equal(U.atoms.selectAtoms('resid 150 and name HA2').atoms[0], 
                     U.atoms[self.ref_E151HA2_index], "Atom selections")

    def test_numatoms(self):
        assert_equal(self.universe.trajectory.numatoms, self.ref_numatoms, "wrong number of atoms")

    def test_numframes(self):
        assert_equal(self.universe.trajectory.numframes, 1, "wrong number of frames")

    def test_coordinates(self):
        A10CA = self.universe.SYSTEM.CA[10]
        assert_almost_equal(A10CA.pos, self.ref_coordinates['A10CA'], self.prec,
                            err_msg="wrong coordinates for A10:CA")
        
    def test_distances(self):
        # NOTE that the prec is only 1 decimal: subtracting two low precision coordinates
        #      low prec: 9.3455122920041109; high prec (from pdb): 9.3513174
        NTERM = self.universe.SYSTEM.N[0]
        CTERM = self.universe.SYSTEM.C[-1]
        d = atom_distance(NTERM, CTERM)
        assert_almost_equal(d, self.ref_distances['endtoend'], self.prec - 1,  # note low prec!! 
                            err_msg="distance between M1:N and G214:C")

    def test_selection(self):
        na = self.universe.selectAtoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size, "Atom selection of last atoms in file")

    def test_unitcell(self):
        assert_array_almost_equal(self.ts.dimensions, self.ref_unitcell, self.prec, 
                                  err_msg="unit cell dimensions (rhombic dodecahedron)")

class TestGROReaderNoConversion(TestCase, RefAdK):
    def setUp(self):
        mda.core.flags['convert_gromacs_lengths'] = False
        self.universe = mda.Universe(GRO)
        self.ts = self.universe.trajectory.ts
        self.prec = 3

    def tearDown(self):
        mda.core.flags['convert_gromacs_lengths'] = True  # default
        del self.universe
        del self.ts

    def test_coordinates(self):
        # note: these are the native coordinates in nm; for the test to succeed:
        assert_equal(mda.core.flags['convert_gromacs_lengths'], False, 
                     "oops, mda.core.flags['convert_gromacs_lengths'] should be False for this test")
        A10CA = self.universe.SYSTEM.CA[10]
        assert_almost_equal(A10CA.pos, RefAdK.ref_coordinates['A10CA']/10.0,  # coordinates in nm 
                            self.prec,
                            err_msg="wrong native coordinates (in nm) for A10:CA")
        
    def test_distances(self):
        # note: these are the native coordinates in nm; for the test to succeed:
        assert_equal(mda.core.flags['convert_gromacs_lengths'], False, 
                     "oops, mda.core.flags['convert_gromacs_lengths'] should be False for this test")

        # 3 decimals on nm in gro but we compare to the distance
        # computed from the pdb file, so the effective precision is 2 again.
        # (Otherwise the distance test fails: 
        #  Arrays are not almost equal distance between M1:N and G214:C
        #    ACTUAL: 0.93455122920041123
        #    DESIRED: 0.93513173999999988
        NTERM = self.universe.SYSTEM.N[0]
        CTERM = self.universe.SYSTEM.C[-1]
        d = atom_distance(NTERM, CTERM)
        assert_almost_equal(d, RefAdK.ref_distances['endtoend']/10.0,  # coordinates in nm 
                            self.prec - 1,
                            err_msg="distance between M1:N and G214:C")

    def test_unitcell(self):
        # lengths in A : convert to nm
        assert_array_almost_equal(self.ts.dimensions[:3], self.ref_unitcell[:3]/10.0, self.prec, 
                                  err_msg="unit cell A,B,C (rhombic dodecahedron)")
        # angles should not have changed
        assert_array_almost_equal(self.ts.dimensions[3:], self.ref_unitcell[3:], self.prec, 
                                  err_msg="unit cell alpha,bet,gamma (rhombic dodecahedron)")
        

class TestPDBReaderBig(TestCase, RefAdK):
    def setUp(self):
        self.universe = mda.Universe(PDB)
        self.prec = 6

    def tearDown(self):
        del self.universe

    @dec.slow
    def test_load_pdb(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_numatoms, "load Universe from big PDB")
        assert_equal(U.atoms.selectAtoms('resid 150 and name HA2').atoms[0], 
                     U.atoms[self.ref_E151HA2_index], "Atom selections")

    @dec.slow
    def test_selection(self):
        na = self.universe.selectAtoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size, "Atom selection of last atoms in file")

    @dec.slow
    def test_numatoms(self):
        assert_equal(self.universe.trajectory.numatoms, self.ref_numatoms, "wrong number of atoms")

    @dec.slow
    def test_numframes(self):
        assert_equal(self.universe.trajectory.numframes, 1, "wrong number of frames")

    @dec.slow
    def test_coordinates(self):
        A10CA = self.universe.SYSTEM.CA[10]
        assert_almost_equal(A10CA.pos, self.ref_coordinates['A10CA'], self.prec,
                            err_msg="wrong coordinates for A10:CA")

    @dec.slow        
    def test_distances(self):
        NTERM = self.universe.SYSTEM.N[0]
        CTERM = self.universe.SYSTEM.C[-1]
        d = atom_distance(NTERM, CTERM)
        assert_almost_equal(d, self.ref_distances['endtoend'], self.prec,
                            err_msg="wrong distance between M1:N and G214:C")

    @dec.slow
    def test_selection(self):
        na = self.universe.selectAtoms('resname NA+')
        assert_equal(len(na), self.ref_Na_sel_size, "Atom selection of last atoms in file")


@attr('issue')
def TestDCD_Issue32():
    """Test for Issue 32: 0-size dcds lead to a segfault: now caught with IOError"""
    assert_raises(IOError, mda.Universe, PSF, DCD_empty)

class _TestDCD(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        self.dcd = self.universe.trajectory
        self.ts = self.universe.coord

    def tearDown(self):
        del self.universe
        del self.dcd
        del self.ts

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

    def test_jump_lastframe_dcd(self):
        self.dcd[-1]
        assert_equal(self.ts.frame, 98, "indexing last frame with dcd[-1]") 

    def test_slice_dcd(self):
        frames = [ts.frame for ts in self.dcd[5:17:3]]
        assert_equal(frames, [6, 9, 12, 15], "slicing dcd [5:17:3]")

    def test_reverse_dcd(self):
        frames = [ts.frame for ts in self.dcd[20:5:-1]]
        assert_equal(frames, range(21,6,-1), "reversing dcd [20:5:-1]")        

    def test_numatoms(self):
        assert_equal(self.universe.trajectory.numatoms, 3341, "wrong number of atoms")

    def test_numframes(self):
        assert_equal(self.universe.trajectory.numframes, 98, "wrong number of frames in dcd")

    def test_dt(self):
        assert_equal(self.universe.trajectory.dt, 1.0, "wrong timestep dt")

    def test_totaltime(self):
        assert_equal(self.universe.trajectory.totaltime, 98.0, "wrong total length of AdK trajectory")

class TestDCDWriter(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        ext = ".dcd"
        fd, self.outfile = tempfile.mkstemp(suffix=ext)
        self.Writer = MDAnalysis.coordinates.DCD.DCDWriter

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass
        del self.universe
        del self.Writer

    @attr('issue')    
    def test_write_trajectory(self):
        """Test writing DCD trajectories (Issue 50)"""
        t = self.universe.trajectory
        W = self.Writer(self.outfile, t.numatoms, delta=t.delta, step=t.skip_timestep)
        for ts in self.universe.trajectory:
            W.write_next_timestep(ts)
        W.close()

        uw = mda.Universe(PSF, self.outfile)

        # check that the coordinates are identical for each time step
        for orig_ts, written_ts in itertools.izip(self.universe.trajectory, uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, 3,
                                      err_msg="coordinate mismatch between original and written trajectory at frame %d (orig) vs %d (written)" % (orig_ts.frame, written_ts.frame))

        
class TestDCDCorrel(_TestDCD):
    def setUp(self):
        # Note: setUp is executed for *every* test !
        super(TestDCDCorrel, self).setUp()
        import MDAnalysis.core.Timeseries as TS
        self.collection = TS.TimeseriesCollection()
        C = self.collection
        all = self.universe.atoms
        ca = self.universe.s4AKE.CA
        ca_termini =  mda.core.AtomGroup.AtomGroup([ca[0], ca[-1]])
        # note that this is not quite phi... HN should be C of prec. residue
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
        C.compute(self.dcd)

    def tearDown(self):
        del self.collection
        super(TestDCDCorrel, self).tearDown()

    def test_correl(self):
        assert_equal(len(self.collection), 10, "Correl: len(collection)")

    def test_Atom(self):
        assert_equal(self.collection[0].shape, (2, 3, 98), 
                     "Correl: Atom positions")

    def test_Bonds(self):
        C = self.collection
        assert_array_equal(C[1].__data__, C[2].__data__, 
                           "Correl: Bonds with lists and AtomGroup")

    def test_Angle(self):
        C = self.collection
        avg_angle = 1.9111695972912988
        assert_almost_equal(C[3].__data__.mean(), avg_angle,
                            err_msg="Correl: average Angle")
        
    def test_Dihedral(self):
        C = self.collection
        avg_phi151 = 0.0088003870749735619        
        assert_almost_equal(C[4].__data__.mean(), avg_phi151,
                            err_msg="Correl: average Dihedral")

    def test_scalarDistance(self):
        C = self.collection        
        avg_dist = 9.7960210987736236
        assert_almost_equal(C[5].__data__.mean(), avg_dist,
                            err_msg="Correl: average scalar Distance")

    def test_CenterOfMass(self):
        C = self.collection
        avg_com_ca  = np.array([ 0.0043688 , -0.27812258, 0.0284051])
        avg_com_all = np.array([-0.10086529, -0.16357276, 0.12724672])
        assert_array_almost_equal(C[6].__data__.mean(axis=1), avg_com_ca,
                                  err_msg="Correl: average CA CenterOfMass")
        assert_almost_equal(C[8].__data__.mean(axis=1), avg_com_all,
                            err_msg="Correl: average all CenterOfMass")

    def test_CenterOfGeometry(self):
        C = self.collection        
        avg_cog_all = np.array([-0.13554797, -0.20521885, 0.2118998])
        assert_almost_equal(C[9].__data__.mean(axis=1), avg_cog_all,
                            err_msg="Correl: average all CenterOfGeometry")

    def test_CA_COMeqCOG(self):
        C = self.collection
        assert_array_almost_equal(C[6].__data__, C[7].__data__,
                                  err_msg="Correl: CA CentreOfMass == CenterOfGeometry")

    def test_clear(self):
        C = self.collection
        C.clear()
        assert_equal(len(C), 0, "Correl: clear()")
   
# notes:
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

class TestChainedReader(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, [DCD, DCD, DCD])
        self.dcd = self.universe.trajectory
        self.ts = self.universe.coord

    def test_next_dcd(self):
        self.dcd.rewind()
        self.dcd.next()
        assert_equal(self.ts.frame, 2, "loading frame 2")

    def test_numatoms(self):
        assert_equal(self.universe.trajectory.numatoms, 3341, "wrong number of atoms")

    def test_numframes(self):
        assert_equal(self.universe.trajectory.numframes, 3*98, "wrong number of frames in chained dcd")

    def test_iteration(self):
        for ts in self.dcd:
            pass # just forward to last frame
        assert_equal(self.dcd.numframes, ts.frame, 
                     "iteration yielded wrong number of frames (%d), should be %d" \
                         % (ts.frame, self.dcd.numframes))

    @dec.knownfailureif(True, "indexing not implemented for chained reader")
    def test_jump_lastframe_dcd(self):
        self.dcd[-1]
        assert_equal(self.ts.frame, self.dcd.numframes, "indexing last frame with dcd[-1]") 

    @dec.knownfailureif(True, "slicing not implemented for chained reader")
    def test_slice_dcd(self):
        frames = [ts.frame for ts in self.dcd[5:17:3]]
        assert_equal(frames, [6, 9, 12, 15], "slicing dcd [5:17:3]")

class _GromacsReader(TestCase):
    # This base class assumes same lengths and dt for XTC and TRR test cases!
    filename = None
    def setUp(self):
        self.universe = mda.Universe(GRO, self.filename) # loading from GRO is 4x faster than the PDB reader
        self.trajectory = self.universe.trajectory
        self.prec = 3
        self.ts = self.universe.coord
    
    @dec.slow
    def test_rewind_xdrtrj(self):
        self.trajectory.rewind()
        assert_equal(self.ts.frame, 1, "rewinding to frame 1")

    @dec.slow
    def test_next_xdrtrj(self):
        self.trajectory.rewind()
        self.trajectory.next()
        assert_equal(self.ts.frame, 2, "loading frame 2")

    @dec.slow
    def test_jump_xdrtrj(self):
        self.trajectory[4]  # index is 0-based but frames are 1-based
        assert_equal(self.ts.frame, 5, "jumping to frame 5")

    @dec.slow
    def test_jump_lastframe_xdrtrj(self):
        self.trajectory[-1]
        assert_equal(self.ts.frame, 10, "indexing last frame with trajectory[-1]") 

    @dec.slow
    def test_slice_xdrtrj(self):
        frames = [ts.frame for ts in self.trajectory[2:9:3]]
        assert_equal(frames,  [3, 6, 9], "slicing xdrtrj [2:9:3]")

    @dec.slow
    @dec.knownfailureif(True, "XTC/TRR reverse slicing not implemented for performance reasons")
    def test_reverse_xdrtrj(self):
        frames = [ts.frame for ts in self.trajectory[::-1]]
        assert_equal(frames, range(10,0,-1), "slicing xdrtrj [::-1]")        

    @dec.slow
    def test_coordinates(self):
        # note: these are the native coordinates in nm; for the test to succeed:
        assert_equal(mda.core.flags['convert_gromacs_lengths'], True, 
                     "oops, mda.core.flags['convert_gromacs_lengths'] should be True")
        ca_nm = np.array([[ 6.043369675,  7.385184479,  1.381425762]], dtype=np.float32)
        # coordinates in the base unit (needed for True)
        ca_Angstrom = ca_nm * 10.0
        U = self.universe
        T = U.trajectory
        T.rewind()
        T.next()
        T.next()
        assert_equal(self.ts.frame, 3, "failed to step to frame 3")
        ca = U.selectAtoms('name CA and resid 122')
        # low precision match (2 decimals in A, 3 in nm) because the above are the trr coords
        assert_array_almost_equal(ca.coordinates(), ca_Angstrom, 2,
                                  err_msg="coords of Ca of resid 122 do not match for frame 3")

    @dec.slow
    @attr('issue')
    def test_unitcell(self):
        """Test that xtc/trr unitcell is read correctly (Issue 34)"""
        self.universe.trajectory.rewind()
        uc = self.ts.dimensions
        ref_uc = np.array([ 80.017,  80.017,  80.017,  90., 60., 60.], dtype=np.float32)
        assert_array_almost_equal(uc, ref_uc, self.prec, err_msg="unit cell dimensions (rhombic dodecahedron)")

    @dec.slow
    def test_dt(self):
        assert_equal(self.universe.trajectory.dt, 100.0, "wrong timestep dt")

    @dec.slow
    def test_totaltime(self):
        assert_equal(self.universe.trajectory.totaltime, 1000.0, "wrong total length of trajectory")


class TestXTCReader(_GromacsReader):
    filename = XTC

class TestTRRReader(_GromacsReader):
    filename = TRR


class _XDRNoConversion(TestCase):
    filename = None
    def setUp(self):
        mda.core.flags['convert_gromacs_lengths'] = False
        self.universe = mda.Universe(PDB, self.filename)
        self.ts = self.universe.trajectory.ts
    def tearDown(self):
        mda.core.flags['convert_gromacs_lengths'] = True  # default
        del self.universe
        del self.ts

    @dec.slow
    def test_coordinates(self):
        # note: these are the native coordinates in nm; for the test to succeed:
        assert_equal(mda.core.flags['convert_gromacs_lengths'], False, 
                     "oops, mda.core.flags['convert_gromacs_lengths'] should be False for this test")
        ca_nm = np.array([[ 6.043369675,  7.385184479,  1.381425762]], dtype=np.float32)
        U = self.universe
        T = U.trajectory
        T.rewind()
        T.next()
        T.next()
        assert_equal(self.ts.frame, 3, "failed to step to frame 3")        
        ca = U.selectAtoms('name CA and resid 122')
        # low precision match because we also look at the trr: only 3 decimals in nm in xtc!
        assert_array_almost_equal(ca.coordinates(), ca_nm, 3,
                                  err_msg="native coords of Ca of resid 122 do not match for frame 3 "\
                                      "with core.flags['convert_gromacs_lengths'] = False")

class TestXTCNoConversion(_XDRNoConversion):
    filename = XTC

class TestTRRNoConversion(_XDRNoConversion):
    filename = TRR

class _GromacsWriter(TestCase):
    infilename = None  # XTC or TRR
    Writers = {'.trr': MDAnalysis.coordinates.TRR.TRRWriter,
               '.xtc': MDAnalysis.coordinates.XTC.XTCWriter,
               }

    def setUp(self):
        self.universe = mda.Universe(GRO, self.infilename)
        ext = os.path.splitext(self.infilename)[1]
        fd, self.outfile = tempfile.mkstemp(suffix=ext)
        self.Writer = self.Writers[ext]

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass
        del self.universe
        del self.Writer

    @dec.slow
    @attr('issue')    
    def test_write_trajectory(self):
        """Test writing Gromacs trajectories (Issue 38)"""
        t = self.universe.trajectory
        W = self.Writer(self.outfile, t.numatoms, delta=t.delta, step=t.skip_timestep)
        for ts in self.universe.trajectory:
            W.write_next_timestep(ts)
        W.close()

        uw = mda.Universe(GRO, self.outfile)

        # check that the coordinates are identical for each time step
        for orig_ts, written_ts in itertools.izip(self.universe.trajectory, uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, 3,
                                      err_msg="coordinate mismatch between original and written trajectory at frame %d (orig) vs %d (written)" % (orig_ts.frame, written_ts.frame))
        
class TestXTCWriter(_GromacsWriter):
    infilename = XTC

class TestTRRWriter(_GromacsWriter):
    infilename = TRR
