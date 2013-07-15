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
from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small, GRO, TRR
import MDAnalysis.core.AtomGroup
from MDAnalysis.core.AtomGroup import Atom, AtomGroup
from MDAnalysis import NoDataError

import numpy
from numpy.testing import *
from numpy import array, float32, rad2deg
from nose.plugins.attrib import attr

import os
import tempfile

try:
    from numpy.testing import assert_
except ImportError:
    # missing in numpy 1.2 but needed here:
    # copied code from numpy.testing 1.5
    def assert_(val, msg='') :
        """
        Assert that works in release mode.

        The Python built-in ``assert`` does not work when executing code in
        optimized mode (the ``-O`` flag) - no byte-code is generated for it.

        For documentation on usage, refer to the Python documentation.

        """
        if not val :
            raise AssertionError(msg)

class TestAtom(TestCase):
    """Tests of Atom."""
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.atom = self.universe.atoms[1000]  # Leu67:CG
        self.known_pos = array([  3.94543672, -12.4060812 ,  -7.26820087], dtype=float32)

    def tearDown(self):
        del self.universe
        del self.atom

    def test_attributes_names(self):
        a = self.atom
        assert_equal(a.name, 'CG')
        assert_equal(a.resname, 'LEU')

    def test_attributes_pos(self):
        # old pos property
        assert_almost_equal(self.atom.pos, self.known_pos)
        def set_pos():
            self.atom.pos = self.known_pos + 3.14
        assert_raises(AttributeError, set_pos)

    def test_attributes_positions(self):
        a = self.atom
        # new position property (mutable)
        assert_almost_equal(a.position, self.known_pos)
        pos = a.position + 3.14
        a.position = pos
        assert_almost_equal(a.position, pos)

    def test_atom_selection(self):
        asel =  self.universe.selectAtoms('atom 4AKE 67 CG').atoms[0]
        assert_equal(self.atom, asel)

    def test_hierarchy(self):
        u = self.universe
        a = self.atom
        assert_equal(a.segment, u.s4AKE)
        assert_equal(a.residue, u.residues[66])

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

    def test_residues(self):
        u = self.universe
        assert_equal(u.residues[100]._atoms,
                     u.selectAtoms('resname ILE and resid 101')._atoms,
                     "Direct selection from residue group does not match expected I101.")

    def test_segments(self):
        u = self.universe
        assert_equal(u.segments.s4AKE._atoms,
                     u.selectAtoms('segid 4AKE')._atoms,
                     "Direct selection of segment 4AKE from segments failed.")

    def test_index_integer(self):
        u = self.universe
        a = u.atoms[100]
        assert_(isinstance(a, Atom), "integer index did not return Atom")

    def test_index_slice(self):
        u = self.universe
        a = u.atoms[100:200:10]
        assert_(isinstance(a, AtomGroup), "slice index did not return AtomGroup")

    def test_index_slice_empty(self):
        u = self.universe
        assert_array_equal(u.atoms[0:0], [], "making an empty AtomGroup failed")

    def test_index_advancedslice(self):
        u = self.universe
        aslice = [0, 10, 20, -1, 10]
        ag = u.atoms[aslice]
        assert_(isinstance(ag, AtomGroup),
                "advanced slicing does not produce a AtomGroup")
        assert_equal(ag[1], ag[-1], "advanced slicing does not preserve order")

    def test_phi_selection(self):
        phisel = self.universe.s4AKE.r10.phi_selection()
        assert_equal(phisel.names(), ['C', 'N', 'CA', 'C'])
        assert_equal(phisel.resids(), [9, 10])
        assert_equal(phisel.resnames(), ['PRO', 'GLY'])

    def test_psi_selection(self):
        psisel = self.universe.s4AKE.r10.psi_selection()
        assert_equal(psisel.names(), ['N', 'CA', 'C', 'N'])
        assert_equal(psisel.resids(), [10, 11])
        assert_equal(psisel.resnames(), ['GLY', 'ALA'])

    def test_omega_selection(self):
        osel =  self.universe.s4AKE.r8.omega_selection()
        assert_equal(osel.names(), ['CA', 'C', 'N', 'CA'])
        assert_equal(osel.resids(), [8, 9])
        assert_equal(osel.resnames(), ['ALA', 'PRO'])

    def test_chi1_selection(self):
        sel =  self.universe.s4AKE.r8.chi1_selection()
        assert_equal(sel, None)  # ALA
        sel =  self.universe.s4AKE.r13.chi1_selection()  # LYS
        assert_equal(sel.names(), ['N', 'CA', 'CB', 'CG'])
        assert_equal(sel.resids(), [13])
        assert_equal(sel.resnames(), ['LYS'])

    def test_dihedral_phi(self):
        u = self.universe
        u.trajectory.rewind()   # just to make sure...
        phisel = u.s4AKE.r10.phi_selection()
        assert_almost_equal(phisel.dihedral(), -168.57384, 3)

    def test_dihedral_psi(self):
        u = self.universe
        u.trajectory.rewind()   # just to make sure...
        psisel = u.s4AKE.r10.psi_selection()
        assert_almost_equal(psisel.dihedral(), -30.064838, 3)

    def test_dihedral_omega(self):
        u = self.universe
        u.trajectory.rewind()   # just to make sure...
        osel =  u.s4AKE.r8.omega_selection()
        assert_almost_equal(osel.dihedral(), -179.93439, 3)

    def test_dihedral_chi1(self):
        u = self.universe
        u.trajectory.rewind()   # just to make sure...
        sel =  u.s4AKE.r13.chi1_selection()  # LYS
        assert_almost_equal(sel.dihedral(), -58.428127, 3)

    def test_dihedral_ValueError(self):
        """test that AtomGroup.dihedral() raises ValueError if not exactly 4 atoms given"""
        nodih = self.universe.selectAtoms("resid 3:10")
        assert_raises(ValueError, nodih.dihedral)
        nodih = self.universe.selectAtoms("resid 3:5")
        assert_raises(ValueError, nodih.dihedral)

    def test_improper(self):
        u = self.universe
        u.trajectory.rewind()   # just to make sure...
        peptbond =  u.selectAtoms("atom 4AKE 20 C", "atom 4AKE 21 CA",
                                  "atom 4AKE 21 N", "atom 4AKE 21 HN")
        assert_almost_equal(peptbond.improper(), 168.52952575683594, 3,
                            "Peptide bond improper dihedral for M21 calculated wrongly.")

    def test_dihedral_equals_improper(self):
        u = self.universe
        u.trajectory.rewind()   # just to make sure...
        peptbond =  u.selectAtoms("atom 4AKE 20 C", "atom 4AKE 21 CA",
                                  "atom 4AKE 21 N", "atom 4AKE 21 HN")
        assert_equal(peptbond.improper(), peptbond.dihedral(),
                     "improper() and proper dihedral() give different results")


    def test_bond(self):
        self.universe.trajectory.rewind()   # just to make sure...
        sel2 = self.universe.s4AKE.r98.selectAtoms("name OE1", "name OE2")
        assert_almost_equal(sel2.bond(), 2.1210737228393555, 3,
                            "distance of Glu98 OE1--OE2 wrong")

    def test_angle(self):
        self.universe.trajectory.rewind()   # just to make sure...
        sel3 = self.universe.s4AKE.r98.selectAtoms("name OE1", "name CD", "name OE2")
        assert_almost_equal(sel3.angle(), 117.46187591552734, 3,
                            "angle of Glu98 OE1-CD-OE2 wrong")

    def test_shapeParameter(self):
        s = self.universe.s4AKE.shapeParameter()
        assert_almost_equal(s, 0.00240753939086033, 6)

    def test_asphericity(self):
        a = self.universe.s4AKE.asphericity()
        assert_almost_equal(a, 0.020227504542775828, 6)

    # TODO: tests for the coordinate manipulation methods
    # - transform
    # - translate
    # - rotate
    # - rotateby

    def test_positions(self):
        ag = self.universe.selectAtoms("bynum 12:42")
        pos = ag.positions + 3.14
        ag.positions = pos
        # shoud work
        assert_almost_equal(ag.coordinates(), pos,
                            err_msg="failed to update atoms 12:42 position to new position")
        def set_badarr(pos=pos):
            # create wrong size array
            badarr = numpy.random.random((pos.shape[0] - 1, pos.shape[1] - 1))
            ag.positions = badarr
        assert_raises(ValueError, set_badarr)

    def test_set_positions(self):
        ag = self.universe.selectAtoms("bynum 12:42")
        pos = ag.get_positions() + 3.14
        ag.set_positions(pos)
        assert_almost_equal(ag.coordinates(), pos,
                            err_msg="failed to update atoms 12:42 position to new position")

    def test_no_velocities_raises_NoDataError(self):
        def get_vel(ag=self.universe.selectAtoms("bynum 12:42")):
            v = ag.get_velocities()
        # trj has no velocities
        assert_raises(NoDataError, get_vel)

class TestAtomGroupVelocities(TestCase):
    """Tests of velocity-related functions in AtomGroup"""
    def setUp(self):
        self.universe = MDAnalysis.Universe(GRO, TRR)
        self.ag = self.universe.selectAtoms("bynum 12:42")

    @dec.slow
    def test_get_velocities(self):
        v = self.ag.get_velocities()
        assert_(numpy.any(numpy.abs(v) > 1e-6), "velocities should be non-zero")

    @dec.slow
    def test_velocities(self):
        ag = self.universe.atoms[42:45]
        ref_v = numpy.array([[ -3.61757946,  -4.9867239 ,   2.46281552],
                             [  2.57792854,   3.25411797,  -0.75065529],
                             [ 13.91627216,  30.17778587, -12.16669178]])
        v = ag.velocities()
        # TODO: for 0.8 change to
        ##v = self.ag.velocities
        assert_almost_equal(v, ref_v, err_msg="velocities were not read correctly")

    @dec.slow
    def test_set_velocities(self):
        ag = self.ag
        v = ag.get_velocities() - 2.7271
        ag.set_velocities(v)
        assert_almost_equal(ag.get_velocities(), v,
                            err_msg="messages were not set to new value")

def test_empty_AtomGroup():
    """Test that a empty AtomGroup can be constructed (Issue 12)"""
    ag = MDAnalysis.core.AtomGroup.AtomGroup([])
    assert_equal(len(ag), 0)

class _WriteAtoms(TestCase):
    """Set up the standard AdK system in implicit solvent."""
    ext = None   # override to test various output writers
    precision = 3

    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        suffix = '.' + self.ext
        fd, self.outfile = tempfile.mkstemp(suffix=suffix)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.universe

    def universe_from_tmp(self):
        return MDAnalysis.Universe(self.outfile, convert_units=True)

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
        assert_equal(len(u2.atoms), len(CA.atoms), "written CA selection does not match original selection")
        assert_almost_equal(CA2.coordinates(), CA.coordinates(), self.precision,
                            err_msg="CA coordinates do not agree with original")

    def test_write_Residue(self):
        G = self.universe.s4AKE.ARG[-2]   # 2nd but last Arg
        G.write(self.outfile)
        u2 = self.universe_from_tmp()
        G2 = u2.selectAtoms('all')   # check EVERYTHING, otherwise we might get false positives!
        assert_equal(len(u2.atoms), len(G.atoms), "written R206 Residue does not match original ResidueGroup")
        assert_almost_equal(G2.coordinates(), G.coordinates(), self.precision,
                            err_msg="Residue R206 coordinates do not agree with original")

    def test_write_ResidueGroup(self):
        G = self.universe.s4AKE.LEU
        G.write(self.outfile)
        u2 = self.universe_from_tmp()
        G2 = u2.selectAtoms('all')   # check EVERYTHING, otherwise we might get false positives!
        assert_equal(len(u2.atoms), len(G.atoms), "written LEU ResidueGroup does not match original ResidueGroup")
        assert_almost_equal(G2.coordinates(), G.coordinates(), self.precision,
                            err_msg="ResidueGroup LEU coordinates do not agree with original")

    def test_write_Segment(self):
        G = self.universe.s4AKE
        G.write(self.outfile)
        u2 = self.universe_from_tmp()
        G2 = u2.selectAtoms('all')   # check EVERYTHING, otherwise we might get false positives!
        assert_equal(len(u2.atoms), len(G.atoms), "written s4AKE segment does not match original segment")
        assert_almost_equal(G2.coordinates(), G.coordinates(), self.precision,
                            err_msg="segment s4AKE coordinates do not agree with original")

    def test_write_Universe(self):
        U = self.universe
        W = MDAnalysis.Writer(self.outfile)
        W.write(U)
        W.close()
        u2 = self.universe_from_tmp()
        assert_equal(len(u2.atoms), len(U.atoms), "written 4AKE universe does not match original universe in size")
        assert_almost_equal(u2.atoms.coordinates(), U.atoms.coordinates(), self.precision,
                            err_msg="written universe 4AKE coordinates do not agree with original")


class TestWritePDB(_WriteAtoms):
    ext = "pdb"
    precision = 3

import MDAnalysis.coordinates
class TestWriteCRD(_WriteAtoms):
    ext = "crd"
    precision = 5

class TestWriteGRO(_WriteAtoms):
    ext = "gro"
    precision = 2

    def test_flag_convert_gromacs_length(self):
        assert_equal(MDAnalysis.core.flags['convert_gromacs_lengths'], True,
                     "The flag convert_gromacs_lengths SHOULD be True by default! "
                     "(If it is not then this might indicate a race condition in the "
                     "testing suite.)")


import MDAnalysis.core.AtomGroup
@attr("issue")
def test_generated_residueselection():
    """Test that a generated residue group always returns a ResidueGroup (Issue 47)"""
    universe = MDAnalysis.Universe(PSF, DCD)
    # only a single Cys in AdK
    cys = universe.s4AKE.CYS
    assert_(isinstance(cys, MDAnalysis.core.AtomGroup.ResidueGroup),
            "Single Cys77 is NOT returned as a ResidueGroup with a single Residue (Issue 47)")

    # multiple Met
    met = universe.s4AKE.MET
    assert_(isinstance(met, MDAnalysis.core.AtomGroup.ResidueGroup),
            "Met selection does not return a ResidueGroup")

    del universe

@attr('issue')
def test_instantselection_termini():
    """Test that instant selections work, even for residues that are also termini (Issue 70)"""
    universe = MDAnalysis.Universe(PSF, DCD)
    assert_equal(universe.residues[20].CA.name, 'CA', "CA of MET21 is not selected correctly")
    del universe


class TestUniverse(TestCase):
    def test_load(self):
        u = MDAnalysis.Universe(PSF, PDB_small)
        assert_equal(len(u.atoms), 3341, "Loading universe failed somehow")

    @attr('issue')
    def test_load_new(self):
        u = MDAnalysis.Universe(PSF, DCD)
        u.load_new(PDB_small)
        assert_equal(len(u.trajectory), 1, "Failed to load_new(PDB)")

    def test_load_new_strict(self):
        u = MDAnalysis.Universe(PSF, DCD)
        u.load_new(PDB_small, permissive=False)
        assert_equal(len(u.trajectory), 1, "Failed to load_new(PDB, permissive=False)")

    def test_load_new_permissive(self):
        u = MDAnalysis.Universe(PSF, DCD)
        u.load_new(PDB_small, permissive=True)
        assert_equal(len(u.trajectory), 1, "Failed to load_new(PDB, permissive=True)")
