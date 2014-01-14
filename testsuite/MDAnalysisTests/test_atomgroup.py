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
from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small, GRO, TRR, \
                                      merge_protein, merge_water, merge_ligand
import MDAnalysis.core.AtomGroup
from MDAnalysis.core.AtomGroup import Atom, AtomGroup
from MDAnalysis import NoDataError

import numpy
from numpy.testing import *
from numpy import array, float32, rad2deg
from nose.plugins.attrib import attr

import os
import tempfile
import itertools

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
        self.dih_prec = 2

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

    def test_coordinates(self):
        assert_array_almost_equal(self.ag.coordinates()[1000:2000:200],
                                  array([[  3.94543672, -12.4060812 ,  -7.26820087],
                                         [ 13.21632767,   5.879035  , -14.67914867],
                                         [ 12.07735443,  -9.00604534,   4.09301519],
                                         [ 11.35541916,   7.0690732 ,  -0.32511973],
                                         [-13.26763439,   4.90658951,  10.6880455 ]], dtype=float32))

    def test_principalAxes(self):
        assert_array_almost_equal(self.ag.principalAxes(),
                                  array([[ -9.99925632e-01,   1.21546132e-02,   9.98264877e-04],
                                         [  1.20986911e-02,   9.98951474e-01,  -4.41539838e-02],
                                         [  1.53389276e-03,   4.41386224e-02,   9.99024239e-01]]))

    def test_totalCharge(self):
        assert_almost_equal(self.ag.totalCharge(), -4.0)

    def test_totalMass(self):
        assert_almost_equal(self.ag.totalMass(), 23582.043)

    def test_indices_ndarray(self):
        assert_equal(isinstance(self.ag.indices(), numpy.ndarray), True)
    def test_indices(self):
        assert_array_equal(self.ag.indices()[:5], array([0, 1, 2, 3, 4]))

    def test_resids_ndarray(self):
        assert_equal(isinstance(self.ag.resids(), numpy.ndarray), True)
    def test_resids(self):
        assert_array_equal(self.ag.resids(), numpy.arange(1,215))

    def test_resnums_ndarray(self):
        assert_equal(isinstance(self.ag.resnums(), numpy.ndarray), True)
    def test_resnums(self):
        assert_array_equal(self.ag.resids(), numpy.arange(1,215))

    def test_resnames_ndarray(self):
        assert_equal(isinstance(self.ag.resnames(), numpy.ndarray), True)
    def test_resnames(self):
        resnames = self.ag.resnames()
        assert_array_equal(resnames[0:3], numpy.array(["MET", "ARG", "ILE"]))

    def test_names_ndarray(self):
        assert_equal(isinstance(self.ag.names(), numpy.ndarray), True)
    def test_names(self):
        names = self.ag.names()
        assert_array_equal(names[0:3], numpy.array(["N", "HT1", "HT2"]))

    def test_segids_ndarray(self):
        assert_equal(isinstance(self.ag.segids(), numpy.ndarray), True)
    def test_segids(self):
        segids = self.ag.segids()
        assert_array_equal(segids[0], numpy.array(["4AKE"]))

    def test_masses_ndarray(self):
        assert_equal(isinstance(self.ag.masses(), numpy.ndarray), True)
    def test_masses(self):
        masses = self.ag.masses()
        assert_array_equal(masses[0:3], numpy.array([14.007, 1.008, 1.008]))

    def test_charges_ndarray(self):
        assert_equal(isinstance(self.ag.charges(), numpy.ndarray), True)
    def test_charges(self):
        assert_array_almost_equal(self.ag.charges()[1000:2000:200],
                                  array([-0.09,  0.09, -0.47,  0.51,  0.09]))

    def test_radii_ndarray(self):
        assert_equal(isinstance(self.ag.radii(), numpy.ndarray), True)
    def test_radii(self):
        radii = self.ag.radii()
        assert_array_equal(radii[0:3], numpy.array([None, None, None]))

    def test_bfactors_ndarray(self):
        assert_equal(isinstance(self.ag.bfactors, numpy.ndarray), True)
    def test_bfactors(self):
        bfactors = self.ag.bfactors   # property, not method!
        assert_array_equal(bfactors[0:3], numpy.array([None, None, None]))

    # TODO: add all other methods except selectAtoms(), see test_selections.py

    # add new methods here...
    def test_packintobox(self):
        """test AtomGroup.packintobox(): Tests application of periodic boundary conditions on coordinates

        Reference system doesn't have dimensions, so an arbitrary box is imposed on the system
        """
        u = self.universe
        u.trajectory.rewind() #just to make sure...
        ag = u.atoms

        def badpack(ag=ag):
            ag.packintobox()
        assert_raises(ValueError, badpack) #This system has no dimensions, so by default can't pack

        ag.packintobox(box=numpy.array([5.,5.,5.])) #Provide arbitrary box
        assert_array_almost_equal(ag.coordinates()[1000:2000:200],
                                  array([[ 3.94543672,  2.5939188 ,  2.73179913],
                                         [ 3.21632767,  0.879035  ,  0.32085133],
                                         [ 2.07735443,  0.99395466,  4.09301519],
                                         [ 1.35541916,  2.0690732 ,  4.67488003],
                                         [ 1.73236561,  4.90658951,  0.6880455 ]], dtype=float32))

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
        assert_almost_equal(phisel.dihedral(), -168.57384, self.dih_prec)

    def test_dihedral_psi(self):
        u = self.universe
        u.trajectory.rewind()   # just to make sure...
        psisel = u.s4AKE.r10.psi_selection()
        assert_almost_equal(psisel.dihedral(), -30.064838, self.dih_prec)

    def test_dihedral_omega(self):
        u = self.universe
        u.trajectory.rewind()   # just to make sure...
        osel =  u.s4AKE.r8.omega_selection()
        assert_almost_equal(osel.dihedral(), -179.93439, self.dih_prec)

    def test_dihedral_chi1(self):
        u = self.universe
        u.trajectory.rewind()   # just to make sure...
        sel =  u.s4AKE.r13.chi1_selection()  # LYS
        assert_almost_equal(sel.dihedral(), -58.428127, self.dih_prec)

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
        assert_almost_equal(peptbond.improper(), 168.52952575683594, self.dih_prec,
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
        # should work
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

    def test_set_resid(self):
        ag = self.universe.selectAtoms("bynum 12:42")
        resid = 999
        ag.set_resid(resid)
        # check individual atoms
        assert_equal([a.resid for a in ag],
                     resid*numpy.ones(ag.numberOfAtoms()),
                     err_msg="failed to set_resid atoms 12:42 to same resid")
        # check residues
        assert_equal(ag.resids(), 999*numpy.ones(ag.numberOfResidues()),
                     err_msg="failed to set_resid of residues belonging to atoms 12:42 to same resid")

    def test_set_resids(self):
        """test_set_resid: set AtomGroup resids on a per-atom basis"""
        ag = self.universe.selectAtoms("bynum 12:42")
        resids = numpy.array([a.resid for a in ag]) + 1000
        ag.set_resid(resids)
        # check individual atoms
        assert_equal([a.resid for a in ag], resids,
                     err_msg="failed to set_resid atoms 12:42 to resids {0}".format(resids))
        # check residues
        assert_equal(ag.resids(), numpy.unique(resids),
                     err_msg="failed to set_resid of residues belonging to atoms 12:42 to same resid")

    def test_merge_residues(self):
        ag = self.universe.selectAtoms("resid 12:14")
        nres_old = self.universe.atoms.numberOfResidues()
        natoms_old = ag.numberOfAtoms()
        ag.set_resid(12)    # merge all into one with resid 12
        nres_new = self.universe.atoms.numberOfResidues()
        r_merged = self.universe.selectAtoms("resid 12:14").residues
        natoms_new = self.universe.selectAtoms("resid 12").numberOfAtoms()
        assert_equal(len(r_merged), 1, err_msg="set_resid failed to merge residues: merged = {0}".format(r_merged))
        assert_equal(nres_new, nres_old - 2, err_msg="set_resid failed to merge residues: merged = {0}".format(r_merged))
        assert_equal(natoms_new, natoms_old, err_msg="set_resid lost atoms on merge".format(r_merged))

    def test_set_mass(self):
        ag = self.universe.selectAtoms("bynum 12:42 and name H*")
        mass = 2.0
        ag.set_mass(mass)
        # check individual atoms
        assert_equal([a.mass for a in ag],
                     mass*numpy.ones(ag.numberOfAtoms()),
                     err_msg="failed to set_mass H* atoms in resid 12:42 to {0}".format(mass))

    def test_set_segid(self):
        u = self.universe
        u.selectAtoms("(resid 1-29 or resid 60-121 or resid 160-214)").set_segid("CORE")
        u.selectAtoms("resid 122-159").set_segid("LID")
        u.selectAtoms("resid 30-59").set_segid("NMP")
        assert_equal(u.atoms.segids(), ["CORE", "NMP", "CORE", "LID", "CORE"],
                     err_msg="failed to change segids = {0}".format(u.atoms.segids()))

class TestResidue(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.res = self.universe.residues[100]

    def test_type(self):
        assert_equal(type(self.res), MDAnalysis.core.AtomGroup.Residue)
        assert_equal(self.res.name, "ILE")
        assert_equal(self.res.id, 101)

    def test_index(self):
        atom = self.res[2]
        assert_equal(type(atom), MDAnalysis.core.AtomGroup.Atom)
        assert_equal(atom.name, "CA")
        assert_equal(atom.number, 1522)
        assert_equal(atom.resid, 101)

    def test_slicing(self):
        atoms = self.res[2:10:2]
        assert_equal(len(atoms), 4)
        assert_equal(type(atoms), MDAnalysis.core.AtomGroup.AtomGroup)

    def test_advanced_slicing(self):
        atoms = self.res[[0, 2, -2, -1]]
        assert_equal(len(atoms), 4)
        assert_equal(type(atoms), MDAnalysis.core.AtomGroup.AtomGroup)
        assert_equal(atoms.names(), ["N", "CA", "C" ,"O"])


class TestResidueGroup(TestCase):
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.rg = self.universe.atoms.residues

    def test_newResidueGroup(self):
        """test that slicing a ResidueGroup returns a new ResidueGroup (Issue 135)"""
        rg = self.universe.atoms.residues
        newrg = rg[10:20:2]
        assert_equal(type(newrg), type(rg), "Failed to make a new ResidueGroup: type mismatch")
        assert_equal(len(newrg), len(rg[10:20:2]))

    def test_numberOfAtoms(self):
        assert_equal(self.rg.numberOfAtoms(), 3341)

    def test_numberOfResidues(self):
        assert_equal(self.rg.numberOfResidues(), 214)

    def test_len(self):
        """testing that len(residuegroup) == residuegroup.numberOfResidues()"""
        assert_equal(len(self.rg), self.rg.numberOfResidues(), "len and numberOfResidues() disagree")

    def test_set_resid(self):
        rg = self.universe.selectAtoms("bynum 12:42").residues
        resid = 999
        rg.set_resid(resid)
        # check individual atoms
        assert_equal([a.resid for a in rg.atoms],
                     resid*numpy.ones(rg.numberOfAtoms()),
                     err_msg="failed to set_resid atoms 12:42 to same resid")
        # check residues
        assert_equal(rg.resids(), resid*numpy.ones(rg.numberOfResidues()),
                     err_msg="failed to set_resid of residues belonging to atoms 12:42 to same resid")

    def test_set_resids(self):
        """test_set_resid: set ResidueGroup resids on a per-residue basis"""
        rg = self.universe.selectAtoms("resid 10:18").residues
        resids = numpy.array(rg.resids()) + 1000
        rg.set_resid(resids)
        # check individual atoms
        for r, resid in itertools.izip(rg, resids):
            assert_equal([a.resid for a in r.atoms],
                         resid*numpy.ones(r.numberOfAtoms()),
                         err_msg="failed to set_resid residues 10:18 to same resid in residue {0}\n"
                         "(resids = {1}\nresidues = {2})".format(r, resids, rg))
        # check residues
        # NOTE: need to create a new selection because underlying Residue objects are not changed;
        #       only Atoms are changed, and Residues are rebuilt from Atoms.
        rgnew = self.universe.selectAtoms("resid 1010:1018").residues
        assert_equal(rgnew.resids(), numpy.unique(resids),
                     err_msg="failed to set_resid of residues belonging to residues 10:18 to new resids")

    def test_set_resids_updates_self(self):
        rg = self.universe.selectAtoms("resid 10:18").residues
        resids = numpy.array(rg.resids()) + 1000
        rg.set_resid(resids)
        #rgnew = self.universe.selectAtoms("resid 1000:1008").residues
        assert_equal(rg.resids(), numpy.unique(resids),
                     err_msg="old selection was not changed in place after set_resid")

    def test_merge_residues(self):
        rg = self.universe.selectAtoms("resid 12:14").residues
        nres_old = self.universe.atoms.numberOfResidues()
        natoms_old = rg.numberOfAtoms()
        rg.set_resid(12)    # merge all into one with resid 12
        nres_new = self.universe.atoms.numberOfResidues()
        r_merged = self.universe.selectAtoms("resid 12:14").residues
        natoms_new = self.universe.selectAtoms("resid 12").numberOfAtoms()
        assert_equal(len(r_merged), 1, err_msg="set_resid failed to merge residues: merged = {0}".format(r_merged))
        assert_equal(nres_new, nres_old - 2, err_msg="set_resid failed to merge residues: merged = {0}".format(r_merged))
        assert_equal(natoms_new, natoms_old, err_msg="set_resid lost atoms on merge".format(r_merged))

        assert_equal(self.universe.residues.numberOfResidues(), self.universe.atoms.numberOfResidues(),
                     err_msg="Universe.residues and Universe.atoms.numberOfResidues() do not agree after residue merge.")

    def test_set_mass(self):
        rg = self.universe.selectAtoms("bynum 12:42 and name H*").residues
        mass = 2.0
        rg.set_mass(mass)
        # check individual atoms
        assert_equal([a.mass for a in rg.atoms],
                     mass*numpy.ones(rg.numberOfAtoms()),
                     err_msg="failed to set_mass H* atoms in resid 12:42 to {0}".format(mass))

class TestSegment(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.universe.residues[:100].set_segid("A")     # make up some segments
        self.universe.residues[100:150].set_segid("B")
        self.universe.residues[150:].set_segid("C")
        self.sB = self.universe.segments[1]

    def test_type(self):
        assert_equal(type(self.sB), MDAnalysis.core.AtomGroup.Segment)
        assert_equal(self.sB.name, "B")

    def test_index(self):
        s = self.sB
        res = s[5]
        assert_equal(type(res), MDAnalysis.core.AtomGroup.Residue)

    def test_slicing(self):
        res = self.sB[5:10]
        assert_equal(len(res), 5)
        assert_equal(type(res), MDAnalysis.core.AtomGroup.ResidueGroup)

    def test_advanced_slicing(self):
        res = self.sB[[3,7,2,4]]
        assert_equal(len(res), 4)
        assert_equal(type(res), MDAnalysis.core.AtomGroup.ResidueGroup)


class TestSegmentGroup(TestCase):
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.g = self.universe.atoms.segments

    def test_newSegmentGroup(self):
        """test that slicing a SegmentGroup returns a new SegmentGroup (Issue 135)"""
        g = self.universe.atoms.segments
        newg = g[:]
        assert_equal(type(newg), type(g), "Failed to make a new SegmentGroup: type mismatch")
        assert_equal(len(newg), len(g))

    def test_numberOfAtoms(self):
        assert_equal(self.g.numberOfAtoms(), 3341)

    def test_numberOfResidues(self):
        assert_equal(self.g.numberOfResidues(), 214)

    def test_set_resid(self):
        g = self.universe.selectAtoms("bynum 12:42").segments
        resid = 999
        g.set_resid(resid)
        # check individual atoms
        assert_equal([a.resid for a in g.atoms],
                     resid*numpy.ones(g.numberOfAtoms()),
                     err_msg="failed to set_resid for segment to same resid")
        # check residues
        assert_equal(g.resids(), resid*numpy.ones(g.numberOfResidues()),
                     err_msg="failed to set_resid of segments belonging to atoms 12:42 to same resid")

    def test_set_resids(self):
        g = self.universe.selectAtoms("resid 10:18").segments
        resid = 999
        g.set_resid(resid*numpy.ones(len(g)))
        # note: all is now one residue... not meaningful but it is the correct behaviour
        assert_equal(g.atoms.resids(), [resid],
                     err_msg="failed to set_resid  in Segment {0}".format(g))

    def test_set_segid(self):
        s = self.universe.selectAtoms('all').segments
        s.set_segid(['ADK'])
        assert_equal(self.universe.segments.segids(), ['ADK'],
                     err_msg="failed to set_segid on segments")

    def test_set_segid_updates_self(self):
        g = self.universe.selectAtoms("resid 10:18").segments
        g.set_segid('ADK')
        assert_equal(g.segids(), ['ADK'],
                     err_msg="old selection was not changed in place after set_segid")

    def test_set_mass(self):
        g = self.universe.selectAtoms("bynum 12:42 and name H*").segments
        mass = 2.0
        g.set_mass(mass)
        # check individual atoms
        assert_equal([a.mass for a in g.atoms],
                     mass*numpy.ones(g.numberOfAtoms()),
                     err_msg="failed to set_mass in segment of  H* atoms in resid 12:42 to {0}".format(mass))

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
        v = ag.velocities
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
        os.close(fd)

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


class TestMerge(TestCase):
    ext = "pdb"
    def setUp(self):
        u1 = MDAnalysis.Universe(merge_protein)
        u2 = MDAnalysis.Universe(merge_ligand)
        u3 = MDAnalysis.Universe(merge_water)
        self.universes = [u1,u2,u3]
        
        suffix = '.' + self.ext
        fd, self.outfile = tempfile.mkstemp(suffix=suffix)
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        for u in self.universes: del u

    
    
    
    def test_merge(self):
        u1,u2,u3 = self.universes
        
        ids_before = [a.number for u in [u1, u2, u3] for a in u.atoms]  
        
        # Do the merge
        u0 = MDAnalysis.Merge(u1.atoms, u2.atoms, u3.atoms)
        
        # Check that the output Universe has the same number of atoms as the 
        # starting AtomGroups
        assert_equal(len(u0.atoms),(len(u1.atoms)+len(u2.atoms)+len(u3.atoms)))
        
        # Make sure that all the atoms in the new universe are assigned to only 
        # one, new Universe
        set0 = set([a.universe for a in u0.atoms])
        assert_equal(len(set0), 1)
        u = list(set0)[0]
        assert_equal(u, u0)
        
        # Make sure that the atom ids of the original universes are unchanged, 
        # ie we didn't make the original Universes 'dirty'
        ids_after = [a.number for u in [u1, u2, u3] for a in u.atoms]
        assert_equal(len(ids_after), (len(u1.atoms)+len(u2.atoms)+len(u3.atoms)))
        assert_equal(ids_before, ids_after)
        
        # Test that we have a same number of atoms in a different way
        ids_new = [a.number for a in u0.atoms]
        assert_equal(len(ids_new), len(ids_before))
    

        u0.atoms.write(self.outfile)
        u = MDAnalysis.Universe(self.outfile)
        ids_new2 = [a.number for a in u.atoms]
        assert_equal(ids_new, ids_new2)

class TestUniverse(TestCase):
    def test_load(self):
        # Universe(top, trj)
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

    def test_load_structure(self):
        # Universe(struct)
        ref = MDAnalysis.Universe(PSF, PDB_small)
        u = MDAnalysis.Universe(PDB_small)
        assert_equal(len(u.atoms), 3341, "Loading universe failed somehow")
        assert_almost_equal(u.atoms.positions, ref.atoms.positions)

    def test_load_multiple_list(self):
        # Universe(top, [trj, trj, ...])
        ref = MDAnalysis.Universe(PSF, DCD)
        u = MDAnalysis.Universe(PSF, [DCD, DCD])
        assert_equal(len(u.atoms), 3341, "Loading universe failed somehow")
        assert_equal(u.trajectory.numframes, 2*ref.trajectory.numframes)

    def test_load_multiple_args(self):
        # Universe(top, trj, trj, ...)
        ref = MDAnalysis.Universe(PSF, DCD)
        u = MDAnalysis.Universe(PSF, DCD, DCD)
        assert_equal(len(u.atoms), 3341, "Loading universe failed somehow")
        assert_equal(u.trajectory.numframes, 2*ref.trajectory.numframes)

