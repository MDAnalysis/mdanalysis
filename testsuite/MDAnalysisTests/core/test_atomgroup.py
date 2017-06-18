# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import division, absolute_import

from glob import glob
import itertools
import os
from os import path
import numpy as np
import warnings

from numpy.testing import (
    dec,
    assert_,
    assert_array_equal,
    assert_almost_equal,
    assert_raises,
    assert_equal,
    assert_array_almost_equal,
    raises
)
from nose.plugins.attrib import attr

import MDAnalysis as mda
from MDAnalysis.lib import transformations
from MDAnalysis.core.topologyobjects import (
    Bond,
    Angle,
    Dihedral,
    ImproperDihedral,
)

from MDAnalysisTests.datafiles import (
    PSF, DCD,
    TRZ_psf, TRZ,
    two_water_gro,
)
from MDAnalysisTests import parser_not_found, tempdir, make_Universe

# I want to catch all warnings in the tests. If this is not set at the start it
# could cause test that check for warnings to fail.
warnings.simplefilter('always')

class TestDeprecationWarnings(object):
    @staticmethod
    def test_AtomGroupUniverse_usage_warning():
        with warnings.catch_warnings(record=True) as warn:
            warnings.simplefilter('always')
            mda.core.AtomGroup.Universe(PSF, DCD)
        assert_equal(len(warn), 1)

    @staticmethod
    def test_old_AtomGroup_init_warns():
        u = make_Universe(('names',))
        at_list = list(u.atoms[:10])
        with warnings.catch_warnings(record=True) as warn:
            warnings.simplefilter('always')
            ag = mda.core.groups.AtomGroup(at_list)
        assert_equal(len(warn), 1)

    @staticmethod
    def test_old_AtomGroup_init_works():
        u = make_Universe(('names',))
        at_list = list(u.atoms[:10])
        ag = mda.core.groups.AtomGroup(at_list)

        assert_(isinstance(ag, mda.core.groups.AtomGroup))
        assert_(len(ag) == 10)
        assert_equal(ag.names, u.atoms[:10].names)

    @staticmethod
    def test_old_ResidueGroup_init_warns():
        u = make_Universe(('resnames',))
        res_list = list(u.residues[:10])
        with warnings.catch_warnings(record=True) as warn:
            warnings.simplefilter('always')
            rg = mda.core.groups.ResidueGroup(res_list)
        assert_equal(len(warn), 1)

    @staticmethod
    def test_old_ResidueGroup_init_works():
        u = make_Universe(('resnames',))
        res_list = list(u.residues[:10])
        rg = mda.core.groups.ResidueGroup(res_list)

        assert_(isinstance(rg, mda.core.groups.ResidueGroup))
        assert_(len(rg) == 10)
        assert_equal(rg.resnames, u.residues[:10].resnames)

    @staticmethod
    def test_old_SegmentGroup_init_warns():
        u = make_Universe(('segids',))
        seg_list = list(u.segments[:3])
        with warnings.catch_warnings(record=True) as warn:
            warnings.simplefilter('always')
            sg = mda.core.groups.SegmentGroup(seg_list)
        assert_equal(len(warn), 1)

    @staticmethod
    def test_old_SegmentGroup_init_works():
        u = make_Universe(('segids',))
        seg_list = list(u.segments[:3])
        sg = mda.core.groups.SegmentGroup(seg_list)

        assert_(isinstance(sg, mda.core.groups.SegmentGroup))
        assert_(len(sg) == 3)
        assert_equal(sg.segids, u.segments[:3].segids)


class TestAtomGroupToTopology(object):
    """Test the conversion of AtomGroup to TopologyObjects"""
    def setUp(self):
        self.u = mda.Universe(PSF, DCD)

    def tearDown(self):
        del self.u

    def test_bond(self):
        ag = self.u.atoms[:2]
        bond = ag.bond
        assert_(isinstance(bond, Bond))

    def test_angle(self):
        ag = self.u.atoms[:3]
        angle = ag.angle
        assert_(isinstance(angle, Angle))

    def test_dihedral(self):
        ag = self.u.atoms[:4]
        dih = ag.dihedral
        assert_(isinstance(dih, Dihedral))

    def test_improper(self):
        ag = self.u.atoms[:4]
        imp = ag.improper
        assert_(isinstance(imp, ImproperDihedral))

    def _check_VE(self, btype):
        ag = self.u.atoms[:10]

        assert_raises(ValueError, getattr, ag, btype)

    def test_VEs(self):
        for btype in ('bond', 'angle', 'dihedral', 'improper'):
            yield self._check_VE, btype


class TestAtomGroupWriting(object):
    def setUp(self):
        self.u = mda.Universe(PSF, DCD)

    def tearDown(self):
        del self.u

    def test_write_no_args(self):
        with tempdir.in_tempdir():
            self.u.atoms.write()
            files = glob('*')
            assert_equal(len(files), 1)

            name = path.splitext(path.basename(DCD))[0]
            assert_equal(files[0], "{}_0.pdb".format(name))

    @raises(ValueError)
    def test_raises(self):
        with tempdir.in_tempdir():
            self.u.atoms.write('useless.format123')

    def test_write_coordinates(self):
        with tempdir.in_tempdir():
            self.u.atoms.write("test.xtc")

    def test_write_selection(self):
        with tempdir.in_tempdir():
            self.u.atoms.write("test.vmd")

    def test_bogus_kwarg_pdb(self):
        # test for resolution of Issue 877
        with tempdir.in_tempdir():
            with assert_raises(TypeError):
                self.u.atoms.write('dummy.pdb', bogus="what?")

class _WriteAtoms(object):
    """Set up the standard AdK system in implicit solvent."""
    ext = None  # override to test various output writers
    precision = 3

    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        suffix = '.' + self.ext
        self.tempdir = tempdir.TempDir()
        self.outfile = os.path.join(self.tempdir.name, 'writeatoms' + suffix)

    def tearDown(self):
        del self.universe
        del self.tempdir

    def universe_from_tmp(self):
        return mda.Universe(self.outfile, convert_units=True)

    def test_write_atoms(self):
        self.universe.atoms.write(self.outfile)
        u2 = self.universe_from_tmp()
        assert_array_almost_equal(
            self.universe.atoms.positions, u2.atoms.positions,
            self.precision,
            err_msg=("atom coordinate mismatch between original and {0!s} file"
                     "".format(self.ext)))

    def test_write_empty_atomgroup(self):
        sel = self.universe.select_atoms('name doesntexist')
        assert_raises(IndexError, sel.write, self.outfile)

    def test_write_selection(self):
        CA = self.universe.select_atoms('name CA')
        CA.write(self.outfile)
        u2 = self.universe_from_tmp()
        # check EVERYTHING, otherwise we might get false positives!
        CA2 = u2.atoms
        assert_equal(len(u2.atoms), len(CA.atoms),
                     "written CA selection does not match original selection")
        assert_almost_equal(
            CA2.positions, CA.positions, self.precision,
            err_msg="CA coordinates do not agree with original")

    def test_write_Residue(self):
        G = self.universe.s4AKE.ARG[-2].atoms  # 2nd but last Arg
        G.write(self.outfile)
        u2 = self.universe_from_tmp()
        # check EVERYTHING, otherwise we might get false positives!
        G2 = u2.atoms
        assert_equal(
            len(u2.atoms), len(G.atoms),
            "written R206 Residue does not match original ResidueGroup")
        assert_almost_equal(
            G2.positions, G.positions, self.precision,
            err_msg="Residue R206 coordinates do not agree with original")

    def test_write_ResidueGroup(self):
        G = self.universe.s4AKE.LEU.atoms
        G.write(self.outfile)
        u2 = self.universe_from_tmp()
        G2 = u2.atoms
        assert_equal(
            len(u2.atoms), len(G.atoms),
            "written LEU ResidueGroup does not match original ResidueGroup")
        assert_almost_equal(
            G2.positions, G.positions, self.precision,
            err_msg="ResidueGroup LEU coordinates do not agree with original")

    def test_write_Segment(self):
        G = self.universe.s4AKE.atoms
        G.write(self.outfile)
        u2 = self.universe_from_tmp()
        G2 = u2.atoms
        assert_equal(len(u2.atoms), len(G.atoms),
                     "written s4AKE segment does not match original segment")
        assert_almost_equal(
            G2.positions, G.positions, self.precision,
            err_msg="segment s4AKE coordinates do not agree with original")

    def test_write_Universe(self):
        U = self.universe
        with mda.Writer(self.outfile) as W:
            W.write(U)
        u2 = self.universe_from_tmp()
        assert_equal(
            len(u2.atoms), len(U.atoms),
            "written 4AKE universe does not match original universe in size")
        assert_almost_equal(
            u2.atoms.positions, U.atoms.positions, self.precision,
            err_msg=("written universe 4AKE coordinates do not"
                     " agree with original"))


class TestWritePDB(_WriteAtoms):
    ext = "pdb"
    precision = 3


class TestWriteGRO(_WriteAtoms):
    ext = "gro"
    precision = 2

    def test_flag_convert_length(self):
        assert_equal(mda.core.flags['convert_lengths'], True,
                     "The flag convert_lengths SHOULD be True by default! "
                     "(If it is not then this might indicate a race condition"
                     " in the testing suite.)")


class TestAtomGroupTransformations(object):
    def setUp(self):
        self.u = mda.Universe(PSF, DCD)
        self.coords = self.u.atoms.positions.copy()
        self.cog = self.u.atoms.center_of_geometry()

    def tearDown(self):
        del self.u

    def test_translate(self):
        disp = np.ones(3)
        ag = self.u.atoms.translate(disp)
        assert_equal(ag, self.u.atoms)

        cog = self.u.atoms.center_of_geometry()
        diff = cog - self.cog
        assert_array_almost_equal(diff, disp, decimal=5)

    def test_rotate(self):
        # ensure that selection isn't centered at 0, 0, 0
        self.u.atoms.translate((1, 1, 1))
        self.coords += 1

        # check identify does nothing
        R = np.eye(3)
        self.u.atoms.rotate(R)
        assert_array_almost_equal(self.u.atoms.positions, self.coords)

        # check default rotation center is at 0, 0, 0. Changing this center
        # will break an unpredictable amount of old code.
        ag = self.u.atoms[:2]
        ag.positions = np.array([[1, 0, 0], [-1, 0, 0]])
        ag.rotate(transformations.rotation_matrix(1, [0, 0, 1])[:3, :3])
        assert_array_almost_equal(ag.positions[0], [np.cos(1), np.sin(1), 0])

        # check general rotation cases
        vec = np.array([[1, 0, 0], [-1, 0, 0]])
        axis = np.array([0, 0, 1])
        for angle in np.linspace(0, np.pi):
            R = transformations.rotation_matrix(angle, axis)
            ag.positions = vec.copy()
            res_ag = ag.rotate(R[:3, :3])
            assert_equal(ag, res_ag)
            assert_array_almost_equal(ag.positions[0], [np.cos(angle),
                                                        np.sin(angle),
                                                        0])

            ag.positions = vec.copy()
            ag.rotate(R[:3, :3], vec[0])
            assert_array_almost_equal(ag.positions[0], vec[0])
            assert_array_almost_equal(ag.positions[1], [- 2 * np.cos(angle) + 1,
                                                        - 2 * np.sin(angle),
                                                        0])

    def test_rotateby(self):
        R = np.eye(3)
        self.u.atoms.rotate(R)
        assert_array_almost_equal(self.u.atoms.positions, self.coords)

        vec = np.array([[1, 0, 0], [-1, 0, 0]])
        axis = np.array([0, 0, 1])

        ag = self.u.atoms[:2]
        ag.positions = vec

        for angle in np.linspace(0, np.pi):
            ag.positions = vec.copy()
            # needs to be rotated about origin
            res_ag = ag.rotateby(np.rad2deg(angle), axis)
            assert_equal(res_ag, ag)
            assert_array_almost_equal(ag.positions[0], [np.cos(angle),
                                                        np.sin(angle),
                                                        0])

            ag.positions = vec.copy()
            ag.rotateby(np.rad2deg(angle), axis, point=vec[0])
            assert_array_almost_equal(ag.positions[0], vec[0])
            assert_array_almost_equal(ag.positions[1], [- 2 * np.cos(angle) + 1,
                                                        - 2 * np.sin(angle),
                                                        0])

    def test_transform_rotation_only(self):
        R = np.eye(3)
        self.u.atoms.rotate(R)
        assert_array_almost_equal(self.u.atoms.positions, self.coords)

        vec = np.array([[1, 0, 0], [-1, 0, 0]])
        axis = np.array([0, 0, 1])

        ag = self.u.atoms[:2]
        ag.positions = vec

        for angle in np.linspace(0, np.pi):
            R = transformations.rotation_matrix(angle, axis)
            ag.positions = vec.copy()
            ag.transform(R)
            assert_array_almost_equal(ag.positions[0], [np.cos(angle),
                                                        np.sin(angle),
                                                        0])

    def test_transform_translation_only(self):
        disp = np.ones(3)
        T = np.eye(4)
        T[:3, 3] = disp
        ag = self.u.atoms.transform(T)
        assert_equal(ag, self.u.atoms)
        cog = self.u.atoms.center_of_geometry()
        diff = cog - self.cog
        assert_array_almost_equal(diff, disp, decimal=5)

    def test_transform_translation_and_rotation(self):
        angle = np.pi / 4
        axis = [0, 0, 1]
        disp = np.ones(3)
        T = transformations.rotation_matrix(angle, axis)
        T[:3, 3] = disp

        ag = self.u.atoms[:2]
        ag.positions = [[1, 0, 0], [-1, 0, 0]]
        ag.transform(T)

        assert_array_almost_equal(ag.positions[0], [np.cos(angle) + 1,
                                                    np.sin(angle) + 1,
                                                    1])

class TestCenter(object):
    def setUp(self):
        self.u = make_Universe(trajectory=True)
        self.ag = self.u.atoms[10:30]

    def tearDown(self):
        del self.u
        del self.ag

    def test_center_1(self):
        weights = np.zeros(self.ag.n_atoms)
        weights[0] = 1
        assert_array_almost_equal(self.ag.center(weights),
                                  self.ag.positions[0])

    def test_center_2(self):
        weights = np.zeros(self.ag.n_atoms)
        weights[:4] = 1. / 4.
        assert_array_almost_equal(self.ag.center(weights),
                                  self.ag.positions[:4].mean(axis=0))

    def test_center_wrong_length(self):
        weights = np.ones(self.ag.n_atoms + 4)

        assert_raises(ValueError, self.ag.center, weights)

    def test_center_wrong_shape(self):
        weights = np.ones((self.ag.n_atoms, 2))

        assert_raises(TypeError, self.ag.center, weights)


class TestSplit(object):
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)

    def tearDown(self):
        del self.universe
    
    def test_split_atoms(self):
        ag = self.universe.select_atoms("resid 1:50 and not resname LYS and "
                                        "(name CA or name CB)")
        sg = ag.split("atom")
        assert_equal(len(sg), len(ag))
        for g, ref_atom in zip(sg, ag):
            atom = g[0]
            assert_equal(len(g), 1)
            assert_equal(atom.name, ref_atom.name)
            assert_equal(atom.resid, ref_atom.resid)

    def test_split_residues(self):
        ag = self.universe.select_atoms("resid 1:50 and not resname LYS and "
                                        "(name CA or name CB)")
        sg = ag.split("residue")
        assert_equal(len(sg), len(ag.residues.resids))
        for g, ref_resname in zip(sg, ag.residues.resnames):
            if ref_resname == "GLY":
                assert_equal(len(g), 1)
            else:
                assert_equal(len(g), 2)
            for atom in g:
                assert_equal(atom.resname, ref_resname)

    def test_split_segments(self):
        ag = self.universe.select_atoms("resid 1:50 and not resname LYS and "
                                        "(name CA or name CB)")
        sg = ag.split("segment")
        assert_equal(len(sg), len(ag.segments.segids))
        for g, ref_segname in zip(sg, ag.segments.segids):
            for atom in g:
                assert_equal(atom.segid, ref_segname)    

    def test_split_VE(self):
        ag = self.universe.atoms[:40]

        assert_raises(ValueError, ag.split, 'something')


class TestWrap(object):
    @dec.skipif(parser_not_found('TRZ'),
                'TRZ parser not available. Are you using python 3?')
    def setUp(self):
        self.u = mda.Universe(TRZ_psf, TRZ)
        self.ag = self.u.atoms[:100]

    def tearDown(self):
        del self.u
        del self.ag

    def test_wrap_comp_fail(self):
        assert_raises(ValueError, self.ag.wrap, compound='strawberries')

    def test_wrap_cent_fail(self):
        assert_raises(ValueError, self.ag.wrap, compound='residues', center='avacado')

    def test_wrap_box_fail(self):
        assert_raises(ValueError, self.ag.wrap, box=np.array([0, 1]))

    def _in_box(self, coords):
        """Check that a set of coordinates are 0.0 <= r <= box"""
        box = self.u.dimensions[:3]

        return (coords >= 0.0).all() and (coords <= box).all()

    def test_wrap_atoms(self):
        ag = self.u.atoms[100:200]
        ag.wrap(compound='atoms')

        assert_equal(self._in_box(ag.positions), True)

    def test_wrap_group(self):
        ag = self.u.atoms[:100]
        ag.wrap(compound='group')

        cen = ag.center_of_mass()

        assert_equal(self._in_box(cen), True)

    def test_wrap_residues(self):
        ag = self.u.atoms[300:400]
        ag.wrap(compound='residues')

        cen = np.vstack([r.atoms.center_of_mass() for r in ag.residues])

        assert_equal(self._in_box(cen), True)

    def test_wrap_segments(self):
        ag = self.u.atoms[1000:1200]
        ag.wrap(compound='segments')

        cen = np.vstack([s.atoms.center_of_mass() for s in ag.segments])

        assert_equal(self._in_box(cen), True)

    def test_wrap_fragments(self):
        ag = self.u.atoms[:250]
        ag.wrap(compound='fragments')

        cen = np.vstack([f.center_of_mass() for f in ag.fragments])

        assert_equal(self._in_box(cen), True)


class TestAtomGroupProperties(object):
    """Test working with the properties of Atoms via AtomGroups

    Check that:
    - getting properties from AG matches the Atom values
    - setting properties from AG changes the Atom
    - setting the property on Atom changes AG
    """
    @staticmethod
    def get_new(att_type):
        """Return enough values to change the small g"""
        if att_type == 'string':
            return ['A', 'B', 'C', 'D', 'E', 'F']
        elif att_type == 'float':
            return np.array([0.001, 0.002, 0.003, 0.005, 0.012, 0.025], dtype=np.float32)
        elif att_type == 'int':
            return [4, 6, 8, 1, 5, 4]

    def _check_ag_matches_atom(self, att, atts, ag):
        """Checking Atomgroup property matches Atoms"""
        # Check that accessing via AtomGroup is identical to doing
        # a list comprehension over AG
        ref = [getattr(atom, att) for atom in ag]

        assert_equal(ref, getattr(ag, atts),
                     err_msg="AtomGroup doesn't match Atoms for property: {0}".format(att))

    def _change_atom_check_ag(self, att, atts, vals, ag):
        """Changing Atom, checking AtomGroup matches this"""
        # Set attributes via Atoms
        for atom, val in zip(ag, vals):
            setattr(atom, att, val)
        # Check that AtomGroup returns new values
        other = getattr(ag, atts)

        assert_equal(vals, other,
                     err_msg="Change to Atoms not reflected in AtomGroup for property: {0}".format(att))

    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_attributes(self):
        u = make_Universe(('names', 'resids', 'segids', 'types', 'altLocs',
                           'charges', 'masses', 'radii', 'bfactors',
                           'occupancies'))
        u.atoms.occupancies = 1.0
        master = u.atoms
        idx = [0, 1, 4, 7, 11, 14]
        ag = master[idx]

        for att, atts, att_type in (
                ('name', 'names', 'string'),
                ('type', 'types', 'string'),
                ('altLoc', 'altLocs', 'string'),
                ('charge', 'charges', 'float'),
                ('mass', 'masses', 'float'),
                ('radius', 'radii', 'float'),
                ('bfactor', 'bfactors', 'float'),
                ('occupancy', 'occupancies', 'float')
        ):
            vals = self.get_new(att_type)
            yield self._check_ag_matches_atom, att, atts, ag
            yield self._change_atom_check_ag, att, atts, vals, ag


class TestOrphans(object):
    """Test moving Universes out of scope and having A/AG persist

    Atoms and AtomGroups from other scopes should work, namely:
      - should have access to Universe
      - should be able to use the Reader (coordinates)
    """
    def test_atom(self):
        u = mda.Universe(two_water_gro)

        def getter():
            u2 = mda.Universe(two_water_gro)
            return u2.atoms[1]

        atom = getter()

        assert_(atom is not u.atoms[1])
        assert_(len(atom.universe.atoms) == len(u.atoms))
        assert_array_almost_equal(atom.position, u.atoms[1].position)

    def test_atomgroup(self):
        u = mda.Universe(two_water_gro)

        def getter():
            u2 = mda.Universe(two_water_gro)
            return u2.atoms[:4]

        ag = getter()
        ag2 = u.atoms[:4]
        assert_(ag is not ag2)
        assert_(len(ag.universe.atoms) == len(u.atoms))
        assert_array_almost_equal(ag.positions, ag2.positions)


class TestCrossUniverse(object):
    """Test behaviour when we mix Universes"""
    def _check_badadd(self, a, b):
        def add(x, y):
            return x + y
        assert_raises(ValueError, add, a, b)

    def test_add_mixed_universes(self):
        # Issue #532
        # Checks that adding objects from different universes
        # doesn't proceed quietly.
        u1 = mda.Universe(two_water_gro)
        u2 = mda.Universe(two_water_gro)

        A = [u1.atoms[:2], u1.atoms[3]]
        B = [u2.atoms[:3], u2.atoms[0]]

        # Checks Atom to Atom, Atom to AG, AG to Atom and AG to AG
        for x, y in itertools.product(A, B):
            yield self._check_badadd, x, y

    def test_adding_empty_ags(self):
        # Check that empty AtomGroups don't trip up on the Universe check
        u = mda.Universe(two_water_gro)

        assert_(len(u.atoms[[]] + u.atoms[:3]) == 3)
        assert_(len(u.atoms[:3] + u.atoms[[]]) == 3)


class TestDihedralSelections(object):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        self.dih_prec = 2

    def tearDown(self):
        del self.universe
        del self.dih_prec

    def test_phi_selection(self):
        phisel = self.universe.s4AKE.r10.phi_selection()
        assert_equal(phisel.names, ['C', 'N', 'CA', 'C'])
        assert_equal(phisel.residues.resids, [9, 10])
        assert_equal(phisel.residues.resnames, ['PRO', 'GLY'])

    def test_psi_selection(self):
        psisel = self.universe.s4AKE.r10.psi_selection()
        assert_equal(psisel.names, ['N', 'CA', 'C', 'N'])
        assert_equal(psisel.residues.resids, [10, 11])
        assert_equal(psisel.residues.resnames, ['GLY', 'ALA'])

    def test_omega_selection(self):
        osel = self.universe.s4AKE.r8.omega_selection()
        assert_equal(osel.names, ['CA', 'C', 'N', 'CA'])
        assert_equal(osel.residues.resids, [8, 9])
        assert_equal(osel.residues.resnames, ['ALA', 'PRO'])

    def test_chi1_selection(self):
        sel = self.universe.s4AKE.r13.chi1_selection()  # LYS
        assert_equal(sel.names, ['N', 'CA', 'CB', 'CG'])
        assert_equal(sel.residues.resids, [13])
        assert_equal(sel.residues.resnames, ['LYS'])

    def test_phi_sel_fail(self):
        sel = self.universe.residues[0].phi_selection()
        assert_equal(sel, None)

    def test_psi_sel_fail(self):
        sel = self.universe.residues[-1].psi_selection()
        assert_equal(sel, None)

    def test_omega_sel_fail(self):
        sel = self.universe.residues[-1].omega_selection()
        assert_equal(sel, None)

    def test_ch1_sel_fail(self):
        sel = self.universe.s4AKE.r8.chi1_selection()
        assert_equal(sel, None)  # ALA

    def test_dihedral_phi(self):
        u = self.universe
        phisel = u.s4AKE.r10.phi_selection()
        assert_almost_equal(phisel.dihedral.value(), -168.57384, self.dih_prec)

    def test_dihedral_psi(self):
        u = self.universe
        psisel = u.s4AKE.r10.psi_selection()
        assert_almost_equal(psisel.dihedral.value(), -30.064838, self.dih_prec)

    def test_dihedral_omega(self):
        u = self.universe
        osel = u.s4AKE.r8.omega_selection()
        assert_almost_equal(osel.dihedral.value(), -179.93439, self.dih_prec)

    def test_dihedral_chi1(self):
        u = self.universe
        sel = u.s4AKE.r13.chi1_selection()  # LYS
        assert_almost_equal(sel.dihedral.value(), -58.428127, self.dih_prec)


class TestPBCFlag(object):
    @dec.skipif(parser_not_found('TRZ'),
                'TRZ parser not available. Are you using python 3?')
    def setUp(self):
        self.prec = 3
        self.universe = mda.Universe(TRZ_psf, TRZ)
        self.ref_noPBC = {
            'COG': np.array([4.23789883, 0.62429816, 2.43123484], dtype=np.float32),
            'COM': np.array([4.1673783, 0.70507009, 2.21175832]),
            'ROG': 119.30368949900134, 'Shape': 0.6690026954813445,
            'Asph': 0.5305456387833748,
            'MOI': np.array([
                [152117.06620921, 55149.54042136, -26630.46034023],
                [55149.54042136, 72869.64061494, 21998.1778074],
                [-26630.46034023, 21998.1778074, 162388.70002471]]),
            'BBox': np.array([[-75.74159241, -144.86634827, -94.47974396], [95.83090973, 115.11561584, 88.09812927]],
                                dtype=np.float32),
            'BSph': (173.40482, np.array([4.23789883, 0.62429816, 2.43123484], dtype=np.float32)),
            'PAxes': np.array([
                [-0.78787867, -0.26771575, 0.55459488],
                [0.40611024, 0.45112859, 0.7947059],
                [0.46294889, -0.85135849, 0.24671249]])
        }
        self.ref_PBC = {
            'COG': np.array([26.82960892, 31.5592289, 30.98238945], dtype=np.float32),
            'COM': np.array([26.67781143, 31.2104336, 31.19796289]),
            'ROG': 27.713008969174918, 'Shape': 0.0017390512580463542,
            'Asph': 0.020601215358731016,
            'MOI': np.array([
                [7333.79167791, -211.8997285, -721.50785456],
                [-211.8997285, 7059.07470427, -91.32156884],
                [-721.50785456, -91.32156884, 6509.31735029]]),
            'BBox': np.array(
                [[1.45964116e-01, 1.85623169e-02, 4.31785583e-02], [5.53314018e+01, 5.54227829e+01, 5.54158211e+01]],
                dtype=np.float32),
            'BSph': (47.923367, np.array([26.82960892, 31.5592289, 30.98238945], dtype=np.float32)),
            'PAxes': np.array([
                [-0.85911708, 0.19258726, 0.4741603],
                [-0.07520116, -0.96394227, 0.25526473],
                [-0.50622389, -0.18364489, -0.84262206]])
        }
        self.ag = self.universe.residues[0:3]

    def tearDown(self):
        mda.core.flags['use_pbc'] = False
        del self.universe
        del self.ref_noPBC
        del self.ref_PBC
        del self.ag

    def test_flag(self):
        # Test default setting of flag
        assert_equal(mda.core.flags['use_pbc'], False)

    def test_default(self):
        # Test regular behaviour
        assert_almost_equal(self.ag.center_of_geometry(), self.ref_noPBC['COG'], self.prec)
        assert_almost_equal(self.ag.center_of_mass(), self.ref_noPBC['COM'], self.prec)
        assert_almost_equal(self.ag.radius_of_gyration(), self.ref_noPBC['ROG'], self.prec)
        assert_almost_equal(self.ag.shape_parameter(), self.ref_noPBC['Shape'], self.prec)
        assert_almost_equal(self.ag.asphericity(), self.ref_noPBC['Asph'], self.prec)
        assert_almost_equal(self.ag.moment_of_inertia(), self.ref_noPBC['MOI'], self.prec)
        assert_almost_equal(self.ag.bbox(), self.ref_noPBC['BBox'], self.prec)
        assert_almost_equal(self.ag.bsphere()[0], self.ref_noPBC['BSph'][0], self.prec)
        assert_almost_equal(self.ag.bsphere()[1], self.ref_noPBC['BSph'][1], self.prec)
        assert_almost_equal(self.ag.principal_axes(), self.ref_noPBC['PAxes'], self.prec)

    def test_pbcflag(self):
        # Test using ag method flag
        assert_almost_equal(self.ag.center_of_geometry(pbc=True), self.ref_PBC['COG'], self.prec)
        assert_almost_equal(self.ag.center_of_mass(pbc=True), self.ref_PBC['COM'], self.prec)
        assert_almost_equal(self.ag.radius_of_gyration(pbc=True), self.ref_PBC['ROG'], self.prec)
        assert_almost_equal(self.ag.shape_parameter(pbc=True), self.ref_PBC['Shape'], self.prec)
        assert_almost_equal(self.ag.asphericity(pbc=True), self.ref_PBC['Asph'], self.prec)
        assert_almost_equal(self.ag.moment_of_inertia(pbc=True), self.ref_PBC['MOI'], self.prec)
        assert_almost_equal(self.ag.bbox(pbc=True), self.ref_PBC['BBox'], self.prec)
        assert_almost_equal(self.ag.bsphere(pbc=True)[0], self.ref_PBC['BSph'][0], self.prec)
        assert_almost_equal(self.ag.bsphere(pbc=True)[1], self.ref_PBC['BSph'][1], self.prec)
        assert_almost_equal(self.ag.principal_axes(pbc=True), self.ref_PBC['PAxes'], self.prec)

    def test_usepbc_flag(self):
        # Test using the core.flags flag
        mda.core.flags['use_pbc'] = True
        assert_almost_equal(self.ag.center_of_geometry(), self.ref_PBC['COG'], self.prec)
        assert_almost_equal(self.ag.center_of_mass(), self.ref_PBC['COM'], self.prec)
        assert_almost_equal(self.ag.radius_of_gyration(), self.ref_PBC['ROG'], self.prec)
        assert_almost_equal(self.ag.shape_parameter(), self.ref_PBC['Shape'], self.prec)
        assert_almost_equal(self.ag.asphericity(), self.ref_PBC['Asph'], self.prec)
        assert_almost_equal(self.ag.moment_of_inertia(), self.ref_PBC['MOI'], self.prec)
        assert_almost_equal(self.ag.bbox(), self.ref_PBC['BBox'], self.prec)
        assert_almost_equal(self.ag.bsphere()[0], self.ref_PBC['BSph'][0], self.prec)
        assert_almost_equal(self.ag.bsphere()[1], self.ref_PBC['BSph'][1], self.prec)
        assert_almost_equal(self.ag.principal_axes(), self.ref_PBC['PAxes'], self.prec)
        mda.core.flags['use_pbc'] = False

    def test_override_flag(self):
        # Test using the core.flags flag, then overriding
        mda.core.flags['use_pbc'] = True
        assert_almost_equal(self.ag.center_of_geometry(pbc=False), self.ref_noPBC['COG'], self.prec)
        assert_almost_equal(self.ag.center_of_mass(pbc=False), self.ref_noPBC['COM'], self.prec)
        assert_almost_equal(self.ag.radius_of_gyration(pbc=False), self.ref_noPBC['ROG'], self.prec)
        assert_almost_equal(self.ag.shape_parameter(pbc=False), self.ref_noPBC['Shape'], self.prec)
        assert_almost_equal(self.ag.asphericity(pbc=False), self.ref_noPBC['Asph'], self.prec)
        assert_almost_equal(self.ag.moment_of_inertia(pbc=False), self.ref_noPBC['MOI'], self.prec)
        assert_almost_equal(self.ag.bbox(pbc=False), self.ref_noPBC['BBox'], self.prec)
        assert_almost_equal(self.ag.bsphere(pbc=False)[0], self.ref_noPBC['BSph'][0], self.prec)
        assert_almost_equal(self.ag.bsphere(pbc=False)[1], self.ref_noPBC['BSph'][1], self.prec)
        assert_almost_equal(self.ag.principal_axes(pbc=False), self.ref_noPBC['PAxes'], self.prec)
        mda.core.flags['use_pbc'] = False


@dec.skipif(parser_not_found('DCD'),
            'DCD parser not available. Are you using python 3?')
def test_instantselection_termini():
    """Test that instant selections work, even for residues that are also termini (Issue 70)"""
    universe = mda.Universe(PSF, DCD)
    assert_equal(universe.residues[20].CA.name, 'CA', "CA of MET21 is not selected correctly")
    del universe


class TestAtomGroup(object):
    """Tests of AtomGroup; selections are tested separately.

    These are from before the big topology rework (aka #363) but are still valid.
    There is likely lots of duplication between here and other tests.
    """
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = mda.Universe(PSF, DCD)
        self.ag = self.universe.atoms  # prototypical AtomGroup
        self.dih_prec = 2

    def test_getitem_int(self):
        assert_equal(self.universe.atoms[0].ix, self.universe.atoms.ix[0])

    def test_getitem_slice(self):
        assert_array_equal(self.universe.atoms[0:4].ix,
                           self.universe.atoms.ix[:4])

    def test_getitem_slice2(self):
        assert_equal(self.universe.atoms[0:8:2].ix,
                     self.universe.atoms.ix[0:8:2])

    def test_getitem_str(self):
        ag1 = self.universe.atoms['HT1']
        # select_atoms always returns an AtomGroup even if single result
        ag2 = self.universe.select_atoms('name HT1')[0]
        assert_equal(ag1, ag2)

    def test_getitem_IE(self):
        d = {'A': 1}
        assert_raises(IndexError, self.universe.atoms.__getitem__, d)

    def test_bad_make(self):
        assert_raises(TypeError, mda.core.groups.AtomGroup, ['these', 'are', 'not', 'atoms'])

    def test_n_atoms(self):
        assert_equal(self.ag.n_atoms, 3341)

    def test_n_residues(self):
        assert_equal(self.ag.n_residues, 214)

    def test_n_segments(self):
        assert_equal(self.ag.n_segments, 1)

    def test_resids_dim(self):
        assert_equal(len(self.ag.resids), len(self.ag))

    def test_resnums_dim(self):
        assert_equal(len(self.ag.resnums), len(self.ag))

    def test_segids_dim(self):
        assert_equal(len(self.ag.segids), len(self.ag))

    def test_len(self):
        """testing that len(atomgroup) == atomgroup.n_atoms"""
        assert_equal(len(self.ag), self.ag.n_atoms, "len and n_atoms disagree")

    def test_center_of_geometry(self):
        assert_array_almost_equal(self.ag.center_of_geometry(),
                                  np.array([-0.04223963, 0.0141824, -0.03505163], dtype=np.float32))

    def test_center_of_mass(self):
        assert_array_almost_equal(self.ag.center_of_mass(),
                                  np.array([-0.01094035, 0.05727601, -0.12885778]))

    def test_coordinates(self):
        assert_array_almost_equal(
            self.ag.positions[1000:2000:200],
            np.array([[3.94543672, -12.4060812, -7.26820087],
                      [13.21632767, 5.879035, -14.67914867],
                      [12.07735443, -9.00604534, 4.09301519],
                      [11.35541916, 7.0690732, -0.32511973],
                      [-13.26763439, 4.90658951, 10.6880455]],
                     dtype=np.float32))

    def test_principal_axes(self):
        assert_array_almost_equal(
            self.ag.principal_axes(),
            np.array([[1.53389276e-03, 4.41386224e-02, 9.99024239e-01],
                      [1.20986911e-02, 9.98951474e-01, -4.41539838e-02],
                      [-9.99925632e-01, 1.21546132e-02, 9.98264877e-04],]))

    def test_total_charge(self):
        assert_almost_equal(self.ag.total_charge(), -4.0, decimal=4)

    def test_total_mass(self):
        assert_almost_equal(self.ag.total_mass(), 23582.043)

    def test_indices_ndarray(self):
        assert_equal(isinstance(self.ag.indices, np.ndarray), True)

    def test_indices(self):
        # pylint: disable=unsubscriptable-object
        assert_array_equal(self.ag.indices[:5], np.array([0, 1, 2, 3, 4]))

    def test_resids_ndarray(self):
        assert_equal(isinstance(self.ag.resids, np.ndarray), True)

    def test_resids(self):
        assert_array_equal(self.ag.residues.resids, np.arange(1, 215))

    def test_resnums_ndarray(self):
        assert_equal(isinstance(self.ag.residues.resnums, np.ndarray), True)

    def test_resnums(self):
        assert_array_equal(self.ag.residues.resnums, np.arange(1, 215))

    def test_resnames_ndarray(self):
        assert_equal(isinstance(self.ag.residues.resnames, np.ndarray), True)

    def test_resnames(self):
        resnames = self.ag.residues.resnames
        assert_array_equal(resnames[0:3], np.array(["MET", "ARG", "ILE"]))

    def test_names_ndarray(self):
        assert_equal(isinstance(self.ag.names, np.ndarray), True)

    def test_names(self):
        names = self.ag.names
        assert_array_equal(names[0:3], np.array(["N", "HT1", "HT2"]))

    def test_segids_ndarray(self):
        assert_equal(isinstance(self.ag.segids, np.ndarray), True)

    def test_segids(self):
        segids = self.ag.segids
        assert_array_equal(segids[0], np.array(["4AKE"]))

    def test_masses_ndarray(self):
        assert_equal(isinstance(self.ag.masses, np.ndarray), True)

    def test_masses(self):
        masses = self.ag.masses
        assert_array_equal(masses[0:3], np.array([14.007, 1.008, 1.008]))

    def test_charges_ndarray(self):
        assert_equal(isinstance(self.ag.charges, np.ndarray), True)

    def test_charges(self):
        assert_array_almost_equal(self.ag.charges[1000:2000:200],
                                  np.array([-0.09, 0.09, -0.47, 0.51, 0.09]))

    def test_bad_add_AG(self):
        def bad_add():
            return self.ag + [1, 2, 3]
        assert_raises(TypeError, bad_add)

    def test_bool_false(self):
        # Issue #304
        ag = self.universe.atoms[[]]
        assert_equal(bool(ag), False)

    def test_bool_true(self):
        # Issue #304
        ag = self.universe.atoms[:3]
        assert_equal(bool(ag), True)

    def test_repr(self):
        # Should make sure that the user facing info stays as expected
        assert_equal(repr(self.ag), "<AtomGroup with 3341 atoms>")

    def test_set_resnum_single(self):
        ag = self.universe.atoms[:3]
        ag.residues.resnums = 5
        for at in ag:
            assert_equal(at.resnum, 5)
        assert_equal(all(ag.resnums == 5), True)

    def test_set_resname_single(self):
        ag = self.universe.atoms[:3]
        new = 'abc'
        ag.residues.resnames = new
        for at in ag:
            assert_equal(at.resname, new)
        assert_equal(all(ag.resnames == new), True)

    def test_packintobox_badshape(self):
        ag = self.universe.atoms[:10]
        box = np.zeros(9, dtype=np.float32).reshape(3, 3)
        assert_raises(ValueError, ag.pack_into_box, box=box)

    def test_packintobox_noshape(self):
        ag = self.universe.atoms[:10]
        assert_raises(ValueError, ag.pack_into_box)

    def test_packintobox(self):
        """test AtomGroup.pack_into_box(): Tests application of periodic boundary
        conditions on coordinates

        Reference system doesn't have dimensions, so an arbitrary box is
        imposed on the system

        """
        u = self.universe
        u.trajectory.rewind()  # just to make sure...
        ag = u.atoms[1000:2000:200]
        # Provide arbitrary box
        ag.pack_into_box(box=np.array([5., 5., 5.], dtype=np.float32))
        assert_array_almost_equal(
            ag.positions,
            np.array([[3.94543672, 2.5939188, 2.73179913],
                      [3.21632767, 0.879035, 0.32085133],
                      [2.07735443, 0.99395466, 4.09301519],
                      [1.35541916, 2.0690732, 4.67488003],
                      [1.73236561, 4.90658951, 0.6880455]], dtype=np.float32))

    def test_residues(self):
        u = self.universe
        assert_equal(u.residues[100].atoms.ix,
                     u.select_atoms('resname ILE and resid 101').atoms.ix,
                     "Direct selection from residue group does not match "
                     "expected I101.")

    def test_segments(self):
        u = self.universe
        assert_equal(len(u.segments.s4AKE.atoms),
                     len(u.select_atoms('segid 4AKE').atoms),
                     "Direct selection of segment 4AKE from segments failed.")

    def test_index_integer(self):
        u = self.universe
        a = u.atoms[100]
        assert_(isinstance(a, mda.core.groups.Atom), "integer index did not return Atom")

    def test_index_slice(self):
        u = self.universe
        a = u.atoms[100:200:10]
        assert_(isinstance(a, mda.core.groups.AtomGroup),
                "slice index did not return AtomGroup")

    def test_index_slice_empty(self):
        u = self.universe
        assert_array_equal(u.atoms[0:0], [],
                           "making an empty AtomGroup failed")

    def test_index_advancedslice(self):
        u = self.universe
        aslice = [0, 10, 20, -1, 10]
        ag = u.atoms[aslice]
        assert_(isinstance(ag, mda.core.groups.AtomGroup),
                "advanced slicing does not produce a AtomGroup")
        assert_equal(ag[1], ag[-1], "advanced slicing does not preserve order")

    def test_boolean_indexing(self):
        # index an array with a sequence of bools
        # issue #282
        sel = np.array([True, False, True])
        ag = self.universe.atoms[10:13]
        ag2 = ag[sel]
        assert_equal(len(ag2), 2)
        for at in [ag[0], ag[2]]:
            assert_equal(at in ag2, True)

    def test_boolean_indexing_2(self):
        # index an array with a sequence of bools
        # issue #282
        sel = [True, False, True]
        ag = self.universe.atoms[10:13]
        ag2 = ag[sel]
        assert_equal(len(ag2), 2)
        for at in [ag[0], ag[2]]:
            assert_equal(at in ag2, True)

    def test_bool_IE(self):
        # indexing with empty list doesn't run foul of bool check
        sel = []
        ag = self.universe.atoms[10:30]
        ag2 = ag[sel]
        assert_equal(len(ag2), 0)

    def test_dihedral_ValueError(self):
        """test that AtomGroup.dihedral() raises ValueError if not exactly
        4 atoms given"""
        nodih = self.universe.select_atoms("resid 3:10")
        assert_raises(ValueError, getattr, nodih, 'dihedral')
        nodih = self.universe.select_atoms("resid 3:5")
        assert_raises(ValueError, getattr, nodih, 'dihedral')

    def test_improper(self):
        u = self.universe
        u.trajectory.rewind()  # just to make sure...
        peptbond = u.select_atoms("atom 4AKE 20 C", "atom 4AKE 21 CA",
                                  "atom 4AKE 21 N", "atom 4AKE 21 HN")
        assert_almost_equal(peptbond.improper.value(), 168.52952575683594,
                            self.dih_prec,
                            "Peptide bond improper dihedral for M21 "
                            "calculated wrongly.")

    def test_dihedral_equals_improper(self):
        u = self.universe
        u.trajectory.rewind()  # just to make sure...
        peptbond = u.select_atoms("atom 4AKE 20 C", "atom 4AKE 21 CA",
                                  "atom 4AKE 21 N", "atom 4AKE 21 HN")
        assert_equal(peptbond.improper.value(), peptbond.dihedral.value(),
                     "improper() and proper dihedral() give different results")

    def test_bond(self):
        self.universe.trajectory.rewind()  # just to make sure...
        sel2 = self.universe.s4AKE.r98.atoms.select_atoms("name OE1", "name OE2")
        assert_almost_equal(sel2.bond.value(), 2.1210737228393555, 3,
                            "distance of Glu98 OE1--OE2 wrong")

    def test_bond_pbc(self):
        self.universe.trajectory.rewind()
        sel2 = self.universe.s4AKE.r98.atoms.select_atoms("name OE1", "name OE2")
        assert_almost_equal(sel2.bond.value(pbc=True), 2.1210737228393555, 3,
                            "distance of Glu98 OE1--OE2 wrong")

    def test_bond_ValueError(self):
        ag = self.universe.atoms[:4]
        assert_raises(ValueError, getattr, ag, 'bond')

    def test_angle(self):
        self.universe.trajectory.rewind()  # just to make sure...
        sel3 = self.universe.s4AKE.r98.atoms.select_atoms("name OE1", "name CD",
                                                          "name OE2")
        assert_almost_equal(sel3.angle.value(), 117.46187591552734, 3,
                            "angle of Glu98 OE1-CD-OE2 wrong")

    def test_angle_ValueError(self):
        ag = self.universe.atoms[:2]
        assert_raises(ValueError, getattr, ag, 'angle')

    def test_shape_parameter(self):
        s = self.universe.s4AKE.atoms.shape_parameter()
        assert_almost_equal(s, 0.00240753939086033, 6)

    def test_asphericity(self):
        a = self.universe.s4AKE.atoms.asphericity()
        assert_almost_equal(a, 0.020227504542775828, 6)

    def test_positions(self):
        ag = self.universe.select_atoms("bynum 12:42")
        pos = ag.positions + 3.14
        ag.positions = pos
        # should work
        assert_almost_equal(ag.positions, pos,
                            err_msg="failed to update atoms 12:42 position "
                            "to new position")

        def set_badarr(pos=pos):
            # create wrong size array
            badarr = np.random.random((pos.shape[0] - 1, pos.shape[1] - 1))
            ag.positions = badarr
        assert_raises(ValueError, set_badarr)

    def test_set_names(self):
        ag = self.universe.atoms[:2]
        names = ['One', 'Two']
        ag.names = names
        for a, b in zip(ag, names):
            assert_equal(a.name, b)

    @attr("issue")
    def test_nonexistent_instantselector_raises_AttributeError(self):
        def access_nonexistent_instantselector():
            self.universe.atoms.NO_SUCH_ATOM
        assert_raises(AttributeError, access_nonexistent_instantselector)

    def test_atom_order(self):
        assert_equal(self.universe.atoms.indices,
                     sorted(self.universe.atoms.indices))


class TestAtomGroupTimestep(object):
    """Tests the AtomGroup.ts attribute (partial timestep)"""

    @dec.skipif(parser_not_found('TRZ'),
                'TRZ parser not available. Are you using python 3?')
    def setUp(self):
        self.universe = mda.Universe(TRZ_psf, TRZ)
        self.prec = 6

    def tearDown(self):
        del self.universe
        del self.prec

    def test_partial_timestep(self):
        ag = self.universe.select_atoms('name Ca')
        idx = ag.indices

        assert_equal(len(ag.ts._pos), len(ag))

        for ts in self.universe.trajectory[0:20:5]:
            assert_array_almost_equal(ts.positions[idx], ag.ts.positions, self.prec,
                                      err_msg="Partial timestep coordinates wrong")
            assert_array_almost_equal(ts.velocities[idx], ag.ts.velocities, self.prec,
                                      err_msg="Partial timestep coordinates wrong")
