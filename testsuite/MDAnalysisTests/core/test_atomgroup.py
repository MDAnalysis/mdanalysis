# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

from glob import glob
from os import path
import numpy as np
import warnings

from numpy.testing import (
    assert_,
    assert_raises,
    assert_equal,
    assert_array_almost_equal,
    raises
)


import MDAnalysis as mda
from MDAnalysis.lib import transformations
from MDAnalysis.core.topologyobjects import (
    Bond,
    Angle,
    Dihedral,
    ImproperDihedral,
)

from MDAnalysisTests.datafiles import (PSF, DCD)
from MDAnalysisTests.core.groupbase import make_Universe
from MDAnalysisTests import tempdir

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
