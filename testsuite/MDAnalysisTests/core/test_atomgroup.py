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
from MDAnalysisTests import tempdir


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
        self.u.atoms.translate(disp)
        cog = self.u.atoms.center_of_geometry()
        diff = cog - self.cog
        assert_array_almost_equal(diff, disp, decimal=5)

    def test_rotate(self):
        R = np.eye(3)
        self.u.atoms.rotate(R)
        assert_array_almost_equal(self.u.atoms.positions, self.coords)

        vec = np.array([1, 0, 0])
        axis = np.array([0, 0, 1])

        ag = self.u.atoms.select_atoms('bynum 1')
        ag.positions = vec

        for angle in np.linspace(0, np.pi):
            R = transformations.rotation_matrix(angle, axis)
            ag.positions = vec.copy()
            ag.rotate(R[:3, :3])
            assert_array_almost_equal(ag.positions[0], [np.cos(angle),
                                                        np.sin(angle),
                                                        0])

    def test_rotateby(self):
        R = np.eye(3)
        self.u.atoms.rotate(R)
        assert_array_almost_equal(self.u.atoms.positions, self.coords)

        vec = np.array([1, 0, 0])
        axis = np.array([0, 0, 1])

        ag = self.u.atoms.select_atoms('bynum 1')
        ag.positions = vec

        for angle in np.linspace(0, np.pi):
            ag.positions = vec.copy()
            # needs to be rotated about origin
            ag.rotateby(np.rad2deg(angle), axis, [0, 0, 0])
            assert_array_almost_equal(ag.positions[0], [np.cos(angle),
                                                        np.sin(angle),
                                                        0])

    def test_transform(self):
        R = np.eye(3)
        self.u.atoms.rotate(R)
        assert_array_almost_equal(self.u.atoms.positions, self.coords)

        vec = np.array([1, 0, 0])
        axis = np.array([0, 0, 1])

        ag = self.u.atoms.select_atoms('bynum 1')
        ag.positions = vec

        for angle in np.linspace(0, np.pi):
            R = transformations.rotation_matrix(angle, axis)
            ag.positions = vec.copy()
            ag.transform(R)
            assert_array_almost_equal(ag.positions[0], [np.cos(angle),
                                                        np.sin(angle),
                                                        0])
