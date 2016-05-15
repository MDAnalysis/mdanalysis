# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

import MDAnalysis
import numpy as np
from numpy.testing import (
    assert_equal, assert_almost_equal, TestCase,
    assert_array_equal,
)
from nose.plugins.attrib import attr

from MDAnalysis.tests.datafiles import GRO_velocity, PDB_xvf, TRR_xvf


class TestAtom_ForceVelocity(TestCase):
    def setUp(self):
        self.u = MDAnalysis.Universe(PDB_xvf, TRR_xvf)
        self.a = self.u.atoms[0]

    def tearDown(self):
        del self.u
        del self.a

    def test_atom_force_get(self):
        assert_equal(self.a.force, self.u.atoms.forces[0])

    def test_atom_velocity_get(self):
        assert_equal(self.a.velocity, self.u.atoms.velocities[0])

    def test_atom_force_set(self):
        ref = np.arange(3)
        self.a.force = ref

        assert_equal(self.a.force, ref)
        assert_equal(self.u.atoms.forces[0], ref)

    def test_atom_velocity_set(self):
        ref = np.arange(3)
        self.a.velocity = ref

        assert_equal(self.a.velocity, ref)
        assert_equal(self.u.atoms.velocities[0], ref)

    def test_pos_iteration(self):
        ag = self.u.atoms[[0]]

        val = np.array([self.a.position for ts in self.u.trajectory])
        ref = np.array([ag.positions[0] for ts in self.u.trajectory])

        assert_array_equal(val, ref)

    def test_vel_iteration(self):
        ag = self.u.atoms[[0]]

        val = np.array([self.a.velocity for ts in self.u.trajectory])
        ref = np.array([ag.velocities[0] for ts in self.u.trajectory])

        assert_array_equal(val, ref)

    def test_for_iteration(self):
        ag = self.u.atoms[[0]]

        val = np.array([self.a.force for ts in self.u.trajectory])
        ref = np.array([ag.forces[0] for ts in self.u.trajectory])

        assert_array_equal(val, ref)


class TestGROVelocities(TestCase):
    def setUp(self):
        #reference velocities for the full 6-atom test case:
        self.reference_velocities = np.array(
            [[-101.227, -0.57999998, 0.43400002],
             [8.08500004, 3.19099998, -7.79099989],
             [-9.04500008, -26.46899986, 13.17999935],
             [2.51899981, 3.1400001, -1.73399997],
             [-10.64100075, -11.34899998, 0.257],
             [19.42700005, -8.21600056, -0.24399999]], dtype=np.float32)
        self.prec = 3

    def testParse_velocities(self):
        #read the velocities from the GRO_velocity file and compare the AtomGroup and individual Atom velocities
        # parsed with the reference values:
        u = MDAnalysis.Universe(GRO_velocity)
        all_atoms = u.select_atoms('all')
        #check for read-in and unit conversion for .gro file velocities for the entire AtomGroup:
        assert_almost_equal(all_atoms.velocities, self.reference_velocities, self.prec,
                            err_msg="problem reading .gro file velocities")
        #likewise for each individual atom (to be robust--in case someone alters the individual atom property code):
        assert_almost_equal(all_atoms[0].velocity, self.reference_velocities[0], self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[1].velocity, self.reference_velocities[1], self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[2].velocity, self.reference_velocities[2], self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[3].velocity, self.reference_velocities[3], self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[4].velocity, self.reference_velocities[4], self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[5].velocity, self.reference_velocities[5], self.prec,
                            err_msg="problem reading .gro file velocities")


class TestTRRForces(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(PDB_xvf, TRR_xvf)
        # extracted protein forces with g_traj into cobrotoxin_protein_forces.xvg.bz2
        # and manually averaged over 918 atoms and 3 time steps
        # native units: kJ/(mol*nm)
        self.reference_mean_protein_force_native = np.array(
            [3.4609879271822823, -0.63302345167392804, -1.0587882545813336], dtype=np.float32)
        # MDAnalysis units of kJ/(mol*A)
        self.reference_mean_protein_force = self.reference_mean_protein_force_native / 10
        self.prec = 6

    def tearDown(self):
        del self.universe

    @attr('slow')
    def testForces(self):
        protein = self.universe.select_atoms("protein")
        assert_equal(len(protein), 918)
        mean_F = np.mean([protein.forces.mean(axis=0) for ts in self.universe.trajectory], axis=0)
        assert_almost_equal(mean_F, self.reference_mean_protein_force, self.prec,
                            err_msg="mean force on protein over whole trajectory does not match")


class TestTRRForcesNativeUnits(TestTRRForces):
    def setUp(self):
        super(TestTRRForcesNativeUnits, self).setUp()
        # get universe without conversion
        self.universe = MDAnalysis.Universe(PDB_xvf, TRR_xvf, convert_units=False)
        # native Gromacs TRR units kJ/(mol*nm)
        self.reference_mean_protein_force = self.reference_mean_protein_force_native
