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

import numpy as np
from numpy.testing import (
    assert_,
    assert_almost_equal,
    assert_array_equal,
    assert_array_almost_equal,
    assert_equal,
    assert_raises,
)
from nose.plugins.attrib import attr

from MDAnalysisTests import make_Universe
from MDAnalysisTests.datafiles import (
    COORDINATES_XYZ, COORDINATES_TRR,
    GRO, TRR,
    GRO_velocity, PDB_xvf, TRR_xvf
)

import MDAnalysis
from MDAnalysis import NoDataError


def assert_not_view(arr):
    assert_(arr.flags['OWNDATA'] == True)

def assert_correct_errormessage(func, var):
    errmsg = "Timestep does not contain {}".format(var)
    try:
        func[0](*func[1:])
    except NoDataError as e:
        assert_(errmsg in e.args[0])
    else:
        raise AssertionError


class TestAtomGroupTrajAccess(object):
    """
    For AtomGroup and Atom access:

    if present:
      - check return type
      - check dtype of array
      - check not a view of original (should always be copy!)
      - check the actual values returned
    if not present in trajectory:
      - check we get a NoDataError
      - check the error message of NDE

    For AtomGroup and Atom setting:

    if present:
      - check AtomGroup value is updated
      - check value in master Timestep object is updated
    if not present, check we get proper NoDataError on setting
    """
    @staticmethod
    def _check_atomgroup_positions_access(u, pos):
        ag = u.atoms[10:20]
        ag_pos = ag.positions

        assert_(isinstance(ag_pos, np.ndarray))
        assert_(ag_pos.dtype == np.float32)
        assert_not_view(ag_pos)
        assert_array_equal(ag_pos, u.trajectory.ts.positions[10:20])

    @staticmethod
    def _check_atomgroup_velocities_access(u, vel):
        ag = u.atoms[10:20]

        if vel:
            ag_vel = ag.velocities

            assert_(isinstance(ag_vel, np.ndarray))
            assert_(ag_vel.dtype == np.float32)
            assert_not_view(ag_vel)
            assert_array_equal(ag_vel, u.trajectory.ts.velocities[10:20])
        else:
            assert_raises(NoDataError, getattr, ag, 'velocities')
            assert_correct_errormessage((getattr, ag, 'velocities'), 'velocities')

    @staticmethod
    def _check_atomgroup_forces_access(u, force):
        ag = u.atoms[10:20]

        if force:
            ag_for = ag.forces

            assert_(isinstance(ag_for, np.ndarray))
            assert_(ag_for.dtype == np.float32)
            assert_not_view(ag_for)
            assert_array_equal(ag_for, u.trajectory.ts.forces[10:20])
        else:
            assert_raises(NoDataError, getattr, ag, 'forces')
            assert_correct_errormessage((getattr, ag, 'forces'), 'forces')

    @staticmethod
    def _check_atom_position_access(u, pos):
        at = u.atoms[55]

        at_pos = at.position

        assert_(isinstance(at_pos, np.ndarray))
        assert_(at_pos.dtype == np.float32)
        assert_not_view(at_pos)
        assert_array_equal(at_pos, u.trajectory.ts.positions[55])

    @staticmethod
    def _check_atom_velocity_access(u, vel):
        at = u.atoms[55]

        if vel:
            at_vel = at.velocity

            assert_(isinstance(at_vel, np.ndarray))
            assert_(at_vel.dtype == np.float32)
            assert_not_view(at_vel)
            assert_array_equal(at_vel, u.trajectory.ts.velocities[55])
        else:
            assert_raises(NoDataError, getattr, at, 'velocity')
            assert_correct_errormessage((getattr, at, 'velocity'), 'velocities')

    @staticmethod
    def _check_atom_force_access(u, force):
        at = u.atoms[55]

        if force:
            at_for = at.force

            assert_(isinstance(at_for, np.ndarray))
            assert_(at_for.dtype == np.float32)
            assert_not_view(at_for)
            assert_array_equal(at_for, u.trajectory.ts.forces[55])
        else:
            assert_raises(NoDataError, getattr, at, 'force')
            assert_correct_errormessage((getattr, at, 'force'), 'forces')

    @staticmethod
    def _check_atomgroup_positions_setting(u, pos):
        ag = u.atoms[[101, 107, 109]]

        new = np.array([[72.4, 64.5, 74.7],
                        [124.6, 15.6, -1.11],
                        [25.2, -66.6, 0]])

        ag.positions = new

        assert_array_almost_equal(ag.positions, new, decimal=5)
        assert_array_almost_equal(u.trajectory.ts.positions[[101, 107, 109]], new, decimal=5)

    @staticmethod
    def _check_atomgroup_velocities_setting(u, vel):
        ag = u.atoms[[101, 107, 109]]

        new = np.array([[72.4, 64.5, 74.7],
                        [124.6, 15.6, -1.11],
                        [25.2, -66.6, 0]]) + 0.1

        if vel:
            ag.velocities = new

            assert_array_almost_equal(ag.velocities, new, decimal=5)
            assert_array_almost_equal(u.trajectory.ts.velocities[[101, 107, 109]], new, decimal=5)
        else:
            assert_raises(NoDataError, setattr, ag, 'velocities', new)
            assert_correct_errormessage((setattr, ag, 'velocities', new), 'velocities')

    @staticmethod
    def _check_atomgroup_forces_setting(u, force):
        ag = u.atoms[[101, 107, 109]]

        new = np.array([[72.4, 64.5, 74.7],
                        [124.6, 15.6, -1.11],
                        [25.2, -66.6, 0]]) + 0.2

        if force:
            ag.forces = new

            assert_array_almost_equal(ag.forces, new, decimal=5)
            assert_array_almost_equal(u.trajectory.ts.forces[[101, 107, 109]], new, decimal=5)
        else:
            assert_raises(NoDataError, setattr, ag, 'forces', new)
            assert_correct_errormessage((setattr, ag, 'forces', new), 'forces')

    @staticmethod
    def _check_atom_position_setting(u, pos):
        at = u.atoms[94]

        new = np.array([58.3, -10.1, 0.001])

        at.position = new

        assert_array_almost_equal(at.position, new, decimal=5)
        assert_array_almost_equal(u.trajectory.ts.positions[94], new, decimal=5)

    @staticmethod
    def _check_atom_velocity_setting(u, vel):
        at = u.atoms[94]

        new = np.array([58.3, -10.1, 0.001]) + 0.1

        if vel:
            at.velocity = new

            assert_array_almost_equal(at.velocity, new, decimal=5)
            assert_array_almost_equal(u.trajectory.ts.velocities[94], new, decimal=5)
        else:
            assert_raises(NoDataError, setattr, at, 'velocity', new)
            assert_correct_errormessage((setattr, at, 'velocity', new), 'velocities')

    @staticmethod
    def _check_atom_force_setting(u, force):
        at = u.atoms[94]

        new = np.array([58.3, -10.1, 0.001]) + 0.2

        if force:
            at.force = new

            assert_array_almost_equal(at.force, new, decimal=5)
            assert_array_almost_equal(u.trajectory.ts.forces[94], new, decimal=5)
        else:
            assert_raises(NoDataError, setattr, at, 'force', new)
            assert_correct_errormessage((setattr, at, 'force', new), 'forces')

    def test_all(self):
        # all combinations of which trajectory attributes we have
        # positions is always present
        for pos, vel, force in (
                (True, False, False),
                (True, True, False),
                (True, False, True),
                (True, True, True),
        ):
            u = make_Universe(trajectory=pos, velocities=vel, forces=force)

            # AtomGroup access
            yield self._check_atomgroup_positions_access, u, pos
            yield self._check_atomgroup_velocities_access, u, vel
            yield self._check_atomgroup_forces_access, u, force
            # Atom access
            yield self._check_atom_position_access, u, pos
            yield self._check_atom_velocity_access, u, vel
            yield self._check_atom_force_access, u, force
            # AtomGroup setting
            yield self._check_atomgroup_positions_setting, u, pos
            yield self._check_atomgroup_velocities_setting, u, vel
            yield self._check_atomgroup_forces_setting, u, force
            # Atom setting
            yield self._check_atom_position_setting, u, pos
            yield self._check_atom_velocity_setting, u, vel
            yield self._check_atom_force_setting, u, force


class TestAtom_ForceVelocity(object):
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


class TestGROVelocities(object):
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


class TestTRRForces(object):
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


class TestAtomGroupVelocities(object):
    """Tests of velocity-related functions in AtomGroup"""
    def setUp(self):
        self.universe = MDAnalysis.Universe(GRO, TRR)
        self.ag = self.universe.select_atoms("bynum 12:42")

    def tearDown(self):
        del self.ag
        del self.universe

    @attr('slow')
    def test_get_velocities(self):
        v = self.ag.velocities
        assert_(np.any(np.abs(v) > 1e-6), "velocities should be non-zero")

    @attr('slow')
    def test_velocities(self):
        ag = self.universe.atoms[42:45]
        ref_v = np.array([
            [-3.61757946, -4.9867239, 2.46281552],
            [2.57792854, 3.25411797, -0.75065529],
            [13.91627216, 30.17778587, -12.16669178]])
        v = ag.velocities
        assert_almost_equal(v, ref_v, err_msg="velocities were not read correctly")

    @attr('slow')
    def test_set_velocities(self):
        ag = self.ag
        v = ag.velocities - 2.7271
        ag.velocities = v
        assert_almost_equal(ag.velocities, v,
                            err_msg="messages were not set to new value")


class TestAtomGroupForces(object):
    """Tests of velocity-related functions in AtomGroup"""
    def setUp(self):
        self.universe = MDAnalysis.Universe(COORDINATES_XYZ, COORDINATES_TRR)
        self.ag = self.universe.atoms[1:4]

    def tearDown(self):
        del self.universe

    @attr('slow')
    def test_get_forces(self):
        v = self.ag.forces
        assert_(np.any(np.abs(v) > 1e-6), "forces should be non-zero")

    @attr('slow')
    def test_forces(self):
        ag = self.universe.atoms[1:4]
        ref_v = np.arange(9).reshape(3, 3) * .01 + .03
        v = ag.forces
        assert_almost_equal(v, ref_v, err_msg="forces were not read correctly")

    @attr('slow')
    def test_set_forces(self):
        ag = self.ag
        v = ag.forces - 2.7271
        ag.forces = v
        assert_almost_equal(ag.forces, v,
                            err_msg="messages were not set to new value")


