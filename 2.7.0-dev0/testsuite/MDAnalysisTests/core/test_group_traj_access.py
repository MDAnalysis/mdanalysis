# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import numpy as np
import pytest
from numpy.testing import (
    assert_almost_equal,
    assert_equal,
)

from MDAnalysisTests import make_Universe
from MDAnalysisTests.datafiles import (
    COORDINATES_XYZ, COORDINATES_TRR,
    GRO, TRR,
    GRO_velocity, PDB_xvf, TRR_xvf
)

import MDAnalysis
from MDAnalysis import NoDataError


def assert_not_view(arr):
    assert arr.flags['OWNDATA'] is True


def assert_correct_errormessage(func, var):
    errmsg = "Timestep has no {}".format(var)
    try:
        func[0](*func[1:])
    except NoDataError as e:
        assert errmsg in e.args[0]
    else:
        pytest.fail()


@pytest.mark.parametrize('pos,vel,force', (
    (True, False, False),
    (True, True, False),
    (True, False, True),
    (True, True, True),
), indirect=True)
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
    @pytest.fixture()
    def pos(self, request):
        return request.param

    @pytest.fixture()
    def force(self, request):
        return request.param

    @pytest.fixture()
    def vel(self, request):
        return request.param

    @pytest.fixture()
    def u(self, pos, vel, force):
        return make_Universe(trajectory=pos, velocities=vel, forces=force)

    def test_atomgroup_positions_access(self, u):
        ag = u.atoms[10:20]
        ag_pos = ag.positions

        assert isinstance(ag_pos, np.ndarray)
        assert ag_pos.dtype == np.float32
        assert_not_view(ag_pos)
        assert_equal(ag_pos, u.trajectory.ts.positions[10:20])

    def test_atomgroup_velocities_access(self, u, vel):
        ag = u.atoms[10:20]

        if vel:
            ag_vel = ag.velocities

            assert isinstance(ag_vel, np.ndarray)
            assert ag_vel.dtype == np.float32
            assert_not_view(ag_vel)
            assert_equal(ag_vel, u.trajectory.ts.velocities[10:20])
        else:
            with pytest.raises(NoDataError):
                getattr(ag, 'velocities')
            assert_correct_errormessage((getattr, ag, 'velocities'),
                                        'velocities')

    def test_atomgroup_forces_access(self, u, force):
        ag = u.atoms[10:20]

        if force:
            ag_for = ag.forces

            assert isinstance(ag_for, np.ndarray)
            assert ag_for.dtype == np.float32
            assert_not_view(ag_for)
            assert_equal(ag_for, u.trajectory.ts.forces[10:20])
        else:
            with pytest.raises(NoDataError):
                getattr(ag, 'forces')
            assert_correct_errormessage((getattr, ag, 'forces'), 'forces')

    def test_atom_position_access(self, u):
        at = u.atoms[55]

        at_pos = at.position

        assert isinstance(at_pos, np.ndarray)
        assert at_pos.dtype == np.float32
        assert_not_view(at_pos)
        assert_equal(at_pos, u.trajectory.ts.positions[55])

    def test_atom_velocity_access(self, u, vel):
        at = u.atoms[55]

        if vel:
            at_vel = at.velocity

            assert isinstance(at_vel, np.ndarray)
            assert at_vel.dtype == np.float32
            assert_not_view(at_vel)
            assert_equal(at_vel, u.trajectory.ts.velocities[55])
        else:
            with pytest.raises(NoDataError):
                getattr(at, 'velocity')
            assert_correct_errormessage(
                (getattr, at, 'velocity'), 'velocities')

    def test_atom_force_access(self, u, force):
        at = u.atoms[55]

        if force:
            at_for = at.force

            assert isinstance(at_for, np.ndarray)
            assert at_for.dtype == np.float32
            assert_not_view(at_for)
            assert_equal(at_for, u.trajectory.ts.forces[55])
        else:
            with pytest.raises(NoDataError):
                getattr(at, 'force')
            assert_correct_errormessage((getattr, at, 'force'), 'forces')

    def test_atomgroup_positions_setting(self, u):
        ag = u.atoms[[101, 107, 109]]

        new = np.array([[72.4, 64.5, 74.7],
                        [124.6, 15.6, -1.11],
                        [25.2, -66.6, 0]])

        ag.positions = new

        assert_almost_equal(ag.positions, new, decimal=5)
        assert_almost_equal(u.trajectory.ts.positions[[101, 107, 109]],
                            new, decimal=5)

    def test_atomgroup_velocities_setting(self, u, vel):
        ag = u.atoms[[101, 107, 109]]

        new = np.array([[72.4, 64.5, 74.7],
                        [124.6, 15.6, -1.11],
                        [25.2, -66.6, 0]]) + 0.1

        if vel:
            ag.velocities = new

            assert_almost_equal(ag.velocities, new, decimal=5)
            assert_almost_equal(
                u.trajectory.ts.velocities[[101, 107, 109]], new, decimal=5)
        else:
            with pytest.raises(NoDataError):
                setattr(ag, 'velocities', new)
            assert_correct_errormessage((setattr, ag, 'velocities', new),
                                        'velocities')

    def test_atomgroup_forces_setting(self, u, force):
        ag = u.atoms[[101, 107, 109]]

        new = np.array([[72.4, 64.5, 74.7],
                        [124.6, 15.6, -1.11],
                        [25.2, -66.6, 0]]) + 0.2

        if force:
            ag.forces = new

            assert_almost_equal(ag.forces, new, decimal=5)
            assert_almost_equal(u.trajectory.ts.forces[[101, 107, 109]],
                                new, decimal=5)
        else:
            with pytest.raises(NoDataError):
                setattr(ag, 'forces', new)
            assert_correct_errormessage((setattr, ag, 'forces', new), 'forces')

    def test_atom_position_setting(self, u):
        at = u.atoms[94]

        new = np.array([58.3, -10.1, 0.001])

        at.position = new

        assert_almost_equal(at.position, new, decimal=5)
        assert_almost_equal(u.trajectory.ts.positions[94], new, decimal=5)

    def test_atom_velocity_setting(self, u, vel):
        at = u.atoms[94]

        new = np.array([58.3, -10.1, 0.001]) + 0.1

        if vel:
            at.velocity = new

            assert_almost_equal(at.velocity, new, decimal=5)
            assert_almost_equal(u.trajectory.ts.velocities[94], new,
                                decimal=5)
        else:
            with pytest.raises(NoDataError):
                setattr(at, 'velocity', new)
            assert_correct_errormessage((setattr, at, 'velocity', new),
                                        'velocities')

    def test_atom_force_setting(self, u, force):
        at = u.atoms[94]

        new = np.array([58.3, -10.1, 0.001]) + 0.2

        if force:
            at.force = new

            assert_almost_equal(at.force, new, decimal=5)
            assert_almost_equal(u.trajectory.ts.forces[94], new,
                                decimal=5)
        else:
            with pytest.raises(NoDataError):
                setattr(at, 'force', new)
            assert_correct_errormessage((setattr, at, 'force', new), 'forces')


class TestAtom_ForceVelocity(object):
    @pytest.fixture()
    def u(self):
        return MDAnalysis.Universe(PDB_xvf, TRR_xvf)

    @pytest.fixture()
    def a(self, u):
        return u.atoms[0]

    def test_atom_force_get(self, u, a):
        assert_equal(a.force, u.atoms.forces[0])

    def test_atom_velocity_get(self, u, a):
        assert_equal(a.velocity, u.atoms.velocities[0])

    def test_atom_force_set(self, u, a):
        ref = np.arange(3)
        a.force = ref

        assert_equal(a.force, ref)
        assert_equal(u.atoms.forces[0], ref)

    def test_atom_velocity_set(self, u, a):
        ref = np.arange(3)
        a.velocity = ref

        assert_equal(a.velocity, ref)
        assert_equal(u.atoms.velocities[0], ref)

    def test_pos_iteration(self, u, a):
        ag = u.atoms[[0]]

        val = np.array([a.position for ts in u.trajectory])
        ref = np.array([ag.positions[0] for ts in u.trajectory])

        assert_equal(val, ref)

    def test_vel_iteration(self, u, a):
        ag = u.atoms[[0]]

        val = np.array([a.velocity for ts in u.trajectory])
        ref = np.array([ag.velocities[0] for ts in u.trajectory])

        assert_equal(val, ref)

    def test_for_iteration(self, u, a):
        ag = u.atoms[[0]]

        val = np.array([a.force for ts in u.trajectory])
        ref = np.array([ag.forces[0] for ts in u.trajectory])

        assert_equal(val, ref)


class TestGROVelocities(object):
    prec = 3

    @pytest.fixture()
    def reference_velocities(self):
        return np.array(
            [[-101.227, -0.57999998, 0.43400002],
             [8.08500004, 3.19099998, -7.79099989],
             [-9.04500008, -26.46899986, 13.17999935],
             [2.51899981, 3.1400001, -1.73399997],
             [-10.64100075, -11.34899998, 0.257],
             [19.42700005, -8.21600056, -0.24399999]], dtype=np.float32)

    def testParse_velocities(self, reference_velocities):
        # read the velocities from the GRO_velocity file and compare the AtomGroup and individual Atom velocities
        # parsed with the reference values:
        u = MDAnalysis.Universe(GRO_velocity)
        all_atoms = u.select_atoms('all')
        # check for read-in and unit conversion for .gro file velocities for the entire AtomGroup:
        assert_almost_equal(all_atoms.velocities, reference_velocities,
                            self.prec,
                            err_msg="problem reading .gro file velocities")
        # likewise for each individual atom (to be robust--in case someone alters the individual atom property code):
        assert_almost_equal(all_atoms[0].velocity, reference_velocities[0],
                            self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[1].velocity, reference_velocities[1],
                            self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[2].velocity, reference_velocities[2],
                            self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[3].velocity, reference_velocities[3],
                            self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[4].velocity, reference_velocities[4],
                            self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[5].velocity, reference_velocities[5],
                            self.prec,
                            err_msg="problem reading .gro file velocities")


class TestTRRForces(object):
    prec = 6

    @pytest.fixture()
    def universe(self):
        return MDAnalysis.Universe(PDB_xvf, TRR_xvf)

    @pytest.fixture()
    def reference_mean_protein_force(self):
        reference_mean_protein_force_native = np.array(
            [3.4609879271822823, -0.63302345167392804, -1.0587882545813336],
            dtype=np.float32)
        return reference_mean_protein_force_native / 10

    def testForces(self, universe, reference_mean_protein_force):
        protein = universe.select_atoms("protein")
        assert_equal(len(protein), 918)
        mean_F = np.mean(
            [protein.forces.mean(axis=0) for ts in universe.trajectory], axis=0)
        assert_almost_equal(mean_F, reference_mean_protein_force, self.prec,
                            err_msg="mean force on protein over whole trajectory does not match")


class TestTRRForcesNativeUnits(TestTRRForces):
    @pytest.fixture()
    def universe(self):
        return MDAnalysis.Universe(PDB_xvf, TRR_xvf, convert_units=False)

    @pytest.fixture()
    def reference_mean_protein_force(self):
        reference_mean_protein_force_native = np.array(
            [3.4609879271822823, -0.63302345167392804, -1.0587882545813336],
            dtype=np.float32)
        return reference_mean_protein_force_native


class TestAtomGroupVelocities(object):
    """Tests of velocity-related functions in AtomGroup"""

    @pytest.fixture()
    def universe(self):
        return MDAnalysis.Universe(GRO, TRR)

    @pytest.fixture()
    def ag(self, universe):
        return universe.select_atoms("bynum 12:42")

    def test_get_velocities(self, ag):
        v = ag.velocities
        assert np.any(np.abs(v) > 1e-6), "velocities should be non-zero"

    def test_velocities(self, universe):
        ag = universe.atoms[42:45]
        ref_v = np.array([
            [-3.61757946, -4.9867239, 2.46281552],
            [2.57792854, 3.25411797, -0.75065529],
            [13.91627216, 30.17778587, -12.16669178]])
        v = ag.velocities
        assert_almost_equal(v, ref_v,
                            err_msg="velocities were not read correctly")

    def test_set_velocities(self, ag):
        ag = ag
        v = ag.velocities - 2.7271
        ag.velocities = v
        assert_almost_equal(ag.velocities, v,
                            err_msg="messages were not set to new value")


class TestAtomGroupForces(object):
    """Tests of velocity-related functions in AtomGroup"""

    @pytest.fixture()
    def universe(self):
        return MDAnalysis.Universe(COORDINATES_XYZ, COORDINATES_TRR)

    @pytest.fixture()
    def ag(self, universe):
        return universe.atoms[1:4]

    def test_get_forces(self, ag):
        v = ag.forces
        assert np.any(np.abs(v) > 1e-6), "forces should be non-zero"

    def test_forces(self, universe):
        ag = universe.atoms[1:4]
        ref_v = np.arange(9).reshape(3, 3) * .01 + .03
        v = ag.forces
        assert_almost_equal(v, ref_v, err_msg="forces were not read correctly")

    def test_set_forces(self, ag):
        v = ag.forces - 2.7271
        ag.forces = v
        assert_almost_equal(ag.forces, v,
                            err_msg="messages were not set to new value")
