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
from __future__ import absolute_import

import numpy as np
import pytest
from numpy.testing import assert_equal
import mock


from MDAnalysisTests.datafiles import PSF, GRO, XTC
from MDAnalysisTests import make_Universe

import MDAnalysis
import MDAnalysis as mda


class TestUpdatingSelection(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(GRO, XTC)

    @pytest.fixture()
    def ag(self, u):
        return u.select_atoms("prop x < 5 and prop y < 5 and prop z < 5")

    @pytest.fixture()
    def ag_updating(self, u):
        return u.select_atoms("prop x < 5 and prop y < 5 and prop z < 5",
                              updating=True)

    @pytest.fixture()
    def ag_updating_compounded(self, u, ag):
        return u.select_atoms("around 2 group sele",
                              sele=ag, updating=True)

    @pytest.fixture()
    def ag_updating_chained(self, u, ag_updating):
        return u.select_atoms("around 2 group sele",
                              sele=ag_updating, updating=True)

    @pytest.fixture()
    def ag_updating_chained2(self, ag_updating):
        return ag_updating.select_atoms("all", updating=True)

    def test_update(self, u, ag, ag_updating):
        assert_equal(ag_updating.indices, ag.indices)
        target_idxs = np.array([4469, 4470, 4472, 6289, 6290, 6291,
                                6292, 31313, 31314, 31315, 31316, 34661,
                                34663, 34664])
        next(u.trajectory)
        assert_equal(ag_updating._lastupdate, 0)
        assert ag_updating.is_uptodate is False
        assert_equal(ag_updating.indices, target_idxs)
        assert ag_updating.is_uptodate
        ag_updating.is_uptodate = False
        assert ag_updating._lastupdate is None

    def test_compounded_update(self, u, ag_updating_compounded):
        target_idxs0 = np.array([3650, 7406, 22703, 31426, 40357,
                                 40360, 41414])
        target_idxs1 = np.array([3650, 8146, 23469, 23472, 31426,
                                 31689, 31692, 34326, 41414])
        assert_equal(ag_updating_compounded.indices,
                     target_idxs0)
        next(u.trajectory)
        assert_equal(ag_updating_compounded.indices,
                     target_idxs1)

    def test_chained_update(self, u, ag_updating_chained,
                            ag_updating_compounded):
        target_idxs = np.array([4471, 7406, 11973, 11975, 34662, 44042])
        assert_equal(ag_updating_chained.indices,
                     ag_updating_compounded.indices)
        next(u.trajectory)
        assert_equal(ag_updating_chained.indices, target_idxs)

    def test_chained_update2(self, u, ag_updating, ag_updating_chained2):
        assert_equal(ag_updating_chained2.indices,
                     ag_updating.indices)
        next(u.trajectory)
        assert_equal(ag_updating_chained2.indices,
                     ag_updating.indices)

    def test_slice_is_static(self, u, ag, ag_updating):
        ag_static1 = ag_updating[:]
        ag_static2 = ag_updating.select_atoms("all")
        assert_equal(ag_static1.indices, ag.indices)
        assert_equal(ag_static2.indices, ag.indices)
        next(u.trajectory)
        assert_equal(ag_static1.indices, ag.indices)
        assert_equal(ag_static2.indices, ag.indices)

    def test_kwarg_check(self, u):
        with pytest.raises(TypeError):
            u.select_atoms("group updating", {"updating": True})

    def test_atoms(self, ag_updating):
        ag = ag_updating.atoms
        assert isinstance(ag, mda.core.groups.AtomGroup)
        assert ag is not ag_updating
        assert ag == ag_updating

    def test_copy(self, ag_updating):
        cgroup = ag_updating.copy()
        assert isinstance(cgroup, mda.core.groups.UpdatingAtomGroup)
        assert cgroup is not ag_updating
        assert cgroup == ag_updating


class TestUpdatingSelectionNotraj(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(PSF)

    @pytest.fixture()
    def ag(self, u):
        return u.select_atoms("name N*")

    @pytest.fixture()
    def ag_updating(self, u):
        return u.select_atoms("name N*", updating=True)

    def test_update(self, ag, ag_updating):
        assert ag_updating.is_uptodate
        assert_equal(ag_updating.indices, ag.indices)
        assert_equal(ag_updating._lastupdate, -1)
        ag_updating.is_uptodate = False
        assert ag_updating._lastupdate is None


class UAGReader(mda.coordinates.base.ReaderBase):
    """
    Positions in this reader are defined as:
    (atom number + frame number, 0, 0)

    Eg::
    Frame 1:
    (0, 0, 0),
    (1, 0, 0),
    etc

    Frame 2:
    (1, 0, 0),
    (2, 0, 0),
    etc

    Whilst quite possible not the best data for molecular simulation,
    it does make for easy to write tests.
    """

    def __init__(self, n_atoms):
        super(UAGReader, self).__init__('UAGReader')
        self._auxs = {}

        self.n_frames = 10
        self.n_atoms = n_atoms
        self.ts = self._Timestep(self.n_atoms)
        self._read_next_timestep()

    def _reopen(self):
        self.ts.frame = -1

    def _read_next_timestep(self):
        ts = self.ts
        ts.frame += 1
        if ts.frame >= self.n_frames:
            raise EOFError

        pos = np.zeros((self.n_atoms, 3))
        pos[:, 0] = np.arange(self.n_atoms) + ts.frame
        ts.positions = pos

        return ts

    def _read_frame(self, frame):
        self.ts.frame = frame - 1  # gets +1'd next
        return self._read_next_frame


class TestUAGCallCount(object):
    # make sure updates are only called when required!
    # 
    # these tests check that potentially expensive selection operations are only
    # done when necessary
    @pytest.fixture()
    def u(self):
        u = make_Universe(('names',))
        u.trajectory = UAGReader(125)
        return u

    @mock.patch.object(MDAnalysis.core.groups.UpdatingAtomGroup,
                       'update_selection',
                       autospec=True,
                       # required to make it get self when called
                       )
    def test_updated_when_creating(self, mock_update_selection, u):
        uag = u.select_atoms('name XYZ', updating=True)

        assert mock_update_selection.call_count == 1

    def test_updated_when_next(self, u):
        uag = u.select_atoms('name XYZ', updating=True)

        # Use mock.patch.object to start inspecting the uag update selection method
        # wraps= keyword makes it still function as normal, just we're spying on it now
        with mock.patch.object(uag, 'update_selection',
                               wraps=uag.update_selection) as mock_update:
            next(u.trajectory)
            assert mock_update.call_count == 0

            # Access many attributes..
            pos = uag.positions
            names = uag.names
            # But check we only got updated once
            assert mock_update.call_count == 1


class TestDynamicUAG(object):
    @pytest.fixture()
    def u(self):
        u = make_Universe(('names',))
        u.trajectory = UAGReader(125)
        return u

    def test_nested_uags(self, u):
        bg = u.atoms[[3, 4]]

        uag1 = u.select_atoms('around 1.5 group bg', bg=bg, updating=True)

        uag2 = u.select_atoms('around 1.5 group uag', uag=uag1, updating=True)

        for ts in u.trajectory:
            assert_equal(len(bg), 2)
            assert_equal(len(uag1), 2)  # around doesn't include bg, so 2
            assert_equal(len(uag2), 4)  # doesn't include uag1

    def test_driveby(self, u):
        uag = u.select_atoms('prop x < 5.5', updating=True)

        n_init = 6
        for i, ts in enumerate(u.trajectory):
            # should initially be 6 atoms with x < 5.5
            n_expected = max(n_init - i, 0)  # floors at 0

            assert_equal(len(uag), n_expected)


def test_representations():
    u = make_Universe()
    ag_updating = u.select_atoms("bynum 0", updating=True)
    rep = repr(ag_updating)
    assert "0 atoms," in rep
    assert "selection " in rep
    assert "bynum 0" in rep
    assert "entire Universe" in rep

    ag_updating = u.select_atoms("bynum 1", updating=True)
    rep = repr(ag_updating)
    assert "1 atom," in rep

    ag_updating = u.atoms[:-1].select_atoms("bynum 1", "bynum 2",
                                            updating=True)
    rep = repr(ag_updating)
    assert "2 atoms," in rep
    assert "selections 'bynum 1' + 'bynum 2'" in rep
    assert "another AtomGroup" in rep


def test_empty_UAG():
    u = make_Universe()

    # technically possible to make a UAG without any selections..
    uag = mda.core.groups.UpdatingAtomGroup(u.atoms, (), '')

    assert isinstance(uag, mda.core.groups.UpdatingAtomGroup)
