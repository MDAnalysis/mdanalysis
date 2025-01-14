# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2 or any higher version
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
import MDAnalysis as mda
import MDAnalysis.tests.datafiles as datafiles
import numpy as np
import pytest
from MDAnalysis import _GUESSERS, _TOPOLOGY_ATTRS
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import Atomnames, Atomtypes, Masses
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.guesser.base import GuesserBase, get_guesser
from numpy.testing import assert_allclose, assert_equal


class TestBaseGuesser:

    def test_get_guesser(self):
        class TestGuesser1(GuesserBase):
            context = "test1"

        class TestGuesser2(GuesserBase):
            context = "test2"

        assert get_guesser(TestGuesser1).context == "test1"
        assert get_guesser("test1").context == "test1"
        assert get_guesser(TestGuesser2()).context == "test2"

    def test_get_guesser_with_universe(self):
        class TestGuesser1(GuesserBase):
            context = "test1"

        u = mda.Universe.empty(n_atoms=5)
        guesser = get_guesser(TestGuesser1(), u, foo=1)

        assert len(guesser._universe.atoms) == 5
        assert "foo" in guesser._kwargs

    def test_guess_invalid_attribute(self):
        with pytest.raises(
            ValueError,
            match="default guesser can not guess "
            "the following attribute: foo",
        ):
            mda.Universe(datafiles.PDB, to_guess=["foo"])

    def test_guess_attribute_with_missing_parent_attr(self):
        names = Atomnames(np.array(["C", "HB", "HA", "O"], dtype=object))
        masses = Masses(
            np.array([np.nan, np.nan, np.nan, np.nan], dtype=np.float64)
        )
        top = Topology(4, 1, 1, attrs=[names, masses])
        u = mda.Universe(top, to_guess=["masses"])
        assert_allclose(
            u.atoms.masses,
            np.array([12.01100, 1.00800, 1.00800, 15.99900]),
            atol=0,
        )

    def test_force_guessing(self):
        names = Atomnames(np.array(["C", "H", "H", "O"], dtype=object))
        types = Atomtypes(np.array(["1", "2", "3", "4"], dtype=object))
        top = Topology(4, 1, 1, attrs=[names, types])
        u = mda.Universe(top, force_guess=["types"])
        assert_equal(u.atoms.types, ["C", "H", "H", "O"])

    def test_partial_guessing(self):
        types = Atomtypes(np.array(["C", "H", "H", "O"], dtype=object))
        masses = Masses(np.array([0, np.nan, np.nan, 0], dtype=np.float64))
        top = Topology(4, 1, 1, attrs=[types, masses])
        u = mda.Universe(top, to_guess=["masses"])
        assert_allclose(
            u.atoms.masses, np.array([0, 1.00800, 1.00800, 0]), atol=0
        )

    def test_force_guess_priority(self):
        "check that passing the attribute to force_guess have higher power"
        types = Atomtypes(np.array(["C", "H", "H", "O"], dtype=object))
        masses = Masses(np.array([0, np.nan, np.nan, 0], dtype=np.float64))
        top = Topology(4, 1, 1, attrs=[types, masses])
        u = mda.Universe(top, to_guess=["masses"], force_guess=["masses"])
        assert_allclose(
            u.atoms.masses,
            np.array([12.01100, 1.00800, 1.00800, 15.99900]),
            atol=0,
        )

    def test_partial_guess_attr_with_unknown_no_value_label(self):
        "trying to partially guess attribute tha doesn't have declared"
        "no_value_label should gives no effect"
        names = Atomnames(np.array(["C", "H", "H", "O"], dtype=object))
        types = Atomtypes(np.array(["", "", "", ""], dtype=object))
        top = Topology(4, 1, 1, attrs=[names, types])
        u = mda.Universe(top, to_guess=["types"])
        assert_equal(u.atoms.types, ["", "", "", ""])

    def test_guess_topology_objects_existing_read(self):
        u = mda.Universe(datafiles.CONECT)
        assert len(u.atoms.bonds) == 72
        assert list(u.bonds[0].indices) == [623, 630]

        # delete some bonds
        u.delete_bonds(u.atoms.bonds[:10])
        assert len(u.atoms.bonds) == 62
        # first bond has changed
        assert list(u.bonds[0].indices) == [1545, 1552]
        # count number of (1545, 1552) bonds
        ag = u.atoms[[1545, 1552]]
        bonds = ag.bonds.atomgroup_intersection(ag, strict=True)
        assert len(bonds) == 1
        assert not bonds[0].is_guessed

        all_indices = [tuple(x.indices) for x in u.bonds]
        assert (623, 630) not in all_indices

        # test guessing new bonds doesn't remove old ones
        u.guess_TopologyAttrs("default", to_guess=["bonds"])
        assert len(u.atoms.bonds) == 1922
        old_bonds = ag.bonds.atomgroup_intersection(ag, strict=True)
        assert len(old_bonds) == 1
        # test guessing new bonds doesn't duplicate old ones
        assert not old_bonds[0].is_guessed

        new_ag = u.atoms[[623, 630]]
        new_bonds = new_ag.bonds.atomgroup_intersection(new_ag, strict=True)
        assert len(new_bonds) == 1
        assert new_bonds[0].is_guessed

    def test_guess_topology_objects_existing_in_universe(self):
        u = mda.Universe(datafiles.CONECT, to_guess=["bonds"])
        assert len(u.atoms.bonds) == 1922
        assert list(u.bonds[0].indices) == [0, 1]

        # delete some bonds
        u.delete_bonds(u.atoms.bonds[:100])
        assert len(u.atoms.bonds) == 1822
        assert list(u.bonds[0].indices) == [94, 99]

        all_indices = [tuple(x.indices) for x in u.bonds]
        assert (0, 1) not in all_indices

        # guess old bonds back
        u.guess_TopologyAttrs("default", to_guess=["bonds"])
        assert len(u.atoms.bonds) == 1922
        # check TopologyGroup contains new (old) bonds
        assert list(u.bonds[0].indices) == [0, 1]

    def test_guess_topology_objects_force(self):
        u = mda.Universe(datafiles.CONECT, force_guess=["bonds"])
        assert len(u.atoms.bonds) == 1922

        with pytest.raises(NoDataError):
            u.atoms.angles

    def test_guess_topology_objects_out_of_order_init(self):
        u = mda.Universe(
            datafiles.PDB_small,
            to_guess=["dihedrals", "angles", "bonds"],
            guess_bonds=False,
        )
        assert len(u.atoms.angles) == 6123
        assert len(u.atoms.dihedrals) == 8921

    def test_guess_topology_objects_out_of_order_guess(self):
        u = mda.Universe(datafiles.PDB_small)
        with pytest.raises(NoDataError):
            u.atoms.angles

        u.guess_TopologyAttrs(
            "default", to_guess=["dihedrals", "angles", "bonds"]
        )
        assert len(u.atoms.angles) == 6123
        assert len(u.atoms.dihedrals) == 8921

    def test_force_guess_overwrites_existing_bonds(self):
        u = mda.Universe(datafiles.CONECT)
        assert len(u.atoms.bonds) == 72

        # This low radius should find no bonds
        vdw = dict.fromkeys(set(u.atoms.types), 0.1)
        u.guess_TopologyAttrs("default", to_guess=["bonds"], vdwradii=vdw)
        assert len(u.atoms.bonds) == 72

        # Now force guess bonds
        u.guess_TopologyAttrs("default", force_guess=["bonds"], vdwradii=vdw)
        assert len(u.atoms.bonds) == 0

    def test_guessing_angles_respects_bond_kwargs(self):
        u = mda.Universe(datafiles.PDB)
        assert not hasattr(u.atoms, "angles")

        # This low radius should find no angles
        vdw = dict.fromkeys(set(u.atoms.types), 0.01)

        u.guess_TopologyAttrs("default", to_guess=["angles"], vdwradii=vdw)
        assert len(u.atoms.angles) == 0

        # set higher radii for lots of angles!
        vdw = dict.fromkeys(set(u.atoms.types), 1)
        u.guess_TopologyAttrs("default", force_guess=["angles"], vdwradii=vdw)
        assert len(u.atoms.angles) == 89466

    def test_guessing_dihedrals_respects_bond_kwargs(self):
        u = mda.Universe(datafiles.CONECT)
        assert len(u.atoms.bonds) == 72

        u.guess_TopologyAttrs("default", to_guess=["dihedrals"])
        assert len(u.atoms.dihedrals) == 3548
        assert not hasattr(u.atoms, "angles")

    def test_guess_invalid_attribute(self):
        default_guesser = get_guesser("default")
        err = "not a recognized MDAnalysis topology attribute"
        with pytest.raises(KeyError, match=err):
            default_guesser.guess_attr("not_an_attribute")

    def test_guess_unsupported_attribute(self):
        default_guesser = get_guesser("default")
        err = "cannot guess this attribute"
        with pytest.raises(ValueError, match=err):
            default_guesser.guess_attr("tempfactors")

    def test_guess_singular(self):
        default_guesser = get_guesser("default")
        u = mda.Universe(datafiles.PDB, to_guess=[])
        assert not hasattr(u.atoms, "masses")

        default_guesser._universe = u
        masses = default_guesser.guess_attr("mass")


def test_Universe_guess_bonds_deprecated():
    with pytest.warns(
        DeprecationWarning, match="`guess_bonds` keyword is deprecated"
    ):
        u = mda.Universe(datafiles.PDB_full, guess_bonds=True)


@pytest.mark.parametrize(
    "universe_input",
    [datafiles.DCD, datafiles.XTC, np.random.rand(3, 3), datafiles.PDB],
)
def test_universe_creation_from_coordinates(universe_input):
    mda.Universe(universe_input)


def test_universe_creation_from_specific_array():
    a = np.array(
        [
            [0.0, 0.0, 150.0],
            [0.0, 0.0, 150.0],
            [200.0, 0.0, 150.0],
            [0.0, 0.0, 150.0],
            [100.0, 100.0, 150.0],
            [200.0, 100.0, 150.0],
            [0.0, 200.0, 150.0],
            [100.0, 200.0, 150.0],
            [200.0, 200.0, 150.0],
        ]
    )
    mda.Universe(a, n_atoms=9)
