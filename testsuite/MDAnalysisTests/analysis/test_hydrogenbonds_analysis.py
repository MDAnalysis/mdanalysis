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
import logging
import pytest

from numpy.testing import (assert_allclose, assert_equal,
                           assert_array_almost_equal, assert_array_equal,
                           assert_almost_equal)

import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (
    HydrogenBondAnalysis)
from MDAnalysis.exceptions import NoDataError
from MDAnalysisTests.datafiles import waterPSF, waterDCD


class TestHydrogenBondAnalysisTIP3P(object):

    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(waterPSF, waterDCD)

    kwargs = {
        'donors_sel': 'name OH2',
        'hydrogens_sel': 'name H1 H2',
        'acceptors_sel': 'name OH2',
        'd_h_cutoff': 1.2,
        'd_a_cutoff': 3.0,
        'd_h_a_angle_cutoff': 120.0
    }

    @pytest.fixture(scope='class')
    def h(self, universe):
        h = HydrogenBondAnalysis(universe, **self.kwargs)
        h.run()
        return h

    def test_hbond_analysis(self, h):

        assert len(np.unique(h.hbonds[:, 0])) == 10
        assert len(h.hbonds) == 32

        reference = {
            'distance': {'mean': 2.7627309, 'std': 0.0905052},
            'angle': {'mean': 158.9038039, 'std': 12.0362826},
        }

        assert_allclose(np.mean(h.hbonds[:, 4]), reference['distance']['mean'])
        assert_allclose(np.std(h.hbonds[:, 4]), reference['distance']['std'])
        assert_allclose(np.mean(h.hbonds[:, 5]), reference['angle']['mean'])
        assert_allclose(np.std(h.hbonds[:, 5]), reference['angle']['std'])

    def test_count_by_time(self, h):

        ref_times = np.arange(0.02, 0.21, 0.02)
        ref_counts = np.array([3, 2, 4, 4, 4, 4, 3, 2, 3, 3])

        counts = h.count_by_time()
        assert_array_almost_equal(h.times, ref_times)
        assert_array_equal(counts, ref_counts)

    def test_count_by_type(self, h):

        # Only one type of hydrogen bond in this system
        ref_count = 32

        counts = h.count_by_type()
        assert int(counts[0, 2]) == ref_count

    def test_count_by_ids(self, h, universe):

        ref_counts = [1.0, 1.0, 0.5, 0.4, 0.2, 0.1]
        unique_hbonds = h.count_by_ids()

        most_common_hbond_ids = [12, 14, 9]
        assert_equal(unique_hbonds[0,:3], most_common_hbond_ids)

        # count_by_ids() returns raw counts
        # convert to fraction of time that bond was observed
        counts = unique_hbonds[:, 3] / len(h.times)

        assert_allclose(counts, ref_counts)


class TestHydrogenBondAnalysisIdeal(object):

    kwargs = {
        'donors_sel': 'name O',
        'hydrogens_sel': 'name H1 H2',
        'acceptors_sel': 'name O',
        'd_h_cutoff': 1.2,
        'd_a_cutoff': 3.0,
        'd_h_a_angle_cutoff': 120.0
    }

    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        # create two water molecules
        """
                       H4
                        \
            O1-H2 .... O2-H3
           /
          H1
        """
        n_residues = 2
        u = MDAnalysis.Universe.empty(
            n_atoms=n_residues*3,
            n_residues=n_residues,
            atom_resindex=np.repeat(range(n_residues), 3),
            residue_segindex=[0] * n_residues,
            trajectory=True,  # necessary for adding coordinates
            )

        u.add_TopologyAttr('name', ['O', 'H1', 'H2'] * n_residues)
        u.add_TopologyAttr('type', ['O', 'H', 'H'] * n_residues)
        u.add_TopologyAttr('resname', ['SOL'] * n_residues)
        u.add_TopologyAttr('resid', list(range(1, n_residues + 1)))
        u.add_TopologyAttr('id', list(range(1, (n_residues * 3) + 1)))

        # Atomic coordinates with a single hydrogen bond between O1-H2---O2
        pos1 = np.array([[0, 0, 0],             # O1
                        [-0.249, -0.968, 0],    # H1
                        [1, 0, 0],              # H2
                        [2.5, 0, 0],            # O2
                        [3., 0, 0],             # H3
                        [2.250, 0.968, 0]       # H4
                        ])

        # Atomic coordinates with no hydrogen bonds
        pos2 = np.array([[0, 0, 0],             # O1
                         [-0.249, -0.968, 0],   # H1
                         [1, 0, 0],             # H2
                         [4.5, 0, 0],           # O2
                         [5., 0, 0],            # H3
                         [4.250, 0.968, 0]      # H4
                         ])

        coordinates = np.empty((3,  # number of frames
                                u.atoms.n_atoms,
                                3))
        coordinates[0] = pos1
        coordinates[1] = pos2
        coordinates[2] = pos1
        u.load_new(coordinates, order='fac')

        return u

    @staticmethod
    @pytest.fixture(scope='class')
    def hydrogen_bonds(universe):
        h = HydrogenBondAnalysis(
            universe,
            **TestHydrogenBondAnalysisIdeal.kwargs
        )
        h.run()
        return h


    def test_no_bond_info_exception(self, universe):

        kwargs = {
            'donors_sel': None,
            'hydrogens_sel': None,
            'acceptors_sel': None,
            'd_h_cutoff': 1.2,
            'd_a_cutoff': 3.0,
            'd_h_a_angle_cutoff': 120.0
        }

        with pytest.raises(NoDataError, match="no bond information"):
            h = HydrogenBondAnalysis(universe, **kwargs)
            h._get_dh_pairs()

    def test_first_hbond(self, hydrogen_bonds):
        assert len(hydrogen_bonds.hbonds) == 2

        frame_no, donor_index, hydrogen_index, acceptor_index, da_dst, angle =\
            hydrogen_bonds.hbonds[0]
        assert_equal(donor_index, 0)
        assert_equal(hydrogen_index, 2)
        assert_equal(acceptor_index, 3)
        assert_almost_equal(da_dst, 2.5)
        assert_almost_equal(angle, 180)

    def test_count_by_time(self, hydrogen_bonds):
        ref_times = np.array([0, 1, 2]) # u.trajectory.dt is 1
        ref_counts = np.array([1, 0, 1])

        counts = hydrogen_bonds.count_by_time()
        assert_array_almost_equal(hydrogen_bonds.times, ref_times)
        assert_array_equal(counts, ref_counts)

    def test_hydrogen_bond_lifetime(self, hydrogen_bonds):
        tau_timeseries, timeseries = hydrogen_bonds.lifetime(tau_max=2)
        assert_array_equal(timeseries, [1, 0, 0])

    def test_hydrogen_bond_lifetime_intermittency(self, hydrogen_bonds):
        tau_timeseries, timeseries = hydrogen_bonds.lifetime(
            tau_max=2, intermittency=1
        )
        assert_array_equal(timeseries, 1)

    def test_no_attr_hbonds(self, universe):
        hbonds = HydrogenBondAnalysis(universe, **self.kwargs)
        # hydrogen bonds are not computed
        with pytest.raises(NoDataError, match=".hbonds attribute is None"):
            hbonds.lifetime(tau_max=2, intermittency=1)

    def test_logging_step_not_1(self, universe, caplog):
        hbonds = HydrogenBondAnalysis(universe, **self.kwargs)
        # using step 2
        hbonds.run(step=2)

        caplog.set_level(logging.WARNING)
        hbonds.lifetime(tau_max=2, intermittency=1)

        warning = \
            "Autocorrelation: Hydrogen bonds were computed with step > 1."
        assert any(warning in rec.getMessage() for rec in caplog.records)


class TestHydrogenBondAnalysisBetween(object):

    kwargs = {
        'donors_sel': 'name O P',
        'hydrogens_sel': 'name H1 H2 PH',
        'acceptors_sel': 'name O P',
        'd_h_cutoff': 1.2,
        'd_a_cutoff': 3.0,
        'd_h_a_angle_cutoff': 120.0
    }

    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        # create two water molecules and two "protein" molecules
        # P1-PH1 are the two atoms that comprise the toy protein PROT1
        """
                       H4
                        \
            O1-H2 .... O2-H3 ... P1-PH1 ... P2-PH2
           /
          H1
        """
        n_residues = 4
        n_sol_residues = 2
        n_prot_residues = 2
        u = MDAnalysis.Universe.empty(
            n_atoms=n_sol_residues*3 + n_prot_residues*2,
            n_residues=n_residues,
            atom_resindex=[0, 0, 0, 1, 1, 1, 2, 2, 3, 3],
            residue_segindex=[0, 0, 1, 1],
            trajectory=True,  # necessary for adding coordinates
            )

        u.add_TopologyAttr(
            'name',
            ['O', 'H1', 'H2'] * n_sol_residues + ['P', 'PH'] * n_prot_residues
        )
        u.add_TopologyAttr(
            'type',
            ['O', 'H', 'H'] * n_sol_residues + ['P', 'PH'] * n_prot_residues
        )
        u.add_TopologyAttr(
            'resname',
            ['SOL'] * n_sol_residues + ['PROT'] * n_prot_residues
        )
        u.add_TopologyAttr(
            'resid',
            list(range(1, n_residues + 1))
        )
        u.add_TopologyAttr(
            'id',
            list(range(1, (n_sol_residues * 3 + n_prot_residues * 2) + 1))
        )

        # Atomic coordinates with hydrogen bonds between:
        #     O1-H2---O2
        #     O2-H3---P1
        #     P1-PH1---P2
        pos = np.array([[0, 0, 0],             # O1
                        [-0.249, -0.968, 0],   # H1
                        [1, 0, 0],             # H2
                        [2.5, 0, 0],           # O2
                        [3., 0, 0],            # H3
                        [2.250, 0.968, 0],     # H4
                        [5.5, 0, 0],           # P1
                        [6.5, 0, 0],           # PH1
                        [8.5, 0, 0],           # P2
                        [9.5, 0, 0],           # PH2
                        ])

        coordinates = np.empty((1,  # number of frames
                                u.atoms.n_atoms,
                                3))
        coordinates[0] = pos
        u.load_new(coordinates, order='fac')

        return u

    def test_between_all(self, universe):
        # don't specify groups between which to find hydrogen bonds
        hbonds = HydrogenBondAnalysis(universe, between=None, **self.kwargs)
        hbonds.run()

        # indices of [donor, hydrogen, acceptor] for each hydrogen bond
        expected_hbond_indices = [
            [0, 2, 3],  # water-water
            [3, 4, 6],  # protein-water
            [6, 7, 8]   # protein-protein
        ]
        assert_array_equal(hbonds.hbonds[:, 1:4], expected_hbond_indices)

    def test_between_PW(self, universe):
        # Find only protein-water hydrogen bonds
        hbonds = HydrogenBondAnalysis(
            universe,
            between=["resname PROT", "resname SOL"],
            **self.kwargs
        )
        hbonds.run()

        # indices of [donor, hydrogen, acceptor] for each hydrogen bond
        expected_hbond_indices = [
            [3, 4, 6]  # protein-water
        ]
        assert_array_equal(hbonds.hbonds[:, 1:4], expected_hbond_indices)

    def test_between_PW_PP(self, universe):
        # Find protein-water and protein-protein hydrogen bonds (not
        # water-water)
        hbonds = HydrogenBondAnalysis(
            universe,
            between=[
                ["resname PROT", "resname SOL"],
                ["resname PROT", "resname PROT"]
            ],
            **self.kwargs
        )
        hbonds.run()

        # indices of [donor, hydrogen, acceptor] for each hydrogen bond
        expected_hbond_indices = [
            [3, 4, 6],  # protein-water
            [6, 7, 8]   # protein-protein
        ]
        assert_array_equal(hbonds.hbonds[:, 1:4], expected_hbond_indices)


class TestHydrogenBondAnalysisTIP3P_GuessAcceptors_GuessHydrogens_UseTopology_(TestHydrogenBondAnalysisTIP3P):
    """Uses the same distance and cutoff hydrogen bond criteria as :class:`TestHydrogenBondAnalysisTIP3P`, so the
    results are identical, but the hydrogens and acceptors are guessed whilst the donor-hydrogen pairs are determined
    via the topology.
    """
    kwargs = {
        'donors_sel': None,
        'hydrogens_sel': None,
        'acceptors_sel': None,
        'd_a_cutoff': 3.0,
        'd_h_a_angle_cutoff': 120.0
    }

    def test_no_hydrogens(self, universe):
        # If no hydrogens are identified at a given frame, check an
        # empty donor atom group is created
        test_kwargs = TestHydrogenBondAnalysisTIP3P.kwargs.copy()
        test_kwargs['donors_sel'] = None         # use topology to find pairs
        test_kwargs['hydrogens_sel'] = "name H"  # no atoms have name H

        h = HydrogenBondAnalysis(universe, **test_kwargs)
        h.run()

        assert h._hydrogens.n_atoms == 0
        assert h._donors.n_atoms == 0
        assert h.hbonds.size == 0

class TestHydrogenBondAnalysisTIP3P_GuessDonors_NoTopology(object):
    """Guess the donor atoms involved in hydrogen bonds using the partial charges of the atoms.
    """

    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(waterPSF, waterDCD)

    kwargs = {
        'donors_sel': None,
        'hydrogens_sel': None,
        'acceptors_sel': None,
        'd_h_cutoff': 1.2,
        'd_a_cutoff': 3.0,
        'd_h_a_angle_cutoff': 120.0
    }

    @pytest.fixture(scope='class')
    def h(self, universe):
        h = HydrogenBondAnalysis(universe, **self.kwargs)
        return h

    def test_guess_donors(self, h):

        ref_donors = "(resname TIP3 and name OH2)"
        donors = h.guess_donors(select='all', max_charge=-0.5)
        assert donors == ref_donors


class TestHydrogenBondAnalysisTIP3P_GuessHydrogens_NoTopology(object):
    """
    Guess the hydrogen atoms involved in hydrogen bonds using the mass and
    partial charge of the atoms.
    """

    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(waterPSF, waterDCD)

    kwargs = {
        'donors_sel': None,
        'hydrogens_sel': None,
        'acceptors_sel': None,
        'd_h_cutoff': 1.2,
        'd_a_cutoff': 3.0,
        'd_h_a_angle_cutoff': 120.0
    }

    @pytest.fixture(scope='class')
    def h(self, universe):
        h = HydrogenBondAnalysis(universe, **self.kwargs)
        return h

    def test_guess_hydrogens(self, h):

        ref_hydrogens = "(resname TIP3 and name H1) or (resname TIP3 and name H2)"
        hydrogens = h.guess_hydrogens(select='all')
        assert hydrogens == ref_hydrogens

    pytest.mark.parametrize(
        "min_mass, max_mass, min_charge",
        [
            (1.05, 1.10, 0.30),
            (0.90, 0.95, 0.30),
            (0.90, 1.10, 1.00)
        ]
    )
    def test_guess_hydrogens_empty_selection(self, h):

        hydrogens = h.guess_hydrogens(select='all', min_charge=1.0)
        assert hydrogens == ""

    def test_guess_hydrogens_min_max_mass(self, h):

        errmsg = r"min_mass is higher than \(or equal to\) max_mass"

        with pytest.raises(ValueError, match=errmsg):

            h.guess_hydrogens(select='all', min_mass=1.1, max_mass=0.9)

class TestHydrogenBondAnalysisTIP3PStartStep(object):
    """Uses the same distance and cutoff hydrogen bond criteria as :class:`TestHydrogenBondAnalysisTIP3P` but starting
    with the second frame and using every other frame in the analysis.
    """

    @staticmethod
    @pytest.fixture(scope='class')
    def universe():
        return MDAnalysis.Universe(waterPSF, waterDCD)

    kwargs = {
        'donors_sel': 'name OH2',
        'hydrogens_sel': 'name H1 H2',
        'acceptors_sel': 'name OH2',
        'd_h_cutoff': 1.2,
        'd_a_cutoff': 3.0,
        'd_h_a_angle_cutoff': 120.0
    }

    @pytest.fixture(scope='class')
    def h(self, universe):
        h = HydrogenBondAnalysis(universe, **self.kwargs)
        h.run(start=1, step=2)
        return h

    def test_hbond_analysis(self, h):

        assert len(np.unique(h.hbonds[:, 0])) == 5
        assert len(h.hbonds) == 15

        reference = {
            'distance': {'mean': 2.73942464, 'std': 0.05867924},
            'angle': {'mean': 157.07768079, 'std': 9.72636682},
        }

        assert_allclose(np.mean(h.hbonds[:, 4]), reference['distance']['mean'])
        assert_allclose(np.std(h.hbonds[:, 4]), reference['distance']['std'])
        assert_allclose(np.mean(h.hbonds[:, 5]), reference['angle']['mean'])
        assert_allclose(np.std(h.hbonds[:, 5]), reference['angle']['std'])

    def test_count_by_time(self, h):

        ref_times = np.array([0.04, 0.08, 0.12, 0.16, 0.20, ])
        ref_counts = np.array([2, 4, 4, 2, 3])

        counts = h.count_by_time()
        assert_array_almost_equal(h.times, ref_times)
        assert_array_equal(counts, ref_counts)

    def test_count_by_type(self, h):

        # Only one type of hydrogen bond in this system
        ref_count = 15

        counts = h.count_by_type()
        assert int(counts[0, 2]) == ref_count
