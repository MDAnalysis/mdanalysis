# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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
import numpy as np
import pytest
import warnings

from MDAnalysisTests.datafiles import waterPSF, waterDCD
from MDAnalysis.analysis.lineardensity import LinearDensity
from numpy.testing import assert_allclose
from MDAnalysis.core._get_readers import get_reader_for
from MDAnalysisTests.util import no_deprecated_call

from MDAnalysis.units import constants



def test_invalid_grouping():
    """Invalid groupings raise AttributeError"""
    universe = mda.Universe(waterPSF, waterDCD)
    sel_string = "all"
    selection = universe.select_atoms(sel_string)
    with pytest.raises(AttributeError):
        # centroid is attribute of AtomGroup, but not valid here
        ld = LinearDensity(selection, grouping="centroid", binsize=5)
        ld.run()


### Creating the Test Universes ###

def make_Universe(coords, charges, masses, n_atoms, n_frames, atomsPerRes, resPerSegs):
    """Generate a reference universe of 100 atoms with uniformly drawn random positions."""
    n_residues = n_atoms // atomsPerRes # Arbitrarily 5 atoms per residue
    n_segments = n_residues // resPerSegs # Arbitrarily 4 residues per segment

    # Indexing atoms into residues & residues into segments
    atom_resindex = np.repeat(np.arange(n_residues), atomsPerRes)
    residue_segindex = np.repeat(np.arange(n_segments), resPerSegs)

    # Creating the universe
    u = mda.Universe.empty(n_atoms=n_atoms,
                           n_residues=n_residues,
                           n_segments=n_segments,
                           atom_resindex=atom_resindex,
                           residue_segindex=residue_segindex)

    # Assigning the Charges & Masses
    u.add_TopologyAttr('charges', values=charges)
    u.add_TopologyAttr('masses', values=masses)

    u.trajectory = get_reader_for(coords)(coords,
                                          order='fac',
                                          n_atoms=n_atoms)

    for ts in u.trajectory:
        ts.dimensions = np.array([1, 1, 1, 90, 90, 90])

    return u


### Defining the Systems ###

# The masses are all set to 1
# The charges are dependent on each system
# The positions are randomly taken from a uniform distribution between 0 and 1
# NOTE: The positions at each time frame are independent on the previous time step

def neutral_Particles(n_atoms, n_frames, atomsPerRes, resPerSegs):
    charges = np.zeros(n_atoms)
    masses = np.ones(n_atoms)

    shape = (n_frames, n_atoms, 3)
    coords = np.random.random(shape)

    return make_Universe(coords, charges, masses, n_atoms, n_frames, atomsPerRes, resPerSegs)

def charged_Particles(n_atoms, n_frames, atomsPerRes, resPerSegs):
    charges = np.random.random(n_atoms) * 2 - np.ones(n_atoms) # Between -1 and 1
    masses = np.ones(n_atoms)

    shape = (n_frames, n_atoms, 3)
    coords = np.random.random(shape)

    return make_Universe(coords, charges, masses, n_atoms, n_frames, atomsPerRes, resPerSegs)

def charged_Dimers(n_dimers=100, n_frames=1, dimersPerRes=5, resPerSegs=4, dimerLength = 0.05):
    n_atoms = 2 * n_dimers
    
    charges = np.random.random(n_atoms) * 2 - np.ones(n_atoms) # Between -1 and 1
    masses = np.ones(n_atoms)

    shape = (n_frames, n_dimers, 3)
    coords = np.random.random(shape) * 0.9 + np.ones(shape) * 0.05
    # Puts in the same coordinate twice per dimer
    coords = np.repeat(coords, 2, axis = 1)

    # Shifts one of the atoms of each dimer by their bondLength in a random direction (defined to be in the box)
    for time in coords:
        for coord in time[::2,:]:
            phi = np.random.random() * 2 * np.pi
            theta = np.random.random() * np.pi
            x = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]) * dimerLength
            coord += np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]) * dimerLength


    return make_Universe(coords, charges, masses, n_atoms, n_frames, dimersPerRes * 2, resPerSegs)



### Calculating the Expected Values ###

def calc_Prop(u, prop = 'masses'): # Property can be 'masses' or 'charges'
    expected_atoms = getattr(u.atoms, prop)
    expected_residues = np.array([sum(getattr(res.atoms, prop)) for res in u.residues])
    expected_segments = np.array([sum(getattr(seg.atoms, prop)) for seg in u.segments])

    return expected_atoms, expected_residues, expected_segments

def find_Centres(groups): # Always based on Centre of Mass (NOT Centre of Charge)
    centres = []
    for group in groups:
        total_Mass = sum(group.atoms.masses)
        if total_Mass != 0:
            centres.append(np.sum(group.atoms.positions.transpose() * abs(group.atoms.masses), axis = 1) / total_Mass)
        # Just in case there are particles with 0 mass (somehow)
        elif total_Mass == 0:
            centres.append(np.sum(group.atoms.positions.transpose() * abs(group.atoms.masses), axis = 1) / len(group.atoms))

    return np.array(centres)

def calc_Densities(u, prop = 'masses', spliceLen = 0.25): # Property can be 'masses' or 'charges'

    propShort = 'mass'
    if prop == 'charges':
        propShort = 'charge'

    ### Atoms
    expected_atoms = np.zeros((3, int(u.dimensions[0] // spliceLen))).astype(float) # Works for cubic Universe
    for atom in u.atoms:
        for i in range(3):
            expected_atoms[i][int(atom.position[i] // spliceLen)] += getattr(atom, propShort)



    _,residue_Totals,segment_Totals = calc_Prop(u, prop) # Total of Charge OR Mass
    ### Residues
    expected_residues = np.zeros((3, int(u.dimensions[0] // spliceLen))).astype(float)
    residue_Centres = find_Centres(u.residues)


    for i in range(len(residue_Centres)):
        for j in range(3):
            expected_residues[j][int(residue_Centres[i][j] // spliceLen)] += residue_Totals[i]

    ### Segments
    expected_segments = np.zeros((3, int(u.dimensions[0] // spliceLen))).astype(float)
    segment_Centres = find_Centres(u.segments)


    for i in range(len(segment_Centres)):
        for j in range(3):
            expected_segments[j][int(segment_Centres[i][j] // spliceLen)] += segment_Totals[i]
    
    
    # Scaling based on splice volumes & converting units
    for i in range(3):
        expected_atoms[i,:] /= spliceLen * u.dimensions[(i + 1) % 3] * u.dimensions[(i + 2) % 3]
        expected_residues[i,:] /= spliceLen * u.dimensions[(i + 1) % 3] * u.dimensions[(i + 2) % 3]
        expected_segments[i,:] /= spliceLen * u.dimensions[(i + 1) % 3] * u.dimensions[(i + 2) % 3]
    expected_atoms /= constants['N_Avogadro'] * 1e-24 # To be consistent with lineardensity.py
    expected_residues /= constants['N_Avogadro'] * 1e-24 # To be consistent with lineardensity.py
    expected_segments /= constants['N_Avogadro'] * 1e-24 # To be consistent with lineardensity.py

    return expected_atoms, expected_residues, expected_segments

### Creating the Expected Values ###

spliceLen = 0.25

universes = []
groupings = ['atoms', 'residues', 'segments']
# Will be [[Atoms, Residues, Segments (of Universe 1)], [... of Universe 2], ...]
expected_masses = [] 
expected_charges = []
expected_xmass = []
expected_xcharge = []

test_Universes = [neutral_Particles, charged_Particles, charged_Dimers]
test_Params = [[100, 1, 5, 4], [100, 1, 5, 4], [100, 1, 5, 4, 0.05]]

for i in range(len(test_Universes)):
    
    universe = test_Universes[i](*test_Params[i])
    universes.append(universe)

    expected_masses.append(calc_Prop(universe, 'masses'))
    expected_charges.append(calc_Prop(universe, 'charges'))
    expected_xmass.append(calc_Densities(universe, 'masses', spliceLen))
    expected_xcharge.append(calc_Densities(universe, 'charges', spliceLen))

### Parametrizing for Pytest ###

# NOTE: There is an extra [0] after the expected_xmass and expected_xcharge as the original data has all 3 co-ords, but only comparing to x here
pytest_Params = [(universes[i], groupings[j], expected_masses[i][j], expected_charges[i][j], expected_xmass[i][j][0], expected_xcharge[i][j][0]) for j in range(len(groupings)) for i in range(len(universes))]
@pytest.mark.parametrize(
    "universe, grouping, expected_masses, expected_charges, expected_xmass, expected_xcharge",
    pytest_Params,
)

def test_lineardensity(
    universe,
    grouping,
    expected_masses,
    expected_charges,
    expected_xmass,
    expected_xcharge,
):
    sel_string = "all"
    selection = universe.select_atoms(sel_string)
    ld = LinearDensity(selection, grouping, binsize=0.25).run()
    assert_allclose(ld.masses, expected_masses)
    assert_allclose(ld.charges, expected_charges)
    # rtol changed here due to floating point imprecision
    assert_allclose(ld.results.x.mass_density, expected_xmass, rtol=1e-06)
    assert_allclose(ld.results.x.charge_density, expected_xcharge)

@pytest.fixture(scope="module")
def testing_Universe():
    """Generate a universe for testing whether LinearDensity works with
    updating atom groups. Also used for parallel analysis test."""
    n_atoms = 3
    u = mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=n_atoms,
        n_segments=n_atoms,
        atom_resindex=np.arange(n_atoms),
        residue_segindex=np.arange(n_atoms),
    )

    for attr in ["charges", "masses"]:
        u.add_TopologyAttr(attr, values=np.ones(n_atoms))

    coords = np.array(
        [
            [[1.0, 1.0, 1.0], [1.0, 2.0, 1.0], [2.0, 1.0, 1.0]],
            [[1.0, 1.0, 2.0], [1.0, 2.0, 1.0], [2.0, 1.0, 1.0]],
            [[1.0, 1.0, 3.0], [1.0, 2.0, 1.0], [2.0, 1.0, 1.0]],
            [[1.0, 1.0, 4.0], [1.0, 2.0, 1.0], [2.0, 1.0, 1.0]],
            [[1.0, 1.0, 5.0], [1.0, 2.0, 1.0], [2.0, 1.0, 1.0]],
        ]
    )

    u.trajectory = get_reader_for(coords)(coords, order="fac", n_atoms=n_atoms)

    for ts in u.trajectory:
        ts.dimensions = np.array([2, 2, 6, 90, 90, 90])

    return u


def test_updating_atomgroup(testing_Universe):
    expected_z_pos = np.array([0.0, 0.91329641, 0.08302695, 0.0, 0.0, 0.0])
    u = testing_Universe
    selection = u.select_atoms("prop z < 3", updating=True)
    ld = LinearDensity(selection, binsize=1).run()
    assert_allclose(ld.results.z.mass_density, expected_z_pos)
    # Test whether histogram bins are saved correctly.
    expected_bin_edges = np.arange(0, 7)
    assert_allclose(ld.results.x.hist_bin_edges, expected_bin_edges)


testdict = {
    "pos": "mass_density",
    "pos_std": "mass_density_stddev",
    "char": "charge_density",
    "char_std": "charge_density_stddev",
}


# TODO: Remove in 3.0.0
def test_old_name_deprecations():
    universe = mda.Universe(waterPSF, waterDCD)
    sel_string = "all"
    selection = universe.select_atoms(sel_string)
    ld = LinearDensity(selection, binsize=5).run()
    with pytest.warns(DeprecationWarning):
        assert_allclose(ld.results.x.pos, ld.results.x.mass_density)
        assert_allclose(ld.results.x.pos_std, ld.results.x.mass_density_stddev)
        assert_allclose(ld.results.x.char, ld.results.x.charge_density)
        assert_allclose(
            ld.results.x.char_std, ld.results.x.charge_density_stddev
        )
        for key in testdict.keys():
            assert_allclose(
                ld.results["x"][key], ld.results["x"][testdict[key]]
            )

    # Check that no DeprecationWarning is raised with new attributes
    with no_deprecated_call():
        ld.results.x.mass_density
        ld.results.x.mass_density_stddev
        ld.results.x.charge_density
        ld.results.x.charge_density_stddev


# TODO: deprecated, remove in 3.0.0
def test_parallel_analysis(testing_Universe):
    """tests _add_other_result() method. Runs LinearDensity for all atoms of
    a universe and for two subsets, then adds the results of the two subsets
    and checks the results are the same."""
    u = testing_Universe
    selection1 = u.select_atoms("prop x < 1.1")
    selection2 = u.select_atoms("prop x >= 1.1")
    selection_whole = u.select_atoms("all")
    ld1 = LinearDensity(selection1, binsize=1).run()
    ld2 = LinearDensity(selection2, binsize=1).run()
    ld_whole = LinearDensity(selection_whole, binsize=1).run()
    with pytest.warns(
        DeprecationWarning, match="`_add_other_results` is deprecated!"
    ):
        ld1._add_other_results(ld2)
    assert_allclose(
        ld1.results.z.mass_density, ld_whole.results.z.mass_density
    )
    assert_allclose(
        ld1.results.x.mass_density, ld_whole.results.x.mass_density
    )

