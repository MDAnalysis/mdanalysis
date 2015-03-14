# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://mdanalysis.googlecode.com
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
from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small, GRO, TRR, \
    TRZ, TRZ_psf, \
    capping_input, capping_output, capping_ace, capping_nma, \
    merge_protein, merge_ligand, merge_water
import MDAnalysis.core.AtomGroup
from MDAnalysis.core.AtomGroup import Atom, AtomGroup
from MDAnalysis import NoDataError

import numpy
from numpy.testing import *
from numpy import array, float32, rad2deg
from nose.plugins.attrib import attr

import os
import tempfile
import itertools

try:
    from numpy.testing import assert_
except ImportError:
    # missing in numpy 1.2 but needed here:
    # copied code from numpy.testing 1.5
    def assert_(val, msg=''):
        """
        Assert that works in release mode.

        The Python built-in ``assert`` does not work when executing code in
        optimized mode (the ``-O`` flag) - no byte-code is generated for it.

        For documentation on usage, refer to the Python documentation.

        """
        if not val:
            raise AssertionError(msg)

from MDAnalysis import Universe, Merge
from MDAnalysis.analysis.align import alignto


def capping(ref, ace, nma, output):
    resids = ref.selectAtoms("all").resids()
    resid_min, resid_max = min(resids), max(resids)

    # There is probably some cache i need to update both the atom and residues
    for a in ace.atoms:
        a.resid += resid_min - max(ace.atoms.resids())
    for r in ace.residues:
        r.id += resid_min - max(ace.atoms.resids())
    for a in nma.atoms:
        a.resid = resid_max
    for r in nma.residues:
        r.id = resid_max

    # TODO pick the first residue in the protein (how should we cap the chains?)
    # TODO consider a case when the protein resid is 1 and all peptide has to be shifted by +1, put that in docs as a
    #  post-processing step
    alignto(ace, ref, select={
        "mobile": "resid {} and backbone".format(resid_min),
        "reference": "resid {} and backbone".format(resid_min)})
    alignto(nma, ref, select={
        "mobile": "resid {} and backbone and not (resname NMA or resname NME)".format(resid_max),
        "reference": "resid {} and (backbone or name OT2)".format(resid_max)})

    #  TODO remove the Hydrogen closest to ACE's oxygen
    u = Merge(ace.selectAtoms("resname ACE"),
              ref.selectAtoms(
                  "not (resid {} and name HT*) and not (resid {} and (name HT* or name OT1))".format(resid_min,
                                                                                                     resid_max)),
              nma.selectAtoms("resname NME or resname NMA"))
    u.trajectory.ts.dimensions = ref.trajectory.ts.dimensions
    u.atoms.write(output)
    return u


class TestCapping(TestCase):
    ext = "pdb"

    def setUp(self):
        suffix = '.' + self.ext
        fd, self.outfile = tempfile.mkstemp(suffix=suffix)
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass

    def test_capping_file(self):
        peptide = MDAnalysis.Universe(capping_input)
        ref = MDAnalysis.Universe(capping_output)
        ace = MDAnalysis.Universe(capping_ace)
        nma = MDAnalysis.Universe(capping_nma)

        u = capping(peptide, ace, nma, self.outfile)
        assert_equal(len(u.selectAtoms("not name H*")), len(ref.selectAtoms("not name H*")))

        u = MDAnalysis.Universe(self.outfile)

        ace = u.selectAtoms("resname ACE")
        nma = u.selectAtoms("resname NMA")

        # Check if the resids were assigned correctly
        assert_equal(ace.resids()[0], 1)
        assert_equal(nma.resids()[0], 15)

        assert_array_equal(peptide.trajectory.ts.dimensions, u.trajectory.ts.dimensions)

    def test_capping_inmemory(self):
        peptide = MDAnalysis.Universe(capping_input)
        ref = MDAnalysis.Universe(capping_output)
        ace = MDAnalysis.Universe(capping_ace)
        nma = MDAnalysis.Universe(capping_nma)

        u = capping(peptide, ace, nma, self.outfile)
        assert_equal(len(u.selectAtoms("not name H*")), len(ref.selectAtoms("not name H*")))

        ace = u.selectAtoms("resname ACE")
        nma = u.selectAtoms("resname NMA")

        # Check if the resids were assigned correctly
        assert_equal(ace.resids()[0], 1)
        assert_equal(nma.resids()[0], 15)

        assert_array_equal(peptide.trajectory.ts.dimensions, u.trajectory.ts.dimensions)


class TestMerge(TestCase):
    ext = "pdb"

    def setUp(self):
        u1 = MDAnalysis.Universe(merge_protein)
        u2 = MDAnalysis.Universe(merge_ligand)
        u3 = MDAnalysis.Universe(merge_water)
        self.universes = [u1, u2, u3]

        suffix = '.' + self.ext
        fd, self.outfile = tempfile.mkstemp(suffix=suffix)
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        for u in self.universes:
            del u

    def test_merge(self):
        u1, u2, u3 = self.universes
        ids_before = [a.number for u in [u1, u2, u3] for a in u.atoms]
        # Do the merge
        u0 = MDAnalysis.Merge(u1.atoms, u2.atoms, u3.atoms)
        # Check that the output Universe has the same number of atoms as the
        # starting AtomGroups
        assert_equal(len(u0.atoms), (len(u1.atoms) + len(u2.atoms) + len(u3.atoms)))

        # Make sure that all the atoms in the new universe are assigned to only
        # one, new Universe
        set0 = set([a.universe for a in u0.atoms])
        assert_equal(len(set0), 1)
        u = list(set0)[0]
        assert_equal(u, u0)

        # Make sure that the atom ids of the original universes are unchanged,
        # ie we didn't make the original Universes 'dirty'
        ids_after = [a.number for u in [u1, u2, u3] for a in u.atoms]
        assert_equal(len(ids_after), (len(u1.atoms) + len(u2.atoms) + len(u3.atoms)))
        assert_equal(ids_before, ids_after)

        # Test that we have a same number of atoms in a different way
        ids_new = [a.number for a in u0.atoms]
        assert_equal(len(ids_new), len(ids_before))

        u0.atoms.write(self.outfile)
        u = MDAnalysis.Universe(self.outfile)
        ids_new2 = [a.number for a in u.atoms]
        assert_equal(ids_new, ids_new2)

    def test_residue_references(self):
        u1, u2, u3 = self.universes
        m = Merge(u1.atoms, u2.atoms)
        assert_equal(m.atoms.residues[0].universe, m,
                     "wrong universe reference for residues after Merge()")

    def test_segment_references(self):
        u1, u2, u3 = self.universes
        m = Merge(u1.atoms, u2.atoms)
        assert_equal(m.atoms.segments[0].universe, m,
                     "wrong universe reference for segments after Merge()")

    def test_empty_ValueError(self):
        assert_raises(ValueError, Merge)

    def test_nonsense_TypeError(self):
        assert_raises(TypeError, Merge, ['1', 2])

    def test_emptyAG_ValueError(self):
        a = AtomGroup([])
        b = AtomGroup([])

        assert_raises(ValueError, Merge, a, b)
