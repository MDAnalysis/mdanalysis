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
from numpy.testing import assert_equal

import pytest
import bz2
import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    PSF,
    PSF_nosegid,
    PSF_notop,
    PSF_NAMD,
    PSF_inscode,
    XYZ_psf,
    XYZ,
)

class PSFBase(ParserBase):
    parser = mda.topology.PSFParser.PSFParser
    expected_attrs = ['ids', 'names', 'types', 'masses',
                      'charges',
                      'resids', 'resnames',
                      'segids',
                      'bonds', 'angles', 'dihedrals', 'impropers']


class TestPSFParser(PSFBase):
    """
    Based on small PDB with AdK (:data:`PDB_small`).
    """
    ref_filename = PSF
    expected_n_atoms = 3341
    expected_n_residues = 214
    expected_n_segments = 1

    @pytest.fixture(params=['uncompressed', 'bz2'])
    def filename(self, request, tmpdir):
        if request.param == 'uncompressed':
            return self.ref_filename
        else:
            fn = str(tmpdir.join('file.psf.bz2'))
            with open(self.ref_filename, 'rb') as f:
                stuff = f.read()
            buf = bz2.compress(stuff)
            with open(fn, 'wb') as out:
                out.write(buf)
            return fn

    def test_bonds_total_counts(self, top):
        assert len(top.bonds.values) == 3365

    def test_bonds_atom_counts(self, filename):
        u = mda.Universe(filename)
        assert len(u.atoms[[0]].bonds) == 4
        assert len(u.atoms[[42]].bonds) == 1

    def test_bonds_identity(self, top):
        vals = top.bonds.values
        for b in ((0, 1), (0, 2), (0, 3), (0, 4)):
            assert (b in vals) or (b[::-1] in vals)

    def test_angles_total_counts(self, top):
        assert len(top.angles.values) == 6123

    def test_angles_atom_counts(self, filename):
        u = mda.Universe(filename)
        assert len(u.atoms[[0]].angles), 9
        assert len(u.atoms[[42]].angles), 2

    def test_angles_identity(self, top):
        vals = top.angles.values
        for b in ((1, 0, 2), (1, 0, 3), (1, 0, 4)):
            assert (b in vals) or (b[::-1] in vals)

    def test_dihedrals_total_counts(self, top):
        assert len(top.dihedrals.values) == 8921

    def test_dihedrals_atom_counts(self, filename):
        u = mda.Universe(filename)
        assert len(u.atoms[[0]].dihedrals) == 14

    def test_dihedrals_identity(self, top):
        vals = top.dihedrals.values
        for b in ((0, 4, 6, 7), (0, 4, 6, 8), (0, 4, 6, 9), (0, 4, 17, 18)):
            assert (b in vals) or (b[::-1] in vals)


class TestNAMDPSFParser(PSFBase):
    """Testfiles provided by JiyongPark77.

    NAMD/VMD XPLOR-style PSF file (using CGENFF residues/atoms).

    https://github.com/MDAnalysis/mdanalysis/issues/107
    """
    ref_filename = PSF_NAMD
    expected_n_atoms = 130
    expected_n_residues = 6
    expected_n_segments = 1


class TestPSFParser2(PSFBase):
    ref_filename = XYZ_psf
    expected_n_atoms = 1284
    expected_n_residues = 152
    expected_n_segments = 4

    def test_as_universe_resids(self):
        u = mda.Universe(XYZ_psf, XYZ)

        # each segment has [380, 381, 382, 383] as the first
        # 4 residues
        for seg in u.segments:
            assert_equal(seg.residues.resids[:4], [380, 381, 382, 383])

class TestPSFParserNoTop(PSFBase):
    ref_filename = PSF_notop
    expected_n_atoms = 3341
    expected_n_residues = 214
    expected_n_segments = 1

    def test_bonds_total_counts(self, top):
        assert len(top.bonds.values) == 0

    def test_angles_total_counts(self, top):
        assert len(top.angles.values) == 0

    def test_dihedrals_total_counts(self, top):
        assert len(top.dihedrals.values) == 0
    
    def test_impropers_total_counts(self, top):
        assert len(top.impropers.values) == 0

def test_psf_nosegid():
    """Issue #121"""
    u = mda.Universe(PSF_nosegid)
    assert isinstance(u, mda.Universe)
    assert u.atoms.n_atoms == 98
    assert_equal(u.segments.segids, ["SYSTEM"])

def test_psf_inscode():
    """Issue #2053 and #4189"""
    u = mda.Universe(PSF_inscode)
    assert_equal(u.residues.resids[:3], [1, 1, 1])