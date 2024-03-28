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
import MDAnalysis as mda
import pytest

from numpy.testing import (
    assert_almost_equal,
    assert_equal,
)

from MDAnalysisTests.coordinates.reference import RefAdKSmall
from MDAnalysisTests.coordinates.base import _SingleFrameReader
from MDAnalysisTests.datafiles import PQR
from MDAnalysisTests import make_Universe


class TestPQRReader(_SingleFrameReader):
    __test__ = True
    def setUp(self):
        self.universe = mda.Universe(PQR)
        # 3 decimals in PDB spec
        # http://www.wwpdb.org/documentation/format32/sect9.html#ATOM
        self.prec = 3

    def test_total_charge(self):
        assert_almost_equal(
            self.universe.atoms.total_charge(), self.ref_charmm_totalcharge, 3,
            "Total charge (in CHARMM) does not match expected value.")

    def test_hydrogenCharges(self):
        assert_almost_equal(self.universe.select_atoms('name H').charges,
                            self.ref_charmm_Hcharges, 3,
                            "Charges for H atoms do not match.")

    # Note that the whole system gets the sysID 'SYSTEM' for the PQR file (when
    # read with a PSF it is 's4AKE')
    def test_ArgCACharges(self):
        ag = self.universe.select_atoms('resname ARG and name CA')
        assert_almost_equal(
            ag.charges, self.ref_charmm_ArgCAcharges,
            3, "Charges for CA atoms in Arg residues do not match.")

    def test_ProNCharges(self):
        ag = self.universe.select_atoms('resname PRO and name N')
        assert_almost_equal(
            ag.charges, self.ref_charmm_ProNcharges, 3,
            "Charges for N atoms in Pro residues do not match.")

    def test_dimensions(self):
        # Issue #3327 - dimensions should always be set to None
        assert self.universe.dimensions is None


class TestPQRWriter(RefAdKSmall):
    @staticmethod
    @pytest.fixture
    def universe():
        return mda.Universe(PQR)

    prec = 3

    @pytest.mark.parametrize('filename',
        ['test.pqr', 'test.pqr.bz2', 'test.pqr.gz'])
    def test_simple_writer_roundtrip(self, universe, filename, tmpdir):
        with tmpdir.as_cwd():
            universe.atoms.write(filename)
            u2 = mda.Universe(filename)
            assert_equal(universe.atoms.positions,
                         u2.atoms.positions)

    def test_writer_noChainID(self, universe, tmpdir):
        outfile = str(tmpdir.join('pqr-test.pqr'))

        assert_equal(universe.segments.segids[0], 'SYSTEM')
        universe.atoms.write(outfile)
        u = mda.Universe(outfile)
        assert_equal(u.segments.segids[0], 'SYSTEM')
        assert_almost_equal(u.atoms.positions,
                            universe.atoms.positions, self.prec,
                            err_msg="Writing PQR file with PQRWriter does "
                            "not reproduce original coordinates")
        assert_almost_equal(u.atoms.charges, universe.atoms.charges,
                            self.prec, err_msg="Writing PQR file with "
                            "PQRWriter does not reproduce original charges")
        assert_almost_equal(u.atoms.radii, universe.atoms.radii,
                            self.prec, err_msg="Writing PQR file with "
                            "PQRWriter does not reproduce original radii")

    # 363 TODO:
    # Not sure if this should be a segid or chainID?
    # Topology system now allows for both of these
    def test_write_withChainID(self, universe, tmpdir):
        outfile = str(tmpdir.join('pqr-test.pqr'))

        universe.segments.segids = 'A'
        assert_equal(universe.segments.segids[0], 'A')  # sanity check
        universe.atoms.write(outfile)
        u = mda.Universe(outfile)
        assert_equal(u.segments.segids[0], 'A')
        assert_almost_equal(u.atoms.positions,
                            universe.atoms.positions, self.prec,
                            err_msg="Writing PQR file with PQRWriter does "
                            "not reproduce original coordinates")
        assert_almost_equal(u.atoms.charges, universe.atoms.charges,
                            self.prec, err_msg="Writing PQR file with "
                            "PQRWriter does not reproduce original charges")
        assert_almost_equal(u.atoms.radii, universe.atoms.radii,
                            self.prec, err_msg="Writing PQR file with "
                            "PQRWriter does not reproduce original radii")

    def test_timestep_not_modified_by_writer(self, universe, tmpdir):
        outfile = str(tmpdir.join('pqr-test.pqr'))

        ts = universe.trajectory.ts
        x = ts.positions.copy()
        universe.atoms.write(outfile)
        assert_equal(ts.positions, x,
                     err_msg="Positions in Timestep were modified by writer.")

    def test_total_charge(self, universe, tmpdir):
        outfile = str(tmpdir.join('pqr-test.pqr'))
        universe.atoms.write(outfile)
        u = mda.Universe(outfile)
        assert_almost_equal(
            u.atoms.total_charge(), self.ref_charmm_totalcharge, 3,
            "Total charge (in CHARMM) does not match expected value.")

class TestPQRWriterMissingAttrs(object):
    # pqr requires names, resids, resnames, segids, radii, charges
    @staticmethod
    @pytest.fixture
    def reqd_attributes():
        return ['names', 'resids', 'resnames', 'radii', 'charges']

    @staticmethod
    @pytest.fixture
    def outfile(tmpdir):
        return str(tmpdir.join('pqr-writer-test.pqr'))

    def test_no_names_writing(self, reqd_attributes, outfile):
        attrs = reqd_attributes
        attrs.remove('names')
        u = make_Universe(attrs, trajectory=True)

        with pytest.warns(UserWarning):
            u.atoms.write(outfile)

        u2 = mda.Universe(outfile)

        assert all(u2.atoms.names == 'X')

    def test_no_resnames_writing(self, reqd_attributes, outfile):
        attrs = reqd_attributes
        attrs.remove('resnames')
        u = make_Universe(attrs, trajectory=True)

        with pytest.warns(UserWarning):
            u.atoms.write(outfile)

        u2 = mda.Universe(outfile)

        assert all(u2.residues.resnames == 'UNK')

    def test_no_radii_writing(self, reqd_attributes, outfile):
        attrs = reqd_attributes
        attrs.remove('radii')
        u = make_Universe(attrs, trajectory=True)

        with pytest.warns(UserWarning):
            u.atoms.write(outfile)

        u2 = mda.Universe(outfile)

        assert all(u2.atoms.radii == 1.0)

    def test_no_charges_writing(self, reqd_attributes, outfile):
        attrs = reqd_attributes
        attrs.remove('charges')
        u = make_Universe(attrs, trajectory=True)

        with pytest.warns(UserWarning):
            u.atoms.write(outfile)

        u2 = mda.Universe(outfile)

        assert all(u2.atoms.charges == 0.0)
