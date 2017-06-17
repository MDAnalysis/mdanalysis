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
from __future__ import absolute_import
import MDAnalysis as mda
import os

from numpy.testing import (
    assert_,
    assert_almost_equal,
    assert_equal,
    assert_warns,
)

from MDAnalysisTests.coordinates.reference import RefAdKSmall
from MDAnalysisTests.coordinates.base import _SingleFrameReader
from MDAnalysisTests.datafiles import PQR
from MDAnalysisTests import tempdir, make_Universe


class TestPQRReader(_SingleFrameReader):
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
        assert_almost_equal(self.universe.atoms.H.charges,
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


class TestPQRWriter(RefAdKSmall):
    def setUp(self):
        self.universe = mda.Universe(PQR)
        self.prec = 3
        ext = ".pqr"
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/pqr-writer-test' + ext

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.universe
        del self.tmpdir

    def test_writer_noChainID(self):
        assert_equal(self.universe.segments.segids[0], 'SYSTEM')
        self.universe.atoms.write(self.outfile)
        u = mda.Universe(self.outfile)
        assert_equal(u.segments.segids[0], 'SYSTEM')
        assert_almost_equal(u.atoms.positions,
                            self.universe.atoms.positions, self.prec,
                            err_msg="Writing PQR file with PQRWriter does "
                            "not reproduce original coordinates")
        assert_almost_equal(u.atoms.charges, self.universe.atoms.charges,
                            self.prec, err_msg="Writing PQR file with "
                            "PQRWriter does not reproduce original charges")
        assert_almost_equal(u.atoms.radii, self.universe.atoms.radii,
                            self.prec, err_msg="Writing PQR file with "
                            "PQRWriter does not reproduce original radii")

    # 363 TODO:
    # Not sure if this should be a segid or chainID?
    # Topology system now allows for both of these
    def test_write_withChainID(self):
        self.universe.segments.segids = 'A'
        assert_equal(self.universe.segments.segids[0], 'A')  # sanity check
        self.universe.atoms.write(self.outfile)
        u = mda.Universe(self.outfile)
        assert_equal(u.segments.segids[0], 'A')
        assert_almost_equal(u.atoms.positions,
                            self.universe.atoms.positions, self.prec,
                            err_msg="Writing PQR file with PQRWriter does "
                            "not reproduce original coordinates")
        assert_almost_equal(u.atoms.charges, self.universe.atoms.charges,
                            self.prec, err_msg="Writing PQR file with "
                            "PQRWriter does not reproduce original charges")
        assert_almost_equal(u.atoms.radii, self.universe.atoms.radii,
                            self.prec, err_msg="Writing PQR file with "
                            "PQRWriter does not reproduce original radii")

    def test_timestep_not_modified_by_writer(self):
        ts = self.universe.trajectory.ts
        x = ts._pos.copy()
        self.universe.atoms.write(self.outfile)
        assert_equal(ts._pos,
                     x,
                     err_msg="Positions in Timestep were modified by writer.")

    def test_total_charge(self):
        self.universe.atoms.write(self.outfile)
        u = mda.Universe(self.outfile)
        assert_almost_equal(
            u.atoms.total_charge(), self.ref_charmm_totalcharge, 3,
            "Total charge (in CHARMM) does not match expected value.")

class TestPQRWriterMissingAttrs(object):
    # pqr requires names, resids, resnames, segids, radii, charges
    def setUp(self):
        self.reqd_attributes = ['names', 'resids', 'resnames', 'radii', 'charges']
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/pqr-writer-test.pqr'

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.tmpdir
        del self.outfile
        del self.reqd_attributes

    @staticmethod
    def assert_writing_warns(u, outfile):
        # write the test universe, and check warning is raised
        assert_warns(UserWarning, u.atoms.write, outfile)

    def test_no_names_writing(self):
        attrs = self.reqd_attributes
        attrs.remove('names')
        u = make_Universe(attrs, trajectory=True)

        self.assert_writing_warns(u, self.outfile)

        u2 = mda.Universe(self.outfile)

        assert_(all(u2.atoms.names == 'X'))

    def test_no_resnames_writing(self):
        attrs = self.reqd_attributes
        attrs.remove('resnames')
        u = make_Universe(attrs, trajectory=True)

        self.assert_writing_warns(u, self.outfile)

        u2 = mda.Universe(self.outfile)

        assert_(all(u2.residues.resnames == 'UNK'))

    def test_no_radii_writing(self):
        attrs = self.reqd_attributes
        attrs.remove('radii')
        u = make_Universe(attrs, trajectory=True)

        self.assert_writing_warns(u, self.outfile)

        u2 = mda.Universe(self.outfile)

        assert_(all(u2.atoms.radii == 1.0))

    def test_no_charges_writing(self):
        attrs = self.reqd_attributes
        attrs.remove('charges')
        u = make_Universe(attrs, trajectory=True)

        self.assert_writing_warns(u, self.outfile)

        u2 = mda.Universe(self.outfile)

        assert_(all(u2.atoms.charges == 0.0))
