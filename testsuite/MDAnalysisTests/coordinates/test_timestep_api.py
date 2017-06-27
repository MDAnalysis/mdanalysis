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

"""Tests for the API definition of Timestep


_TestTimestepInterface tests the Readers are correctly using Timesteps
"""
from __future__ import absolute_import

from numpy.testing import assert_equal, dec
from MDAnalysisTests import module_not_found, parser_not_found

import MDAnalysis as mda
from MDAnalysisTests.datafiles import (PSF, XYZ_five, INPCRD, DCD, DLP_CONFIG,
                                       DLP_HISTORY, DMS, GMS_ASYMOPT, GRO, XTC,
                                       TRR, LAMMPSdata, LAMMPSdata2,
                                       LAMMPSdcd2, mol2_molecules, PDB_small,
                                       PDBQT_input, PQR, PRM, TRJ, PRMncdf,
                                       NCDF, TRZ_psf, TRZ)

from MDAnalysisTests.coordinates.base import BaseTimestepTest


# Can add in custom tests for a given Timestep here!
class TestBaseTimestep(BaseTimestepTest):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def test_other_timestep(self):
        # use a subclass to base.Timestep to check it works
        ts1 = mda.coordinates.base.Timestep(10)
        ts1.positions = self._get_pos()

        # can't do TRR Timestep here as it always has vels and forces
        # so isn't actually equal to a position only timestep
        for otherTS in [mda.coordinates.DCD.Timestep,
                        mda.coordinates.TRJ.Timestep,
                        mda.coordinates.DMS.Timestep,
                        mda.coordinates.GRO.Timestep,
                        mda.coordinates.TRZ.Timestep,
                        ]:
            ts2 = otherTS(10)
            ts2.positions = self._get_pos()
            yield (self._check_ts_equal, ts1, ts2,
                   "Failed on {0}".format(otherTS))


# TODO: Merge this into generic Reader tests
# These tests are all included in BaseReaderTest
# Once Readers use that TestClass, delete this one
class BaseTimestepInterfaceTest(object):
    """Test the Timesteps created by Readers

    This checks that Readers are creating correct Timestep objects,
    ensuring a consistent interface when using Timestep objects.

    Failures here are the Reader's fault

    See Issue #250 for discussion
    """
    def tearDown(self):
        del self.ts
        del self.u

    def test_frame(self):
        assert_equal(self.ts.frame, 0)

    def test_dt(self):
        assert_equal(self.u.trajectory.dt, self.ts.dt)


class TestCRD(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(XYZ_five, INPCRD)
        self.ts = u.trajectory.ts


class TestDCD(BaseTimestepInterfaceTest):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        u = self.u = mda.Universe(PSF, DCD)
        self.ts = u.trajectory.ts


class TestDLPConfig(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(DLP_CONFIG, format='CONFIG')
        self.ts = u.trajectory.ts


class TestDLPHistory(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(DLP_HISTORY, format='HISTORY')
        self.ts = u.trajectory.ts


class TestDMS(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(DMS)
        self.ts = u.trajectory.ts


class TestGMS(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(GMS_ASYMOPT, GMS_ASYMOPT,
                                  format='GMS', topology_format='GMS')
        self.ts = u.trajectory.ts


class TestGRO(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(GRO)
        self.ts = u.trajectory.ts


class TestINPCRD(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(XYZ_five, INPCRD)
        self.ts = u.trajectory.ts


class TestLAMMPS(BaseTimestepInterfaceTest):
    @dec.skipif(parser_not_found('LAMMPS'),
                'LAMMPS parser not available. Are you using python 3?')
    def setUp(self):
        u = self.u = mda.Universe(LAMMPSdata)
        self.ts = u.trajectory.ts


class TestLAMMPSDCD(BaseTimestepInterfaceTest):
    @dec.skipif(parser_not_found('LAMMPS'),
                'LAMMPS parser not available. Are you using python 3?')
    def setUp(self):
        u = self.u = mda.Universe(LAMMPSdata2, LAMMPSdcd2,
                                  format='LAMMPS', topology_format='DATA',
                                  timeunit='fs')
        self.ts = u.trajectory.ts


class TestMOL2(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(mol2_molecules)
        self.ts = u.trajectory.ts


class TestPDB(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(PDB_small)
        self.ts = u.trajectory.ts


class TestPDBQT(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(PDBQT_input)
        self.ts = u.trajectory.ts


class TestPQR(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(PQR)
        self.ts = u.trajectory.ts


class TestTRJ(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(PRM, TRJ)
        self.ts = u.trajectory.ts


class TestNCDF(BaseTimestepInterfaceTest):
    @dec.skipif(module_not_found("netCDF4"),
                "Test skipped because netCDF is not available.")
    def setUp(self):
        u = self.u = mda.Universe(PRMncdf, NCDF)
        self.ts = u.trajectory.ts


class TestTRR(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(GRO, TRR)
        self.ts = u.trajectory.ts


class TestTRZ(BaseTimestepInterfaceTest):
    @dec.skipif(parser_not_found('TRZ'),
                'TRZ parser not available. Are you using python 3?')
    def setUp(self):
        u = self.u = mda.Universe(TRZ_psf, TRZ)
        self.ts = u.trajectory.ts


class TestXTC(BaseTimestepInterfaceTest):
    def setUp(self):
        u = self.u = mda.Universe(GRO, XTC)
        self.ts = u.trajectory.ts
