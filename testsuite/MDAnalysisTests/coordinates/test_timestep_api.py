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

"""Tests for the API definition of Timestep


_TestTimestepInterface tests the Readers are correctly using Timesteps
"""
from __future__ import absolute_import

import numpy as np
from numpy.testing import assert_equal

import MDAnalysis as mda
from MDAnalysisTests.datafiles import (PSF, XYZ_five, INPCRD, DCD, DLP_CONFIG,
                                       DLP_HISTORY, DMS, GMS_ASYMOPT, GRO, XTC,
                                       TRR, LAMMPSdata, LAMMPSdata2,
                                       LAMMPSdcd2, mol2_molecules, PDB_small,
                                       PDBQT_input, PQR, PRM, TRJ, PRMncdf,
                                       NCDF, TRZ_psf, TRZ)

from MDAnalysisTests.coordinates.base import BaseTimestepTest, assert_timestep_equal
import pytest

# Can add in custom tests for a given Timestep here!


class TestBaseTimestep(BaseTimestepTest):
    @pytest.mark.parametrize('otherTS', [
        mda.coordinates.TRJ.Timestep,
        mda.coordinates.DMS.Timestep,
        mda.coordinates.GRO.Timestep,
        mda.coordinates.TRZ.Timestep,
    ])
    def test_other_timestep(self, otherTS):
        # use a subclass to base.Timestep to check it works
        ts1 = mda.coordinates.base.Timestep(10)
        ts1.positions = self._get_pos()
        ts2 = otherTS(10)
        ts2.positions = self._get_pos()
        assert_timestep_equal(ts1, ts2, "Failed on {0}".format(otherTS))


# TODO: Merge this into generic Reader tests
# These tests are all included in BaseReaderTest
# Once Readers use that TestClass, delete this one

class TestBaseTimestepInterface(object):
    """Test the Timesteps created by Readers

    This checks that Readers are creating correct Timestep objects,
    ensuring a consistent interface when using Timestep objects.

    Failures here are the Reader's fault

    See Issue #250 for discussion
    """
    @pytest.fixture(params=(
        (XYZ_five, INPCRD, None, None),
        (PSF, DCD, None, None),
        (DLP_CONFIG, None, 'CONFIG', None),
        (DLP_HISTORY, None, 'HISTORY', None),
        (DMS, None, None, None),
        (GRO, None, None, None),
        (XYZ_five, INPCRD, None, None),
        (LAMMPSdata, None, None, None),
        (mol2_molecules, None, None, None),
        (PDB_small, None, None, None),
        (PQR, None, None, None),
        (PDBQT_input, None, None, None),
        (PRM, TRJ, None, None),
        (GRO, XTC, None, None),
        (TRZ_psf, TRZ, None, None),
        (GRO, TRR, None, None),
        (GMS_ASYMOPT, GMS_ASYMOPT, 'GMS', 'GMS'),
        (LAMMPSdata2, LAMMPSdcd2, 'LAMMPS', 'DATA'),
        (PRMncdf, NCDF, None, None),
    ))
    def universe(self, request):
        topology, trajectory, trajectory_format, topology_format = request.param
        if trajectory_format is not None and topology_format is not None:
            return mda.Universe(topology, trajectory, format=trajectory_format,
                                topology_format=topology_format)

        if trajectory is not None:
            return mda.Universe(topology, trajectory)
        else:
            return mda.Universe(topology, format=trajectory_format)

    def test_frame(self, universe):
        assert_equal(universe.trajectory.ts.frame, 0)

    def test_dt(self, universe):
        assert_equal(universe.trajectory.dt, universe.trajectory.ts.dt)


@pytest.mark.parametrize('uni', [
    [(PSF, DCD), {}],  # base Timestep
    [(DLP_CONFIG,), {'format': 'CONFIG'}],  # DLPoly
    [(DMS,), {}],  # DMS
    [(GRO,), {}],  # GRO
    [(np.zeros((1, 30, 3)),), {}],  # memory
    [(PRM, TRJ), {}],  # TRJ
    [(TRZ_psf, TRZ), {}],  # TRZ
])
def test_atomgroup_dims_access(uni):
    uni_args, uni_kwargs = uni
    # check that AtomGroup.dimensions always returns a copy
    u = mda.Universe(*uni_args, **uni_kwargs)

    ag = u.atoms[:10]

    dims = ag.dimensions
    ts = u.trajectory.ts

    # dimensions from AtomGroup should be equal
    assert_equal(dims, ts.dimensions)
    # but not identical
    assert dims is not ts.dimensions
