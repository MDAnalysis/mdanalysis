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
import multiprocessing

import numpy as np
import os
import shutil
import pytest
import pickle
from numpy.testing import assert_equal

import MDAnalysis as mda
import MDAnalysis.coordinates
from MDAnalysis.coordinates.core import get_reader_for
from MDAnalysis.coordinates.GSD import HAS_GSD
from MDAnalysis.analysis.rms import RMSD

from MDAnalysisTests.datafiles import (
    CRD,
    PSF, DCD,
    DMS,
    DLP_CONFIG,
    DLP_HISTORY,
    FHIAIMS,
    INPCRD,
    GMS_ASYMOPT,
    GMS_SYMOPT,
    GRO,
    GSD,
    GSD_long,
    LAMMPSdata_mini,
    LAMMPSDUMP,
    mol2_molecules,
    MMTF,
    NCDF,
    PDB, PDB_small, PDB_multiframe,
    PDBQT_input,
    PQR,
    TRC_PDB_VAC, TRC_TRAJ1_VAC, TRC_TRAJ2_VAC,
    TRR,
    TRJ,
    TRZ,
    TXYZ,
    XTC,
    XPDB_small,
    XYZ_mini, XYZ, XYZ_bz2,
)


@pytest.fixture(params=[
    (PSF, DCD),
    (GRO, XTC),
    (PDB_multiframe,),
    (XYZ,),
    (XYZ_bz2,),  # .bz2
    (GMS_SYMOPT,),  # .gms
    (GMS_ASYMOPT,),  # .gz
    pytest.param(
        (GSD_long,),
        marks=pytest.mark.skipif(not HAS_GSD, reason='gsd not installed')
    ),
    (NCDF,),
    (np.arange(150).reshape(5, 10, 3).astype(np.float64),),
    (GRO, [GRO, GRO]),
    (PDB, [PDB, PDB]),
    (GRO, [XTC, XTC]),
    (TRC_PDB_VAC, TRC_TRAJ1_VAC),
    (TRC_PDB_VAC, [TRC_TRAJ1_VAC, TRC_TRAJ2_VAC]),
])
def u(request):
    if len(request.param) == 1:
        f = request.param[0]
        return mda.Universe(f)
    else:
        top, trj = request.param
        return mda.Universe(top, trj)

@pytest.fixture(scope="function")
def temp_xtc(tmp_path):
    fresh_xtc = tmp_path / "testing.xtc"
    shutil.copy(XTC, fresh_xtc)
    # In principle there is no possibility that the offset for
    # fresh_xtc exists as it is a fresh copy for every test. However,
    # the code is left for documentation.
    try:
        os.remove(MDAnalysis.coordinates.XDR.offsets_filename(fresh_xtc))
    except OSError:
        pass
    return fresh_xtc


# Define target functions here
# inside test functions doesn't work
def cog(u, ag, frame_id):
    u.trajectory[frame_id]

    return ag.center_of_geometry()


def test_multiprocess_COG(u):
    ag = u.atoms[2:5]

    ref = np.array([cog(u, ag, i)
                    for i in range(2)])

    p = multiprocessing.Pool(2)
    res = np.array([p.apply(cog, args=(u, ag, i))
                    for i in range(2)])
    p.close()
    assert_equal(ref, res)


def getnames(u, ix):
    # Check topology stuff works
    return u.atoms[ix].name


def test_universe_unpickle_in_new_process():
    u = mda.Universe(GRO, XTC)
    ref = [getnames(u, i)
           for i in range(3)]

    p = multiprocessing.Pool(2)
    res = [p.apply(getnames, args=(u, i))
           for i in range(3)]
    p.close()

    assert_equal(ref, res)


def test_creating_multiple_universe_without_offset(temp_xtc, ncopies=3):
    #  test if they can be created without generating
    #  the offset simultaneously.
    #  The tested XTC file is way too short to induce a race scenario but the
    #  test is included as documentation for the scenario that used to create
    #  a problem (see PR #3375 and issues #3230, #1988)

    args = (GRO, str(temp_xtc))
    with multiprocessing.Pool(2) as p:
        universes = [p.apply_async(mda.Universe, args) for i in range(ncopies)]
        universes = [universe.get() for universe in universes]


    assert_equal(universes[0].trajectory._xdr.offsets,
                 universes[1].trajectory._xdr.offsets)


@pytest.fixture(params=[
    # formatname, filename
    ('CRD', CRD, dict()),
    ('DATA', LAMMPSdata_mini, dict(n_atoms=1)),
    ('DCD', DCD, dict()),
    ('DMS', DMS, dict()),
    ('CONFIG', DLP_CONFIG, dict()),
    ('FHIAIMS', FHIAIMS, dict()),
    ('HISTORY', DLP_HISTORY, dict()),
    ('INPCRD', INPCRD, dict()),
    ('LAMMPSDUMP', LAMMPSDUMP, dict()),
    ('GMS', GMS_ASYMOPT, dict()),
    ('GRO', GRO, dict()),
    pytest.param(
        ('GSD', GSD, dict()),
        marks=pytest.mark.skipif(not HAS_GSD, reason='gsd not installed')
    ),
    ('MMTF', MMTF, dict()),
    ('MOL2', mol2_molecules, dict()),
    ('PDB', PDB_small, dict()),
    ('PQR', PQR, dict()),
    ('PDBQT', PDBQT_input, dict()),
    ('TRR', TRR, dict()),
    ('TRZ', TRZ, dict(n_atoms=8184)),
    ('TRJ', TRJ, dict(n_atoms=252)),
    ('XTC', XTC, dict()),
    ('XPDB', XPDB_small, dict()),
    ('XYZ', XYZ_mini, dict()),
    ('NCDF', NCDF, dict()),
    ('TXYZ', TXYZ, dict()),
    ('memory', np.arange(60).reshape(2, 10, 3).astype(np.float64), dict()),
    ('TRC', TRC_TRAJ1_VAC, dict()),
    ('CHAIN', [TRC_TRAJ1_VAC, TRC_TRAJ2_VAC], dict()),
    ('CHAIN', [GRO, GRO], dict()),
    ('CHAIN', [PDB, PDB], dict()),
    ('CHAIN', [XTC, XTC], dict()),
])
def ref_reader(request):
    fmt_name, filename, extras = request.param
    r = get_reader_for(filename, format=fmt_name)(filename, **extras)
    try:
        yield r
    finally:
        # make sure file handle is closed afterwards
        r.close()


def test_readers_pickle(ref_reader):
    ps = pickle.dumps(ref_reader)
    reanimated = pickle.loads(ps)
    assert len(ref_reader) == len(reanimated)
    try:
        ref_reader[2]
        ref_reader[0]
    except IndexError:
        # single frame files
        pass
    assert_equal(reanimated.ts, ref_reader.ts)


def test_analysis_pickle():
    u = mda.Universe(PSF, DCD)
    rmsd = RMSD(u.atoms, u.atoms)
    rmsd.run()
    rmsd_p = pickle.dumps(rmsd)
    rmsd_new = pickle.loads(rmsd_p)
