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
from __future__ import absolute_import
import re

import numpy as np
import pytest
from numpy.testing import assert_equal, assert_almost_equal

import MDAnalysis as mda
import MDAnalysis.analysis.helanal
from MDAnalysisTests.datafiles import (GRO, XTC, PSF, DCD, PDB_small,
                                       HELANAL_BENDING_MATRIX,
                                       HELANAL_BENDING_MATRIX_SUBSET)

# reference data from a single PDB file:
#   data = MDAnalysis.analysis.helanal.helanal_main(PDB_small,
#                    select="name CA and resnum 161-187")
HELANAL_SINGLE_DATA = {
    'Best fit tilt': 1.3309656332019535,
    'Height': np.array([ 1.5286051 ,  0.19648294,  0.11384312],
                       dtype=np.float32),
    'Local bending angles':
        np.array([  3.44526005,   4.85425806,   4.69548464,   2.39473653,
                    3.56172442,   3.97128344,   3.41563916,   1.86140978,
                    5.22997046,   5.41381264,  27.49601364,  39.69839478,
                    35.05921936,  21.78928566,   9.8632431 ,   8.80066967,
                    5.5344553 ,   6.14356709,  10.15450764,  11.07686138,
                    9.23541832], dtype=np.float32),
    'Residues/turn': np.array([ 3.64864163,  0.152694  ,  0.1131402 ]),
    'Rotation Angles':
        np.array([  87.80540079, -171.86019984,  -75.2341296 ,   24.61695962,
                    121.77104796, -134.94786976,  -35.07857424,   58.9621866 ,
                    159.40210233, -104.83368122,   -7.54816243,   87.40202629,
                    -176.13071955,  -89.13196878,   17.17321345,  105.26627814,
                    -147.00075298,  -49.36850775,   54.24557615,  156.06486532,
                    -110.82698327,   -5.72138626,   85.36050546, -167.28218858,
                    -68.23076936]),
    'Twist': np.array([ 98.83011627,   4.08171701,   3.07253003],
                      dtype=np.float32),
    'Unit twist angles':
        np.array([  97.23709869,   99.09676361,   97.25350952,  101.76019287,
                    100.42689514,   97.08784485,   97.07430267,   98.33553314,
                    97.86578369,   95.45792389,   97.10089111,   95.26415253,
                    87.93136597,  108.38458252,   95.27167511,  104.01931763,
                    100.7199707 ,  101.48034668,   99.64170074,   94.78946686,
                    102.50147247,   97.25154877,  104.54204559,  101.42829895],
                 dtype=np.float32),
    }


def read_bending_matrix(fn):
    """Read helanal_bending_matrix.dat into dict of numpy arrays.

    This is a quick and dirty parser with no error checking.

    Format::

      Mean
         0.0    5.7   10.9 ....
         .

      SD
        ....
        .

      ABDEV
       ....
       .

    Returns
    -------
       data : dict
           dictionary with keys 'Mean', 'SD', 'ABDEV' and NxN matrices.
    """
    data = {}
    current = None
    with open(fn) as bendingmatrix:
        for line in bendingmatrix:
            line = line.strip()
            if not line:
                continue
            if re.match("\D+", line):
                label = line.split()[0]
                current = data[label] = []
            else:
                current.append([float(x) for x in line.split()])
    for k,v in data.items():
        data[k] = np.array(v)
    return data


def test_helanal_trajectory(tmpdir, reference=HELANAL_BENDING_MATRIX,
                            outfile="helanal_bending_matrix.dat"):
    u = mda.Universe(PSF, DCD)
    with tmpdir.as_cwd():
        # Helix 8: 161 - 187 http://www.rcsb.org/pdb/explore.do?structureId=4AKE
        MDAnalysis.analysis.helanal.helanal_trajectory(
            u, select="name CA and resnum 161-187")
        bendingmatrix = read_bending_matrix(outfile)
        ref = read_bending_matrix(reference)
        assert_equal(sorted(bendingmatrix.keys()), sorted(ref.keys()),
                     err_msg="different contents in bending matrix data file")
        for label in ref.keys():
            assert_almost_equal(
                bendingmatrix[label], ref[label], err_msg="bending matrix "
                "stats for {0} mismatch".format(label))


def test_helanal_trajectory_slice(tmpdir,
                                  reference=HELANAL_BENDING_MATRIX_SUBSET,
                                  outfile="helanal_bending_matrix.dat"):
    """Simple regression test to validate that begin/finish work as
    intended. In this case, frames 10 (time: 10.999999031120204) through to
    79 (time: 79.99999295360149) should be picked up."""
    u = mda.Universe(PSF, DCD)
    with tmpdir.as_cwd():
        # Helix 8: 161 - 187 http://www.rcsb.org/pdb/explore.do?structureId=4AKE
        MDAnalysis.analysis.helanal.helanal_trajectory(
            u, select="name CA and resnum 161-187", begin=10, finish=80
        )
        bendingmatrix = read_bending_matrix(outfile)
        ref = read_bending_matrix(reference)
        assert_equal(sorted(bendingmatrix.keys()), sorted(ref.keys()),
                     err_msg="different contents in bending matrix data file")
        for label in ref.keys():
            assert_almost_equal(
                bendingmatrix[label], ref[label], err_msg="bending matrix "
                "starts for {0} mismatch".format(label))


def test_helanal_main(reference=HELANAL_SINGLE_DATA):
    u = mda.Universe(PDB_small)
    # Helix 8: 161 - 187 http://www.rcsb.org/pdb/explore.do?structureId=4AKE
    data = MDAnalysis.analysis.helanal.helanal_main(
        PDB_small, select="name CA and resnum 161-187")
    ref = reference
    assert_equal(sorted(data.keys()), sorted(ref.keys()),
                     err_msg="different contents in data dict")
    for label in ref.keys():
        assert_almost_equal(data[label], ref[label], decimal=4,
                                  err_msg="data[{0}] mismatch".format(label))


def test_exceptions(tmpdir):
    """Testing exceptions which can be raised"""
    u = MDAnalysis.Universe(GRO, XTC)
    u.trajectory[1]

    # Testing xtc striding: Check for resolution of Issue #188
    with tmpdir.as_cwd():
        with pytest.raises(ValueError):
            MDAnalysis.analysis.helanal.helanal_trajectory(
                u, select="name CA", finish=5
            )

    with tmpdir.as_cwd():
        with pytest.raises(ValueError):
            MDAnalysis.analysis.helanal.helanal_trajectory(
                u, select="name CA", begin=1, finish=0
            )

    with tmpdir.as_cwd():
        with pytest.raises(ValueError):
            MDAnalysis.analysis.helanal.helanal_trajectory(
                u, select="name CA", begin=99999
            )


def test_warnings(tmpdir):
    """Testing that a warning which can be raised"""
    u = MDAnalysis.Universe(GRO, XTC)
    u.trajectory[1]

    with tmpdir.as_cwd():
        with pytest.warns(UserWarning) as record:
            MDAnalysis.analysis.helanal.helanal_trajectory(
                u, select="name CA", begin=-1, finish=99999
            )

            assert len(record) == 3
            wmsg1 = ("The input begin time (-1 ps) precedes the starting "
                   "trajectory time --- Setting starting frame to 0".format(
                    -1,0))
            assert str(record[1].message.args[0]) == wmsg1
            wmsg2 = ("The input finish time ({0} ps) occurs after the end of "
                     "the trajectory ({1} ps). Finish time will be set to "
                     "the end of the trajectory".format(
                      99999,1000.0000762939453))
            assert str(record[2].message.args[0]) == wmsg2
