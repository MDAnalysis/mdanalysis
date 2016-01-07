# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
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
import os
import re
import tempdir

import numpy as np
from numpy.testing import (assert_raises, assert_,
                           assert_equal, assert_array_almost_equal)
from six.moves import zip

import MDAnalysis as mda
import MDAnalysis.analysis.helanal
from MDAnalysis import FinishTimeException
from MDAnalysisTests.datafiles import GRO, XTC, PSF, DCD, HELANAL_BENDING_MATRIX


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


def test_helanal_trajectory(reference=HELANAL_BENDING_MATRIX,
                            outfile="helanal_bending_matrix.dat"):
    u = mda.Universe(PSF, DCD)
    with tempdir.in_tempdir():
        # Helix 8: 161 - 187 http://www.rcsb.org/pdb/explore.do?structureId=4AKE
        MDAnalysis.analysis.helanal.helanal_trajectory(u,
                                                       selection="name CA and resnum 161-187")
        bendingmatrix = read_bending_matrix(outfile)
        ref = read_bending_matrix(reference)
        assert_equal(sorted(bendingmatrix.keys()), sorted(ref.keys()),
                     err_msg="different contents in bending matrix data file")
        for label in ref.keys():
            assert_array_almost_equal(bendingmatrix[label], ref[label],
                                      err_msg="bending matrix stats for {0} mismatch".format(label))


def test_xtc_striding():
    """testing MDAnalysis.analysis.helanal xtc striding: Check for resolution of Issue #188."""
    u = MDAnalysis.Universe(GRO, XTC)
    u.trajectory[1]

    with tempdir.in_tempdir():
        assert_raises(FinishTimeException,
                      MDAnalysis.analysis.helanal.helanal_trajectory,
                      u, selection="name CA", finish=5)

    #with assert_raises(FinishTimeException):
    #    try:
    #        MDAnalysis.analysis.helanal.helanal_trajectory(u, selection=sel, finish=5)
    #   except IndexError:
    #       self.fail("IndexError consistent with Issue 188.")


