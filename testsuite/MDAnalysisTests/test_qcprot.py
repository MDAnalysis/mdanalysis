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

"""
:Author:   Joshua L. Adelman, University of Pittsburgh
:Contact:  jla65@pitt.edu

Sample code to use the routine for fast RMSD & rotational matrix calculation.
For the example provided below, the minimum least-squares RMSD for the two
7-atom fragments should be 0.719106 A.

    And the corresponding 3x3 rotation matrix is:

    [[ 0.72216358 -0.52038257 -0.45572112]
     [ 0.69118937  0.51700833  0.50493528]
     [-0.0271479  -0.67963547  0.73304748]]

"""

import numpy as np

import MDAnalysis.lib.qcprot as qcp

from numpy.testing import assert_almost_equal
from nose.plugins.attrib import attr


# Calculate rmsd after applying rotation
def rmsd(a, b):
    """Returns RMSD between two coordinate sets a and b."""
    return np.sqrt(np.sum(np.power(a - b, 2)) / a.shape[1])


def test_CalcRMSDRotationalMatrix():
    # Setup coordinates
    frag_a = np.zeros((3, 7), dtype=np.float64)
    frag_b = np.zeros((3, 7), dtype=np.float64)
    N = 7

    frag_a[0][0] = -2.803
    frag_a[1][0] = -15.373
    frag_a[2][0] = 24.556
    frag_a[0][1] = 0.893
    frag_a[1][1] = -16.062
    frag_a[2][1] = 25.147
    frag_a[0][2] = 1.368
    frag_a[1][2] = -12.371
    frag_a[2][2] = 25.885
    frag_a[0][3] = -1.651
    frag_a[1][3] = -12.153
    frag_a[2][3] = 28.177
    frag_a[0][4] = -0.440
    frag_a[1][4] = -15.218
    frag_a[2][4] = 30.068
    frag_a[0][5] = 2.551
    frag_a[1][5] = -13.273
    frag_a[2][5] = 31.372
    frag_a[0][6] = 0.105
    frag_a[1][6] = -11.330
    frag_a[2][6] = 33.567

    frag_b[0][0] = -14.739
    frag_b[1][0] = -18.673
    frag_b[2][0] = 15.040
    frag_b[0][1] = -12.473
    frag_b[1][1] = -15.810
    frag_b[2][1] = 16.074
    frag_b[0][2] = -14.802
    frag_b[1][2] = -13.307
    frag_b[2][2] = 14.408
    frag_b[0][3] = -17.782
    frag_b[1][3] = -14.852
    frag_b[2][3] = 16.171
    frag_b[0][4] = -16.124
    frag_b[1][4] = -14.617
    frag_b[2][4] = 19.584
    frag_b[0][5] = -15.029
    frag_b[1][5] = -11.037
    frag_b[2][5] = 18.902
    frag_b[0][6] = -18.577
    frag_b[1][6] = -10.001
    frag_b[2][6] = 17.996

    # Allocate rotation array
    rot = np.zeros((9,), dtype=np.float64)

    # Calculate center of geometry
    comA = np.sum(frag_a, axis=1) / N
    comB = np.sum(frag_b, axis=1) / N

    # Center each fragment
    frag_a = frag_a - comA.reshape(3, 1)
    frag_b = frag_b - comB.reshape(3, 1)

    # Calculate rmsd and rotation matrix
    qcp_rmsd = qcp.CalcRMSDRotationalMatrix(frag_a, frag_b, N, rot, None)

    #print 'qcp rmsd = ',rmsd
    #print 'rotation matrix:'
    #print rot.reshape((3,3))

    # rotate frag_b to obtain optimal alignment
    frag_br = frag_b.T * np.matrix(rot.reshape((3, 3)))
    aligned_rmsd = rmsd(frag_br.T, frag_a)
    #print 'rmsd after applying rotation: ',rmsd

    assert_almost_equal(aligned_rmsd, 0.719106, 6, "RMSD between fragments A and B does not match excpected value.")

    expected_rot = np.array([
        [0.72216358, -0.52038257, -0.45572112],
        [0.69118937, 0.51700833, 0.50493528],
        [-0.0271479, -0.67963547, 0.73304748]])
    assert_almost_equal(rot.reshape((3, 3)), expected_rot, 6,
                        "Rotation matrix for aliging B to A does not have expected values.")
