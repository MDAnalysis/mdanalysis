# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2014 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see AUTHORS for the full list)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

cimport c_numpy

cdef extern from "clapack.h":
    int dsyev_(char *jobz, char *uplo, int *n, double *fa, int *lda, double *w, double *work, int *lwork, int *info)


cdef extern from "string.h":
    void* memcpy(void *dst, void *src, int len)

cdef eigen4x4sym(double* mat, double* eval, double* evec):
    cdef double work[264]
    cdef int n, lwork, info

    # Copy mat data to evec since dsyev places the eigenvectors in the original matrix
    memcpy(<void*>evec, <void*>mat, sizeof(double)*16)

    n = 4
    lwork = 264
    info = -1

    dsyev_("V", "U", &n, evec, &n, eval, work, &lwork, &info)
    return info

c_numpy.import_array()

import numpy as np
import warnings
def rms_rotation_matrix(c_numpy.ndarray conf, c_numpy.ndarray ref, c_numpy.ndarray weights):
    ''' Computes the weighted RMS rotation matrix between a reference set of coordinates and a comparison set.
        conf - Nx3 array of comparison coordinates (size float - 32 bits)
        ref - Nx3 array of reference coordinates (size float - 32 bits)
        weights - N array of masses (size double - 64 bits)

        Returns a numpy.matrix corresponding to the rotation matrix
        Note: Due to numpy's broadcasting rules, you should  multiply the vector with the matrix
        ie   coor * R, where coor is the coordinate array and R the rotation matrix'''
    cdef c_numpy.ndarray ref_cms, pos, cross, k, mat
    cdef c_numpy.ndarray e, v
    cdef int i, n_atoms
    cdef double possq

    warnings.warn("rms_rotation_matrix() is deprecated and will be removed in 0.8. Use MDAnalysis.core.qcprot.CalcRMSDRotationalMatrix()",
                  category=DeprecationWarning)

    if conf.dimensions[0] != ref.dimensions[0]:
        raise Exception("Error: RMS fit - conformation and reference don't have the same number of atoms")
    n_atoms = conf.dimensions[0]
    ref_cms = np.zeros((3,), np.float64)
    for i from 0 <= i < n_atoms:
        ref_cms = ref_cms + ref[i]
    ref_cms = ref_cms / n_atoms
    pos = np.zeros((3,), np.float64)
    possq = 0.
    cross = np.zeros((3,3), np.float64)
    for i from 0 <= i < n_atoms:
        r = conf[i]
        r_ref = ref[i]-ref_cms
        w = weights[i]
        pos = pos + w*r
        possq = possq + w*np.add.reduce(r*r) + w*np.add.reduce(r_ref*r_ref)
        cross = cross + w*r[:, np.newaxis]*r_ref[np.newaxis, :]
    k = np.zeros((4,4), np.float64)
    k[0, 0] = -cross[0, 0]-cross[1, 1]-cross[2, 2]
    k[0, 1] = cross[1, 2]-cross[2, 1]
    k[0, 2] = cross[2, 0]-cross[0, 2]
    k[0, 3] = cross[0, 1]-cross[1, 0]
    k[1, 1] = -cross[0, 0]+cross[1, 1]+cross[2, 2]
    k[1, 2] = -cross[0, 1]-cross[1, 0]
    k[1, 3] = -cross[0, 2]-cross[2, 0]
    k[2, 2] = cross[0, 0]-cross[1, 1]+cross[2, 2]
    k[2, 3] = -cross[1, 2]-cross[2, 1]
    k[3, 3] = cross[0, 0]+cross[1, 1]-cross[2, 2]
    for i from 1 <= i < 4:
        for j from 0 <= j < i:
            k[i, j] = k[j, i]
    k = 2.*k
    for i from 0 <= i < 4:
        k[i, i] = k[i, i] + possq - np.add.reduce(pos*pos)
    e = np.zeros((4,), np.float64)
    v = np.zeros((4,4), np.float64)
    i = eigen4x4sym(<double*>k.data, <double*>e.data, <double*>v.data)
    v = v[0]
    # Quaternion def: v = qw + iqx + jqy + kqz
    mat = np.zeros((3,3), np.float64)
    mat[0,0] = 1 - 2*v[2]*v[2] - 2*v[3]*v[3]
    mat[0,1] = 2*v[1]*v[2] - 2*v[0]*v[3]
    mat[0,2] = 2*v[1]*v[3] + 2*v[0]*v[2]
    mat[1,0] = 2*v[1]*v[2] + 2*v[0]*v[3]
    mat[1,1] = 1 - 2*v[1]*v[1] - 2*v[3]*v[3]
    mat[1,2] = 2*v[2]*v[3] - 2*v[0]*v[1]
    mat[2,0] = 2*v[1]*v[3] - 2*v[0]*v[2]
    mat[2,1] = 2*v[2]*v[3] + 2*v[0]*v[1]
    mat[2,2] = 1 - 2*v[1]*v[1] - 2*v[2]*v[2]
    return np.matrix(mat, copy=False)
