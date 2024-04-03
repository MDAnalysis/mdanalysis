# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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

import numpy as np
cimport numpy as cnp

cnp.import_array()
from libc.stdint cimport int64_t

from libc.stdio cimport SEEK_SET, SEEK_CUR, SEEK_END

cdef extern from 'include/xdrfile.h':
    ctypedef struct XDRFILE:
        pass

    XDRFILE* xdrfile_open (char * path, char * mode)
    int xdrfile_close (XDRFILE * xfp)
    int xdr_seek(XDRFILE *xfp, int64_t pos, int whence)
    int64_t xdr_tell(XDRFILE *xfp)
    ctypedef float matrix[3][3]
    ctypedef float rvec[3]


cdef extern from 'include/xdrfile_xtc.h':
    int read_xtc_natoms(char * fname, int * natoms)
    int read_xtc(XDRFILE * xfp, int natoms, int * step, float * time, matrix box,
                 rvec * x, float * prec)
    int write_xtc(XDRFILE * xfp, int natoms, int step, float time, matrix box,
                  rvec * x, float prec)



cdef extern from 'include/xdrfile_trr.h':
    int read_trr_natoms(char *fname, int *natoms)
    int read_trr(XDRFILE *xfp, int natoms, int *step, float *time, float *_lambda,
                 matrix box, rvec *x, rvec *v, rvec *f, int *has_prop)
    int write_trr(XDRFILE *xfp, int natoms, int step, float time, float _lambda,
                  matrix box, rvec *x, rvec *v, rvec *f)


cdef extern from 'include/xtc_seek.h':
    int read_xtc_n_frames(char *fn, int *n_frames, int *est_nframes, int64_t **offsets)


cdef extern from 'include/trr_seek.h':
    int read_trr_n_frames(char *fn, int *n_frames, int *est_nframes, int64_t **offsets)


cdef class _XDRFile:
    """Base python wrapper for gromacs-xdr formats

    This class can be similar to the normal file objects in python. The read()
    function will return a frame and all information in it instead of a single
    line. Additionally the context-manager protocoll is supported as well.

    Note
    ----
    This class can't be initialized; use one of the subclasses :class:`XTCFile` or :class:`TRRFile`.
    """
    # number of atoms
    cdef readonly int n_atoms
    # whether the file is open
    cdef int is_open
    # whether we have reached the end of the file
    cdef int reached_eof
    # the XDR file pointer
    cdef XDRFILE *xfp
    # the name of the xdr file 
    cdef readonly fname
    # the current frame in the file
    cdef int current_frame
    # the file mode
    cdef str mode
    # the simulation box
    cdef cnp.ndarray box
    # numpy array of offsets into the fle 
    cdef cnp.ndarray _offsets
    # whether we have the offsets
    cdef readonly int _has_offsets


cdef class XTCFile(_XDRFile):
    """File-like wrapper for gromacs XTC files."""
    # precision of the XTC file
    cdef float precision



cdef class TRRFile(_XDRFile):
    """File-like wrapper for gromacs TRR files"""
