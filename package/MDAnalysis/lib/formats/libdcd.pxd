# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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


from libc.stdlib cimport free
from libc.stdint cimport uintptr_t
from libc.stdio cimport SEEK_SET, SEEK_CUR, SEEK_END
import cython

cimport numpy as cnp

cnp.import_array()


# Tell cython about the off_t type. It doesn't need to match exactly what is
# defined since we don't expose it to python but the cython compiler needs to
# know about it.
cdef extern from 'sys/types.h':
    ctypedef int off_t

ctypedef int fio_fd
ctypedef off_t fio_size_t

cdef extern from 'include/fastio.h':
    int fio_open(const char * filename, int mode, fio_fd * fd)
    int fio_fclose(fio_fd fd)
    fio_size_t fio_ftell(fio_fd fd)
    fio_size_t fio_fseek(fio_fd fd, fio_size_t offset, int whence)

cdef extern from 'include/readdcd.h':
    int read_dcdheader(fio_fd fd, int * natoms, int * nsets, int * istart,
                       int * nsavc, double * delta, int * nfixed, int ** freeind,
                       float ** fixedcoords, int * reverse_endian, int * charmm,
                       char ** remarks, int * len_remarks)
    void close_dcd_read(int * freeind, float * fixedcoords)
    int read_dcdstep(fio_fd fd, int natoms, float * X, float * Y, float * Z,
                     double * unitcell, int num_fixed,
                     int first, int * indexes, float * fixedcoords,
                     int reverse_endian, int charmm)
    int read_dcdsubset(fio_fd fd, int natoms, int lowerb, int upperb,
                       float * X, float * Y, float * Z,
                       double * unitcell, int num_fixed,
                       int first, int * indexes, float * fixedcoords,
                       int reverse_endian, int charmm)
    int write_dcdheader(fio_fd fd, const char * remarks, int natoms,
                        int istart, int nsavc, double delta, int with_unitcell,
                        int charmm)
    int write_dcdstep(fio_fd fd, int curframe, int curstep,
                      int natoms, const float * x, const float * y, const float * z,
                      const double * unitcell, int charmm)

cdef class DCDFile:
    """DCDFile(fname, mode='r')

    File like wrapper for DCD files

    This class can be similar to the normal file objects in python. The read()
    function will return a frame and all information in it instead of a single
    line. Additionally the context-manager protocol is supported as well.
    """
    # DCD file pointer
    cdef fio_fd fp
    # File name
    cdef readonly fname
    # Starting timestep of dcd file
    cdef int istart
    # Timesteps between dcd saves
    cdef int nsavc
    # Trajectory timestep
    cdef double delta
    # Number of atoms
    cdef int natoms
    # Number of fixed atoms
    cdef int nfixed
    # Free indices
    cdef int * freeind
    # Fixed coordinates
    cdef float * fixedcoords
    # Are we reverse endian?
    cdef int reverse_endian
    # Is the DCD file CHARMM style
    cdef int charmm
    # Is the file periodic
    cdef readonly int is_periodic
    # String data in the file
    cdef remarks
    # File mode
    cdef str mode
    # Number of dimensions in 
    cdef readonly int ndims
    # The number of frames
    cdef readonly int n_frames
    # Flag to indicate if header has been read
    cdef bint b_read_header
    # The current DCD frame
    cdef int current_frame
    # size of the first DCD frame
    cdef readonly fio_size_t _firstframesize
    # Size of a DCD frame 
    cdef readonly fio_size_t _framesize
    # Size of the DCD header
    cdef readonly fio_size_t _header_size
    # Is the file open?
    cdef int is_open
    # Have we reached the end of the file
    cdef int reached_eof
    # Have we written the header?
    cdef int wrote_header
    # Are buffers set up?
    cdef int _buffers_setup


    # buffer for reading coordinates
    cdef cnp.ndarray _coordinate_buffer
    # buffer for reading unitcell     
    cdef cnp.ndarray _unitcell_buffer

    # fortran contiguious memoryviews of the buffers to pass to the C code
    cdef float[::1] xview
    cdef float[::1] yview
    cdef float[::1] zview
    cdef double[::1] unitcellview

    cdef void _setup_buffers(self)

    cdef void _read_header(self)
    # Estimate the number of frames
    cdef int _estimate_n_frames(self)
    # Helper to read current DCD frame
    cdef int c_readframes_helper(self, float[::1] x,
                                 float[::1] y, float[::1] z,
                                 double[::1] unitcell, int first_frame)

# Helper in readframes to copy given a specific memory layout
cdef void copy_in_order(float[:, :] source, float[:, :, :] target, int order, int index)
