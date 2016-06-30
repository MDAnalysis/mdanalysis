# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

cdef extern from 'include/fastio.h':
    ctypedef fio_fd:
    static int fio_open(const char *filename, int mode, fio_fd *fd)
    static int fio_fclose(fio_fd fd)

#define FIO_READ  0x01
#define FIO_WRITE 0x02

cdef enum:
    FIO_READ = 0x01
    FIO_WRITE = 0x02

cdef class DCDFile:
    cdef fio_fd fp
    cdef readonly fname
    cdef int is_open

    def __cinit__(self, fname, mode='r'):
        self.fname = fname.encode('utf-8')
        self.is_open = False
        self.open(self.fname, mode)


    def __dealloc__(self):
        # call a close_dcd_read
        self.close()

    def open(self, filename, mode):
        # NOTE: to make handling easier lets disallow read/write mode
        if mode == 'r':
            fio_mode = FIO_READ
        elif mode == 'w':
            fio_mode = FIO_WRITE
        else:
            raise IOError("unkown mode '{}', use either r or w".format(mode))
        ok = fio_open(self.fname, fio_mode, self.fp)
        if ok != 0:
            raise IOError("couldn't open file: {}".format(filename))

    def close(self):
        ok = fio_close(self.fp)
        if ok != 0:
            raise IOError("couldn't close file: {}".format(self.fname))
