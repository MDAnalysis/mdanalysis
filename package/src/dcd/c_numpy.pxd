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

cdef extern from "Python.h":
                ctypedef int Py_intptr_t

cdef extern from "numpy/arrayobject.h":
                ctypedef class numpy.ndarray [object PyArrayObject]:
                                cdef char *data
                                cdef int nd
                                cdef Py_intptr_t *dimensions
                                cdef Py_intptr_t *strides
                                cdef object base
                                # descr not implemented yet here...
                                cdef int flags
                                cdef int itemsize
                                cdef object weakreflist

                cdef void import_array()
