# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
cdef extern from *:
    ctypedef char const_char "const char"

cdef extern from "stdio.h":
    int printf(const_char *, ...)

cdef extern from "stdlib.h":
    ctypedef struct time_t:
        pass
    time_t time(time_t*)
    void*   malloc (size_t, size_t)
    void*   realloc (void*, size_t)
    void    srand(unsigned int)
    int     rand()

cdef extern from "math.h":
    double  sqrt(double)

cdef extern from "spe.h":
    ctypedef struct IVWrapper:
        pass
    ctypedef void* empty

    int trmIndex(int, int)
    double ed(double*, int, int, int)
    double stress(double, double, int, int)
    double neighbours_stress(double, double, int, int, double)
    int neighbours(double, int, double, int*, int*, int*)
    int* nearest_neighbours(double*, int, int)
    int cmp_ivwrapper(void*,void*)
    double CStochasticProximityEmbedding(double*, double*, double, int, int, double, double, int, int, int)
