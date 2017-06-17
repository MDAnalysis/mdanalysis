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
    #extern int printf (__const char *__restrict __format, ...);
    int printf(const_char *, ...)

cdef extern from "stdlib.h":
    void *calloc (size_t COUNT, size_t ELTSIZE)

cdef extern from "float.h":
    enum:  FLT_MAX

cdef extern from "ap.h":
    int trmIndex(int, int)
    int sqmIndex(int, int, int)
    float pwmax(float, float)
    float pwmin(float, float)
    float min(float*, int)
    float max(float*, int)
    int CAffinityPropagation(float*, int, float, int, int, bint, long*)
