/* -*- tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
  vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

  MDAnalysis --- https://www.mdanalysis.org
  Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
  (see the file AUTHORS for the full list of names)

  Released under the Lesser GNU Public Licence, v2.1 or any higher version

  Please cite your use of MDAnalysis in published work:

  R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
  D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
  MDAnalysis: A Python package for the rapid analysis of molecular dynamics
  simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
  Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.

  N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
  MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
  J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
*/
int trmIndex(int, int);

int sqmIndex(int, int, int);

#ifndef _WIN32
    float min(float*, int);
    float max(float*, int);
#endif


int CAffinityPropagation(float*, int, float, int, int, int, long*);
