/*
 MDAnalysis --- https://www.mdanalysis.org
 Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
 (see the file AUTHORS for the full list of names)

 Released under the Lesser GNU Public Licence, v2.1 or any higher version

 Please cite your use of MDAnalysis in published work:

 R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
 D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
 MDAnalysis: A Python package for the rapid analysis of molecular dynamics
 simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
 Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
 doi: 10.25080/majora-629e541a-00e

 N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
 MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
 J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
*/

#ifndef __DISTANCES_HPP__
#define __DISTANCES_HPP__

#include <limits>

template <typename T>
void _minimize_vectors_ortho(T* dx, T* box, T* inverse_box)
{
  /*
    Minimise dx to be the shortest vector
    Operates in-place on dx!

    This version is near identical to minimum_image, with the
    difference being that this version uses templatized types.
  */
  int i;
  T s;
  T eps = std::numeric_limits<T>::epsilon();
  for (i=0; i<3; i++) {
    if (box[i] > eps) {
      s = inverse_box[i] * dx[i];
      dx[i] = box[i] * (s - round(s));
    }
  }
}

template <typename T>
static void _calc_minimize_vectors_ortho(T* vectors, uint64_t numvectors, T* box, T* output)
{
  T inverse_box[3];
  inverse_box[0] = 1.0 / box[0];
  inverse_box[1] = 1.0 / box[1];
  inverse_box[2] = 1.0 / box[2];

#ifdef PARALLEL
#pragma omp parallel for shared(vectors)
#endif
  for (uint64_t i = 0; i < numvectors; i++) {
    T dx[3];
    dx[0] = vectors[i * 3];
    dx[1] = vectors[i * 3 + 1];
    dx[2] = vectors[i * 3 + 2];
    _minimize_vectors_ortho(dx, box, inverse_box);
    output[i * 3] = dx[0];
    output[i * 3 + 1] = dx[1];
    output[i * 3 + 2] = dx[2];
  }
}

template <typename T>
void _minimize_vectors_triclinic(T* dx, T* box, T* inverse_box)
{
  /*
    Minimise dx to be the shortest vector
    Operates in-place on dx!

    This version is near identical to minimum_image_triclinic, with the
    difference being that this version does not assume that coordinates are at
    most a single box length apart and uses templatized types.
  */
  T dx_min[3] = {0.0, 0.0, 0.0};
  T s;
  T dsq;
  T dsq_min = DBL_MAX;
  T rx;
  T ry[2];
  T rz[3];
  int ix, iy, iz;

  // Shift into primary unit cell
  s = round(inverse_box[2] * dx[2]);
  dx[0] -= s * box[6];
  dx[1] -= s * box[7];
  dx[2] -= s * box[8];
  s = round(inverse_box[1] * dx[1]);
  dx[0] -= s * box[3];
  dx[1] -= s * box[4];
  s = round(inverse_box[0] * dx[0]);
  dx[0] -= s * box[0];

  for (ix = -1; ix < 2; ++ix) {
    rx = dx[0] + box[0] * ix;
    for (iy = -1; iy < 2; ++iy) {
      ry[0] = rx + box[3] * iy;
      ry[1] = dx[1] + box[4] * iy;
      for (iz = -1; iz < 2; ++iz) {
        rz[0] = ry[0] + box[6] * iz;
        rz[1] = ry[1] + box[7] * iz;
        rz[2] = dx[2] + box[8] * iz;
        dsq = rz[0] * rz[0] + rz[1] * rz[1] + rz[2] * rz[2];
        if (dsq < dsq_min) {
            dsq_min = dsq;
            dx_min[0] = rz[0];
            dx_min[1] = rz[1];
            dx_min[2] = rz[2];
        }
      }
    }
  }
  dx[0] = dx_min[0];
  dx[1] = dx_min[1];
  dx[2] = dx_min[2];
}

template <typename T>
static void _calc_minimize_vectors_triclinic(T* vectors, uint64_t numvectors, T* box, T* output)
{
  T inverse_box[3];
  inverse_box[0] = 1.0 / box[0];
  inverse_box[1] = 1.0 / box[4];
  inverse_box[2] = 1.0 / box[8];

#ifdef PARALLEL
#pragma omp parallel for shared(vectors)
#endif
  for (uint64_t i = 0; i < numvectors; i++) {
    T dx[3];
    dx[0] = vectors[i * 3];
    dx[1] = vectors[i * 3 + 1];
    dx[2] = vectors[i * 3 + 2];
    _minimize_vectors_triclinic(dx, box, inverse_box);
    output[i * 3] = dx[0];
    output[i * 3 + 1] = dx[1];
    output[i * 3 + 2] = dx[2];
  }
}

#endif
