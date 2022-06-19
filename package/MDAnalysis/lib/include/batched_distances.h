/* -*- Mode: C; tab-width: 4; indent-tabs-mode:nil; -*- */
/* vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 */
/*
  MDAnalysis --- http://mdanalysis.googlecode.com

  Copyright (c) 2006-2022 Naveen Michaud-Agrawal,
                Elizabeth J. Denning, Oliver Beckstein,
                and contributors (see AUTHORS for the full list)
  Released under the GNU Public Licence, v2 or any higher version

  Please cite your use of MDAnalysis in published work:

      N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
      O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
      Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
      in press.
*/

#ifndef __BATCHED_DISTANCES_H
#define __BATCHED_DISTANCES_H

#include <cstdint>
#include <numeric>
#include "iterators.h"

/* Functions in this header file calculate distances using batched algorithms
 with stack allocated memory owned by the function. Testing indicated that
 batching the calculation was much more performant than using a `.next()`
 style API, with contiguity of memory access of paramount importance.
 Optimal sizing of the batching is yet to be fully explored.

 Functions in this header accept the iterators defined in `iterators.h`
 (currently _AtomGroupIterator and _ArrayIterator) that provide a homogenous
 interface for iterating through their respective datatypes
 (AtomGroups and NumPy ndarrays).

 Templates are used to allow iterators at either position, avoiding writing
 multiple overloads for functions that accept many iterators.

*/

/* _calc_distance_array_batched uses the following batched algorithm:
   - figure out how many perfect tiles of size (batchsize, batchsize) we can do
   - figure out how overhanging our remaining distances are
   - figure out the maximum number of coordinates at a time for overhang (gcd)
   - calculate distances in tiles
   - if overhanging in ref dimension, do overhang NOTE: non contiguous access
   - if overhanging in conf dimension, do overhang

   NOTE: possible poor performance on overhanging portion when gcd is low (ie a
   prime number of distances). Should be investigated.
 */

// implement gcd because we don't compile with C++17
template <typename T>
T _gcd(T a, T b)
{
    if (b == 0)
        return a;
    return _gcd(b, a % b);
}

template <typename T, typename U>
void _calc_distance_array_batched(T ref, U conf, double *distances, uint64_t batchsize)
{
    // batchsize is in particles, so allocate buffers of 3x batchsize.
    const uint64_t atom_bufsize = 3 * batchsize;
    float ref_buffer[atom_bufsize];
    float conf_buffer[atom_bufsize];

    // slight indirection required here to homogenize interface between
    // _AtomGroupIterator and _ArrayIterator where passing stack allocated
    // array as float*& does not decay to a float* as desired.
    float *ref_buffer_ = ref_buffer;
    float *conf_buffer_ = conf_buffer;

    uint64_t nref = ref.n_atoms;
    uint64_t nconf = conf.n_atoms;

    // check if batchsize is larger than the number of coordinates and if so
    // truncate.
    uint64_t bsize_ref = std::min(nref, batchsize);
    uint64_t bsize_conf = std::min(nconf, batchsize);

    // internals
    double rsq;
    float dx[3];
    uint64_t iter_ref;
    uint64_t iter_conf;

    // what is modulo number of particles?
    uint64_t ref_overhang = nref % bsize_ref;
    uint64_t conf_overhang = nconf % bsize_conf;

    // what is the gcd to do iteration of overhang?
    // using a gcd allows us to use the maximum size of memory buffer that fits
    // remaining iteration in (gcd_ref, gcd_conf) sized tiles
    int64_t gcd_conf = _gcd(nconf, bsize_conf);
    int64_t gcd_ref = _gcd(nref, bsize_ref);

    for (iter_ref = 0; iter_ref < nref - ref_overhang; iter_ref += bsize_ref)
    {
        ref.load_into_external_buffer(ref_buffer_, bsize_ref);

        for (iter_conf = 0; iter_conf < nconf - conf_overhang; iter_conf += bsize_conf)
        {
            conf.load_into_external_buffer(conf_buffer_, bsize_conf);

            for (uint64_t i = 0; i < bsize_ref; i++)
            {
                for (uint64_t j = 0; j < bsize_conf; j++)
                {
                    dx[0] = conf_buffer_[3 * j] - ref_buffer_[3 * i];
                    dx[1] = conf_buffer_[3 * j + 1] - ref_buffer_[3 * i + 1];
                    dx[2] = conf_buffer_[3 * j + 2] - ref_buffer_[3 * i + 2];
                    rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
                    *(distances + iter_ref * nconf + iter_conf + i * nconf + j) = sqrt(rsq);
                }
            }
        }

        conf.reset_iteration();
    }

    if (ref_overhang) // technically not required but may help with branch prediction
    {
        // deal with overhang in ref dimension, contiguous dim
        // if ref is overhanging we enter this block with just the overhang left
        // we enter this block with conf already reset
        ref.load_into_external_buffer(ref_buffer_, ref_overhang);

        for (uint64_t i = 0; i < nconf; i += gcd_conf)
        {
            conf.load_into_external_buffer(conf_buffer_, gcd_conf);
            for (uint64_t j = 0; j < ref_overhang; j++)
            {
                for (uint64_t k = 0; k < gcd_conf; k++)
                {
                    dx[0] = conf_buffer_[3 * k] - ref_buffer_[3 * j];
                    dx[1] = conf_buffer_[3 * k + 1] - ref_buffer_[3 * j + 1];
                    dx[2] = conf_buffer_[3 * k + 2] - ref_buffer_[3 * j + 2];
                    rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
                    *(distances + iter_ref * nconf + j * nconf + i + k) = sqrt(rsq);
                }
            }
        }
    }

    if (conf_overhang) // technically not required but may help with branch prediction
    {
        // deal with overhang in the conf dimension, strided dim
        // we need to rewind ref
        ref.reset_iteration();
        // we are at the beginning or end of conf depending if we had an
        // overhang or not and need to rewind
        conf.seek(nconf - conf_overhang);
        conf.load_into_external_buffer(conf_buffer_, conf_overhang);

        for (uint64_t i = 0; i < nref; i += gcd_ref)
        {
            ref.load_into_external_buffer(ref_buffer_, gcd_ref);
            for (uint64_t j = 0; j < conf_overhang; j++)
            {
                for (uint64_t k = 0; k < gcd_ref; k++)
                {
                    dx[0] = conf_buffer_[3 * j] - ref_buffer_[3 * k];
                    dx[1] = conf_buffer_[3 * j + 1] - ref_buffer_[3 * k + 1];
                    dx[2] = conf_buffer_[3 * j + 2] - ref_buffer_[3 * k + 2];
                    rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
                    *(distances + iter_conf + i * nconf + k * nconf + j) = sqrt(rsq);
                }
            }
        }
    }
}

#endif // __BATCHED_DISTANCES_H