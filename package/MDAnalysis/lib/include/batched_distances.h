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

    // global counters
    uint64_t iter_ref = 0;
    uint64_t iter_conf = 0;
    
    // internals
    uint64_t i, j;
    double rsq;
    float dx[3];

    // what is modulo number of particles?  
    int ref_overhang = nref % bsize_ref;
    int conf_overhang = nconf % bsize_conf;

    for (; iter_ref < nref - ref_overhang; iter_ref += bsize_ref)
    {
        ref.load_into_external_buffer(ref_buffer_, bsize_ref);

        for (; iter_conf < nconf - conf_overhang; iter_conf += bsize_conf)
        {
            conf.load_into_external_buffer(conf_buffer_, bsize_conf);

            for (i = 0; i < bsize_ref; i++)
            {
                for (j = 0; j < bsize_conf; j++)
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
        iter_conf = 0;
    }
}

#endif // __BATCHED_DISTANCES_H