/* -*- Mode: C; tab-width: 4; indent-tabs-mode:nil; -*- */
/* vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 */
/*
  MDAnalysis --- http://mdanalysis.googlecode.com

  Copyright (c) 2006-2014 Naveen Michaud-Agrawal,
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
#include "wrapper_classes.h"

template <typename T, typename U>
void _calc_distance_array_batched(T ref, U conf, double *distances, uint64_t batchsize)
{

    const uint64_t atom_bufsize = 3 * batchsize;
    float ref_buffer[atom_bufsize];
    float conf_buffer[atom_bufsize];

    uint64_t nref = ref->n_atoms;
    uint64_t nconf = conf->n_atoms;
    uint64_t bsize_ref = std::min(nref, batchsize); // is  our batchsize larger than the number of coords?
    uint64_t bsize_conf = std::min(nconf, batchsize);


    uint64_t iter_ref = 0;
    uint64_t iter_conf = 0;
    uint64_t i, j;
    double rsq;
    float dx[3];
    int ref_overhang = nref % bsize_ref;
    int conf_overhang = nconf % bsize_conf;

    for (; iter_ref < nref - ref_overhang; iter_ref += bsize_ref)
    {
        ref.load_into_external_buffer(ref_buffer, bsize_ref);

        for (; iter_conf < nconf - conf_overhang; iter_conf += bsize_conf)
        {
            conf.load_into_external_buffer(conf_buffer, bsize_conf);

            for (i = 0; i < bsize_ref; i++)
            {
                for (j = 0; j < bsize_conf; j++)
                {
                    dx[0] = conf_buffer[3 * j] - ref_buffer[3 * i];
                    dx[1] = conf_buffer[3 * j + 1] - ref_buffer[3 * i + 1];
                    dx[2] = conf_buffer[3 * j + 2] - ref_buffer[3 * i + 2];
                    rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
                    *(distances + iter_ref * nconf + iter_conf + i * nconf + j) = sqrt(rsq);
                }
            }
        }

        reset_iteration(conf);
        iter_conf = 0;
    }
}


#endif // __BATCHED_DISTANCES_H