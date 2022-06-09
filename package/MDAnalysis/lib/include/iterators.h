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

#ifndef __WRAPPER_CLASSES_H
#define __WRAPPER_CLASSES_H

#include <cstdint>
#include <vector>


class _AtomGroupIterator
{
public:
    uint64_t n_atoms;
    std::vector<uint64_t> ix;
    uint64_t i;
    float *ptr;

    _AtomGroupIterator() {}

    explicit _AtomGroupIterator(uint64_t n_atoms) : n_atoms(n_atoms), i(0), ptr(nullptr)
    {
    }

    void inline reset_iteration()
    {
        i = 0;
    }

    void seek(int i)
    {
        i = i;
    }

    void load_into_external_buffer(float *buffer, uint64_t n_idx)
    {
        for (uint64_t i = 0; i < n_idx; i++)
        {
            // buffer[3 * i] = coords[3 * ix[i_preload]];
            // buffer[3 * i + 1] = coords[3 * ix[i_preload] + 1];
            // buffer[3 * i + 2] = coords[3 * ix[i_preload] + 2];
            // i += 1;
        }
    }
};

class _ArrayIterator
{
public:
    uint64_t n_atoms;
    uint64_t i;
    float *ptr;

    _ArrayIterator() {}

    explicit _ArrayIterator(uint64_t n_atoms) : n_atoms(n_atoms), i(0), ptr(nullptr)
    {
    }

    void inline reset_iteration()
    {
        i = 0;
    }

    void seek(int i)
    {
        i = i;
    }

    void load_into_external_buffer(float *buffer, uint64_t n_idx)
    {
        buffer = ptr;
    }
};
#endif //__WRAPPER_CLASSES_H
