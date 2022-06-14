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
#include <cstdio>

class _AtomGroupIterator
{
public:
    int64_t n_atoms;
    std::vector<int64_t> ix;
    int64_t i;
    float *ptr;

    _AtomGroupIterator() : n_atoms(0), i(0), ptr(nullptr) {}

    explicit _AtomGroupIterator(int64_t n_atoms) : n_atoms(n_atoms), i(0), ptr(nullptr)
    {
    }

    // can possibly use std::reference wrapper here to avoid the copy
    void copy_ix(const int64_t *source)
    {
        std::vector<int64_t> tmp(source, source + n_atoms);
        ix = tmp;
    }

    void print_ix()
    {
        for (int64_t n = 0; n < n_atoms; n++)
        {
            printf("ix %ld val %ld \n", n, ix[n]);
        }
    }

    void inline reset_iteration()
    {
        i = 0;
    }

    void seek(int64_t i)
    {
        i = i;
    }


    void load_into_external_buffer(float *buffer, int64_t n_idx)
    {
        for (int64_t n = 0; n < n_idx; n++)
        {
            buffer[3 * n] = ptr[3 * ix[i]];
            buffer[3 * n + 1] = ptr[3 * ix[i] + 1];
            buffer[3 * n + 2] = ptr[3 * ix[i] + 2];
            i += 1;
        }
    }
};

class _ArrayIterator
{
public:
    int64_t n_atoms;
    int64_t i;
    float *ptr;

    _ArrayIterator() : n_atoms(0), i(0), ptr(nullptr) {}

    explicit _ArrayIterator(int64_t n_atoms) : n_atoms(n_atoms), i(0), ptr(nullptr)
    {
    }

    void inline reset_iteration()
    {
        i = 0;
    }

    void seek(int64_t i)
    {
        i = i;
    }

    void load_into_external_buffer(float *buffer, int64_t n_idx)
    {
        buffer = ptr;
    }
};
#endif //__WRAPPER_CLASSES_H
