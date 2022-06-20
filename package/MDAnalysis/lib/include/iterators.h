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

#ifndef __ITERATORS_H
#define __ITERATORS_H

#include <cstdint>
#include <vector>
#include <cstdio>

/* This header provides classes to iterate groups (eg AtomGroup) and NumPy
   arrays using the same API. This is done by using thin wrapper classes with
   the same interface methods:

   - reset_iteration() -> rewind to start of group/array
   - load_into_external_buffer -> pass values from group/array into an external
     buffer (likely a stack allocated array)
   - seek -> go to a particular atom in the group/array

   Classes in this header do not inherit from an abstract base as this can
   become quite complicated when wrapped with Cython. Instead they are designed
   to be passed into templated C++ functions. For example:

   <template typename T,U,V,W>
   void calc_dihedrals(T g1, U g2, V g3, W g4)

   This has the added benefit of being able to accept groups or arrays at any
   position without writing overloads. However, in return all of the interface
   methods listed above must have the exact same signature and calling
   convention for all iterators.
*/

class _AtomGroupIterator
{
public:
    // number of atoms in the AtomGroup
    int64_t n_atoms;
    // array of indices in the AtomGroup
    std::vector<int64_t> ix;
    // internal index for iteration
    int64_t i;
    // pointer for coordinates, hooked onto coordinate memoryview in Cython layer
    float *ptr;

    // nullary constructor so is stack allocatable
    _AtomGroupIterator() : n_atoms(0), i(0), ptr(nullptr) {}

    explicit _AtomGroupIterator(int64_t n_atoms) : n_atoms(n_atoms), i(0), ptr(nullptr)
    {
    }
    // copy indices from AtomGroup
    // can possibly use pointer or std::reference wrapper here to avoid the copy
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
    // rewind to start of AtomGroup
    void inline reset_iteration()
    {
        i = 0;
    }
    // seek to atom i
    void seek(int64_t i_)
    {
        i = i_;
    }
    // load n_idx coordinate values into external buffer of size 3*n_atoms.
    // No checking done for maximal performance, callee's responsibility to not
    // overrun buffer or coordinate pointer (ptr).
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
    // number of atoms in the AtomGroup
    int64_t n_atoms;
    // internal index for iteration
    int64_t i;
    // pointer for coordinates, hooked onto coordinate memoryview in Cython layer
    float *ptr;

    // nullary constructor so is stack allocatable
    _ArrayIterator() : n_atoms(0), i(0), ptr(nullptr) {}

    explicit _ArrayIterator(int64_t n_atoms) : n_atoms(n_atoms), i(0), ptr(nullptr)
    {
    }
    // rewind to start of array
    void inline reset_iteration()
    {
        i = 0;
    }
    // seek to atom i
    void seek(int64_t i_)
    {
        i = i_;
    }
    // load coordinate values into external buffer. For an array this is done by
    // passing incoming pointer by reference and setting it equal to the correct
    // location of the coordinate pointer (ptr). No checking done for maximal
    // performance, callee's responsibility to not overrun buffer or
    // coordinate pointer (ptr).
    void load_into_external_buffer(float *&buffer, int64_t n_idx)
    {
        buffer = ptr + i * 3;
        i += n_idx;
    }
};
#endif //__ITERATORS_H
