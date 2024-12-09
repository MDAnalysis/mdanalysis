/* -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
 * vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
 *
 * MDAnalysis --- http://www.MDAnalysis.org
 * Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
 * Beckstein and contributors (see AUTHORS for the full list)
 *
 * Released under the Lesser GNU Public Licence, v2.1 or any higher version
 *
 * Please cite your use of MDAnalysis in published work:
 *
 * N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
 * MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
 * J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
 */

#ifndef _xtc_seek_h
#define _xtc_seek_h

#ifdef __cplusplus
extern "C" {
#endif

#include "xdrfile.h"

/* Seek through trajectory counting and indexing frames */
extern int read_xtc_n_frames(char *fn, int *n_frames, int *est_nframes,
                             int64_t **offsets);

/* XTC header fields until coord floats: *** only for trajectories of less than
 * 10 atoms! ***  */
/* magic natoms step time DIM*DIM_box_vecs natoms */
#define XTC_SHORTHEADER_SIZE (20 + DIM * DIM * 4)
/* Short XTCs store each coordinate as a 32-bit float. */
#define XTC_SHORT_BYTESPERATOM 12
/* XTC header fields until frame bytes: *** only for trajectories of more than 9
 * atoms! ***  */
/* magic natoms step time DIM*DIM_box_vecs natoms prec DIM_min_xyz DIM_max_xyz
 * smallidx */
#define XTC_HEADER_SIZE (DIM * DIM * 4 + DIM * 2 + 46)

#ifdef __cplusplus
}
#endif

#endif
