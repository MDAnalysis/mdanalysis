/* -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
 * vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
 *
 * MDAnalysis --- http://www.MDAnalysis.org
 * Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
 * Beckstein and contributors (see AUTHORS for the full list)
 *
 * Released under the GNU Public Licence, v2 or any higher version
 *
 * Please cite your use of MDAnalysis in published work:
 *
 * N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
 * MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
 * J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
 */

#ifndef _trr_seek_h
#define _trr_seek_h

#ifdef __cplusplus
extern "C" {
#endif

#include "xdrfile.h"

/* Skip through trajectory, reading headers, obtain the total number of frames
 * in the trr */
extern int read_trr_n_frames(char *fn, int *n_frames, int *est_nframes,
                             int64_t **offsets);

/* Minimum TRR header size. It can have 8 bytes more if we have double time and
 * lambda. */
#define TRR_MIN_HEADER_SIZE 54
#define TRR_DOUBLE_XTRA_HEADER 8

#ifdef __cplusplus
}
#endif

#endif
