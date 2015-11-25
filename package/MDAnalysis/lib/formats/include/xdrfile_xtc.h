/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*-
 *
 * $Id$
 *
 /* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*-
 *
 * $Id$
 *
 * Copyright (c) 2009-2014, Erik Lindahl & David van der Spoel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 this
 * list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _xdrfile_xtc_h
#define _xdrfile_xtc_h

#ifdef __cplusplus
extern "C" {
#endif

#include "xdrfile.h"

/* All functions return exdrOK if succesfull.
 * (error codes defined in xdrfile.h).
 */

/* This function returns the number of atoms in the xtc file in *natoms */
extern int read_xtc_natoms(char *fn, int *natoms);

/* Seek through trajectory counting and indexing frames */
extern int read_xtc_n_frames(char *fn, int *n_frames, int *est_nframes,
                             int64_t **offsets);

  /* Read one frame of an open xtc file */
  extern int read_xtc(XDRFILE * xd, int natoms, int *step, float *time,
                      matrix box, rvec *x, float *prec);

  /* Write a frame to xtc file */
  extern int write_xtc(XDRFILE * xd, int natoms, int step, float time,
                       matrix box, rvec *x, float prec);

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
