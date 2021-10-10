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
#ifndef _xdrfile_trr_h
#define _xdrfile_trr_h

#ifdef __cplusplus
extern "C" {
#endif

#include "xdrfile.h"

/* All functions return exdrOK if succesfull.
 * (error codes defined in xdrfile.h).
 */

typedef struct /* This struct describes the order and the	*/
    /* sizes of the structs in a trjfile, sizes are given in bytes.	*/
    {
  mybool bDouble; /* Double precision?                            */
  int ir_size;    /* Backward compatibility		        */
  int e_size;     /* Backward compatibility		        */
  int box_size;   /* Non zero if a box is present			*/
  int vir_size;   /* Backward compatibility		        */
  int pres_size;  /* Backward compatibility		        */
  int top_size;   /* Backward compatibility		        */
  int sym_size;   /* Backward compatibility		        */
  int x_size;     /* Non zero if coordinates are present		*/
  int v_size;     /* Non zero if velocities are present		*/
  int f_size;     /* Non zero if forces are present		*/

  int natoms;     /* The total number of atoms			*/
  int step;       /* Current step number				*/
  int nre;        /* Backward compatibility		        */
  float tf;       /* Current time					*/
  float lambdaf;  /* Current value of lambda			*/
  double td;      /* Current time					*/
  double lambdad; /* Current value of lambda			*/
} t_trnheader;

/* This function returns the number of atoms in the trr file in *natoms */
extern int read_trr_natoms(char *fn, int *natoms);

/* Read one frame of an open trr file. If either of x,v,f,box are
   NULL the arrays will be read from the file but not used.  */
extern int read_trr(XDRFILE *xd, int natoms, int *step, float *t, float *lambda,
                    matrix box, rvec *x, rvec *v, rvec *f, int *has_prop);

/* Write a frame to trr file */
extern int write_trr(XDRFILE *xd, int natoms, int step, float t, float lambda,
                     matrix box, rvec *x, rvec *v, rvec *f);

/* unpack header of current frame */
extern int do_trnheader(XDRFILE *xd, mybool bRead, t_trnheader *sh);

/* Flags to signal the update of pos/vel/forces */
#define HASX 1
#define HASV 2
#define HASF 4

#ifdef __cplusplus
}
#endif

#endif
