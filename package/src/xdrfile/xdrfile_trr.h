/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- 
 *
 * $Id$
 *
 * Copyright (c) Erik Lindahl, David van der Spoel 2003,2004.
 * Copyright (c) Manuel Melo <manuel.nuno.melo@gmail.com> 2013,2014.
 * Coordinate compression (c) by Frans van Hoesel. 
 * XTC/TRR seeking and indexing (c) Manuel Melo.
 *
 *    This file is part of libxdrfile2.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
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
   
  /* This function returns the number of atoms in the xtc file in *natoms */
  extern int read_trr_natoms(char *fn,int *natoms);

  /* Skip through trajectory, reading headers, obtain the total number of frames in the trr */ 
  extern int read_trr_n_frames(char *fn, int *n_frames, int64_t **offsets);

  /* Read one frame of an open trr file. If either of x,v,f,box are
     NULL the arrays will be read from the file but not used.  */
  extern int read_trr(XDRFILE *xd,int natoms,int *step,float *t,float *lambda,
		      matrix box,rvec *x,rvec *v,rvec *f, int *has_prop);

  /* Write a frame to trr file */
  extern int write_trr(XDRFILE *xd,int natoms,int step,float t,float lambda,
		       matrix box,rvec *x,rvec *v,rvec *f);

/* Minimum TRR header size. It can have 8 bytes more if we have double time and lambda. */
#define TRR_MIN_HEADER_SIZE 54
#define TRR_DOUBLE_XTRA_HEADER 8
  
/* Flags to signal the update of pos/vel/forces */
#define HASX 1
#define HASV 2
#define HASF 4

#ifdef __cplusplus
}
#endif

#endif
