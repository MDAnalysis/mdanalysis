/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- 
 *
 * $Id$
 *
 * Copyright (c) Erik Lindahl, David van der Spoel 2003,2004.
 * Coordinate compression (c) by Frans van Hoesel. 
 * Copyright (c) Manuel Melo 2013,2014.
 * XTC/TRR seeking and indexing (c) Manuel Melo.
 *
 *    This file is part of libxdrfile2.0.
 * 
 *  libxdrfile2.0 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  libxdrfile2.0 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libxdrfile2.0.  If not, see <http://www.gnu.org/licenses/>.
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
  extern int read_xtc_natoms(char *fn,int *natoms);

  /* Seek through trajectory counting and indexing frames */
  extern int read_xtc_numframes(char *fn, int *numframes, int64_t **offsets);
  
  /* Read one frame of an open xtc file */
  extern int read_xtc(XDRFILE *xd,int natoms,int *step,float *time,
		      matrix box,rvec *x,float *prec);

  /* Write a frame to xtc file */
  extern int write_xtc(XDRFILE *xd,
		       int natoms,int step,float time,
		       matrix box,rvec *x,float prec);
  
/* XTC header fields until coord floats: *** only for trajectories of less than 10 atoms! ***  */
/* magic natoms step time DIM*DIM_box_vecs natoms */
#define XTC_SHORTHEADER_SIZE (20 + DIM*DIM*4)
/* Short XTCs store each coordinate as a 32-bit float. */
#define XTC_SHORT_BYTESPERATOM 12
/* XTC header fields until frame bytes: *** only for trajectories of more than 9 atoms! ***  */
/* magic natoms step time DIM*DIM_box_vecs natoms prec DIM_min_xyz DIM_max_xyz smallidx */
#define XTC_HEADER_SIZE (DIM*DIM*4 + DIM*2 + 46)

#ifdef __cplusplus
}
#endif

#endif
