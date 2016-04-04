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
 * this
 * list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include <stdio.h>
#include <stdlib.h>

#define MAGIC 1995

enum { FALSE, TRUE };

int xtc_header(XDRFILE *xd, int *natoms, int *step, float *time,
                      mybool bRead) {
  int result, magic, n = 1;

  /* Note: read is same as write. He he he */
  magic = MAGIC;
  if ((result = xdrfile_write_int(&magic, n, xd)) != n) {
    if (bRead)
      return exdrENDOFFILE;
    else
      return exdrINT;
  }
  if (magic != MAGIC)
    return exdrMAGIC;
  if ((result = xdrfile_write_int(natoms, n, xd)) != n)
    return exdrINT;
  if ((result = xdrfile_write_int(step, n, xd)) != n)
    return exdrINT;
  if ((result = xdrfile_write_float(time, n, xd)) != n)
    return exdrFLOAT;

  return exdrOK;
}

static int xtc_coord(XDRFILE *xd, int *natoms, matrix box, rvec *x, float *prec,
                     mybool bRead) {
  int result;

  /* box */
  result = xdrfile_read_float(box[0], DIM * DIM, xd);
  if (DIM * DIM != result)
    return exdrFLOAT;
  else {
    if (bRead) {
      result = xdrfile_decompress_coord_float(x[0], natoms, prec, xd);
      if (result != *natoms)
        return exdr3DX;
    } else {
      result = xdrfile_compress_coord_float(x[0], *natoms, *prec, xd);
      if (result != *natoms)
        return exdr3DX;
    }
  }
  return exdrOK;
}

int read_xtc_natoms(char *fn, int *natoms) {
  XDRFILE *xd;
  int step, result;
  float time;

  xd = xdrfile_open(fn, "r");
  if (NULL == xd)
    return exdrFILENOTFOUND;
  result = xtc_header(xd, natoms, &step, &time, TRUE);
  xdrfile_close(xd);

  return result;
}

int read_xtc(XDRFILE *xd, int natoms, int *step, float *time, matrix box,
             rvec *x, float *prec)
/* Read subsequent frames */
{
  int result;

  if ((result = xtc_header(xd, &natoms, step, time, TRUE)) != exdrOK)
    return result;

  if ((result = xtc_coord(xd, &natoms, box, x, prec, 1)) != exdrOK)
    return result;

  return exdrOK;
}


int write_xtc(XDRFILE *xd, int natoms, int step, float time, matrix box,
              rvec *x, float prec)
/* Write a frame to xtc file */
{
  int result;

  if ((result = xtc_header(xd, &natoms, &step, &time, FALSE)) != exdrOK)
    return result;

  if ((result = xtc_coord(xd, &natoms, box, x, &prec, 0)) != exdrOK)
    return result;

  return exdrOK;
}
