/* -*- Mode: C; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*- */
/* vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 */
/*

   MDAnalysis --- https://www.mdanalysis.org
   Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
   (see the file AUTHORS for the full list of names)

   Released under the Lesser GNU Public Licence, v2.1 or any higher version

   Please cite your use of MDAnalysis in published work:

   R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
   D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
   MDAnalysis: A Python package for the rapid analysis of molecular dynamics
   simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
   Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.

   N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
   MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
   J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787

*/

#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "xtc_seek.h"
#include <stdio.h>
#include <stdlib.h>

enum { FALSE, TRUE };

int read_xtc_n_frames(char *fn, int *n_frames, int *est_nframes,
                      int64_t **offsets) {
  XDRFILE *xd;
  int framebytes, natoms, step;
  float time;
  int64_t filesize;

  if ((xd = xdrfile_open(fn, "r")) == NULL)
    return exdrFILENOTFOUND;

  if (xtc_header(xd, &natoms, &step, &time, TRUE) != exdrOK) {
    xdrfile_close(xd);
    return exdrHEADER;
  }

  if (xdr_seek(xd, 0L, SEEK_END) != exdrOK) {
    xdrfile_close(xd);
    return exdrNR;
  }
  filesize = xdr_tell(xd);

  /* Case of fewer than 10 atoms. Framesize known. */
  if (natoms < 10) {
    int i;
    xdrfile_close(xd);
    framebytes = XTC_SHORTHEADER_SIZE + XTC_SHORT_BYTESPERATOM * natoms;
    *n_frames = filesize / framebytes; /* Should we complain if framesize
                                          doesn't divide filesize? */
    /* Allocate memory for the frame index array */
    if ((*offsets = malloc(sizeof(int64_t) * (*n_frames))) == NULL)
      return exdrNOMEM;
    for (i = 0; i < *n_frames; i++) {
      (*offsets)[i] = i * framebytes;
    }
    *est_nframes = *n_frames;
    return exdrOK;
  } else /* No easy way out. We must iterate. */
  {
    /* Estimation of number of frames, with 20% allowance for error. */
    if (xdr_seek(xd, (int64_t)XTC_HEADER_SIZE, SEEK_SET) != exdrOK) {
      xdrfile_close(xd);
      return exdrNR;
    }
    if (xdrfile_read_int(&framebytes, 1, xd) == 0) {
      xdrfile_close(xd);
      return exdrENDOFFILE;
    }
    framebytes =
        (framebytes + 3) & ~0x03; // Rounding to the next 32-bit boundary
    *est_nframes =
        (int)(filesize / ((int64_t)(framebytes + XTC_HEADER_SIZE)) +
              1); // add one because it'd be easy to underestimate low
                  // frame numbers.
    *est_nframes += *est_nframes / 5;

    /* Allocate memory for the frame index array */
    if ((*offsets = malloc(sizeof(int64_t) * *est_nframes)) == NULL) {
      xdrfile_close(xd);
      return exdrNOMEM;
    }
    (*offsets)[0] = 0L;
    *n_frames = 1;
    while (1) {
      if (xdr_seek(xd, (int64_t)(framebytes + XTC_HEADER_SIZE), SEEK_CUR) !=
          exdrOK) {
        free(*offsets);
        xdrfile_close(xd);
        return exdrNR;
      }
      if (xdrfile_read_int(&framebytes, 1, xd) == 0)
        break;
      /* Read was successful; this is another frame */
      /* Check if we need to enlarge array */
      if (*n_frames == *est_nframes) {
        *est_nframes += *est_nframes / 5 + 1; // Increase in 20% stretches
        if ((*offsets = realloc(*offsets, sizeof(int64_t) * *est_nframes)) ==
            NULL) {
          free(*offsets);
          xdrfile_close(xd);
          return exdrNOMEM;
        }
      }
      (*offsets)[*n_frames] =
          xdr_tell(xd) - 4L - (int64_t)(XTC_HEADER_SIZE); // Account for the
                                                          // header and the
                                                          // nbytes bytes we
                                                          // read.
      (*n_frames)++;
      framebytes =
          (framebytes + 3) & ~0x03; // Rounding to the next 32-bit boundary
    }
    xdrfile_close(xd);
    return exdrOK;
  }
}
