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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xdrfile.h"
#include "xdrfile_trr.h"
#include "trr_seek.h"

int read_trr_n_frames(char *fn, int *n_frames, int *est_nframes,
                      int64_t **offsets) {
  XDRFILE *xd;
  t_trnheader sh;
  float time, lambda;
  int result, framebytes, totalframebytes;
  int64_t filesize, frame_offset;

  if ((xd = xdrfile_open(fn, "r")) == NULL)
    return exdrFILENOTFOUND;
  if (xdr_seek(xd, 0L, SEEK_END) != exdrOK) {
    xdrfile_close(xd);
    return exdrNR;
  }
  filesize = xdr_tell(xd);
  if (xdr_seek(xd, 0L, SEEK_SET) != exdrOK) {
    xdrfile_close(xd);
    return exdrNR;
  }

  if ((result = do_trnheader(xd, 1, &sh)) != exdrOK) {
    xdrfile_close(xd);
    return result;
  }

  framebytes = sh.ir_size + sh.e_size + sh.box_size + sh.vir_size +
               sh.pres_size + sh.top_size + sh.sym_size + sh.x_size +
               sh.v_size + sh.f_size;

  *est_nframes =
      (int)(filesize / ((int64_t)(framebytes + TRR_MIN_HEADER_SIZE)) +
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
    if (xdr_seek(xd, (int64_t)(framebytes), SEEK_CUR) != exdrOK) {
      free(*offsets);
      xdrfile_close(xd);
      return exdrNR;
    }
    frame_offset = xdr_tell(xd); /* Store it now, before we read the header */
    if ((result = do_trnheader(xd, 1, &sh)) != exdrOK) /* Interpreting as EOF */
      break;
    /* Read was successful; this is another frame */
    /* Check if we need to enlarge array */
    if (*n_frames == *est_nframes) {
      *est_nframes += *est_nframes / 5 + 1; // Increase in 20% stretches
      if ((*offsets = realloc(*offsets, sizeof(int64_t) * *est_nframes)) ==
          NULL) {
        xdrfile_close(xd);
        return exdrNOMEM;
      }
    }
    (*offsets)[*n_frames] = frame_offset;
    (*n_frames)++;
    /* Calculate how much to skip this time */
    framebytes = sh.ir_size + sh.e_size + sh.box_size + sh.vir_size +
                 sh.pres_size + sh.top_size + sh.sym_size + sh.x_size +
                 sh.v_size + sh.f_size;
  }
  xdrfile_close(xd);
  return exdrOK;
}
