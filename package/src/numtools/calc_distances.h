/* -*- Mode: C; tab-width: 4; indent-tabs-mode:nil; -*- */
/* vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 */
/*
  MDAnalysis --- http://mdanalysis.googlecode.com
  Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
                Elizabeth J. Denning, Oliver Beckstein,
                and contributors (see website for details)
  Released under the GNU Public Licence, v2 or any higher version

  Please cite your use of MDAnalysis in published work:

      N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
      O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
      Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
      in press.
*/

#ifndef __DISTANCES_H
#define __DISTANCES_H

#include <math.h>

#include <float.h>
typedef float coordinate[3];

static void minimum_image(double *x, float *box, float *inverse_box)
{
    int i;
    double s;
    for (i = 0; i < 3; i++) {
        if (box[i] > FLT_EPSILON) {
            s = inverse_box[i] * x[i];
            x[i] = box[i] * (s - round(s));
        }
    }
}
static void minimum_image_triclinic(double *dx, coordinate* box, float* box_half)
{
  // Minimum image convention for triclinic systems, modelled after domain.cpp in LAMMPS
  // Assumes that there is a maximum separation of 1 box length (enforced in dist functions
  // by moving all particles to inside the box before calculating separations)
  // Requires a box 
  // Assumes box having zero values for box[0][1], box[0][2] and box [1][2]

  // z
  if (fabs(dx[2]) > box_half[2]) {
    if (dx[2] < 0.0 ) {
      dx[2] += box[2][2];
      dx[1] += box[2][1];
      dx[0] += box[2][0];
    } else {
      dx[2] -= box[2][2];
      dx[1] -= box[2][1];
      dx[0] -= box[2][0];
    }
  }
  // y
  if (fabs(dx[1]) > box_half[1]) {
    if (dx[1] < 0.0) {
      dx[1] += box[1][1];
      dx[0] += box[1][0];
    } else {
      dx[1] -= box[1][1];
      dx[0] -= box[1][0];
    }
  }
  // x
  if ( fabs(dx[0]) > box_half[0]) {
    if (dx[0] < 0.0) {
      dx[0] += box[0][0];
    } else {
      dx[0] -= box[0][0];
    }
  }
}
static void calc_distance_array(coordinate* ref, int numref, coordinate* conf, int numconf, float* box, double* distances)
{
	int i, j;
	double dx[3];
	float inverse_box[3];
	double rsq;

	inverse_box[0] = 1.0/box[0];
	inverse_box[1] = 1.0/box[1];
	inverse_box[2] = 1.0/box[2];
	
	for (i=0; i < numref; i++) {
		for (j=0; j < numconf; j++) {
            dx[0] = conf[j][0]-ref[i][0];
			dx[1] = conf[j][1]-ref[i][1];
			dx[2] = conf[j][2]-ref[i][2];
			// Periodic boundaries
			minimum_image(dx, box, inverse_box);
			rsq = (dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]);
			*(distances+i*numconf+j) = sqrt(rsq);
		}
	}
}

static void calc_distance_array_noPBC(coordinate* ref, int numref, coordinate* conf, int numconf, double* distances)
{
	int i, j;
	double dx[3];
	double rsq;

	for (i=0; i < numref; i++) {
		for (j=0; j < numconf; j++) {
			dx[0] = conf[j][0]-ref[i][0];
			dx[1] = conf[j][1]-ref[i][1];
			dx[2] = conf[j][2]-ref[i][2];
			rsq = (dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]);
			*(distances+i*numconf+j) = sqrt(rsq);
		}
	}
}


static void calc_self_distance_array(coordinate* ref, int numref, float* box, double* distances, int distnum)
{
	int i, j, distpos;
	double dx[3];
	float inverse_box[3];
	double rsq;

	inverse_box[0] = 1.0/box[0];
	inverse_box[1] = 1.0/box[1];
	inverse_box[2] = 1.0/box[2];
	
	distpos = 0;
	for (i=0; i < numref; i++) {
		for (j=i+1; j < numref; j++) {
			dx[0] = ref[j][0]-ref[i][0];
			dx[1] = ref[j][1]-ref[i][1];
			dx[2] = ref[j][2]-ref[i][2];
			// Periodic boundaries
			minimum_image(dx, box, inverse_box);
			rsq = (dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]);
			*(distances+distpos) = sqrt(rsq);
			distpos += 1;
		}
	}
}


static void calc_self_distance_array_noPBC(coordinate* ref, int numref, double* distances, int distnum)
{
	int i, j, distpos;
	double dx[3];
	double rsq;

	distpos = 0;
	for (i=0; i < numref; i++) {
		for (j=i+1; j < numref; j++) {
			dx[0] = ref[j][0]-ref[i][0];
			dx[1] = ref[j][1]-ref[i][1];
			dx[2] = ref[j][2]-ref[i][2];
			rsq = (dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]);
			*(distances+distpos) = sqrt(rsq);
			distpos += 1;
		}
	}
}

static void coord_transform(coordinate* coords, int numCoords, coordinate* box)
{
  int i, j, k;
  float new[3];
  // Matrix multiplication inCoords * box = outCoords  
  // Multiplication done in place using temp array 'new'
  // Used to transform coordinates to/from S/R space in trilinic boxes
  for (i=0; i < numCoords; i++){
    new[0] = 0.0;
    new[1] = 0.0;
    new[2] = 0.0;    
    for (j=0; j < 3; j++){
      for (k=0; k < 3; k++){
        new[j] += coords[i][k] * box[k][j];
      }
    }
    coords[i][0] = new[0];
    coords[i][1] = new[1];
    coords[i][2] = new[2];
  }
}
static void triclinic_pbc(coordinate* coords, int numcoords, coordinate* box, float* box_inverse){
  int i, s;
  // Moves all coordinates to within the box boundaries for a triclinic box
  // Assumes box having zero values for box[0][1], box[0][2] and box [1][2]
  for (i=0; i < numcoords; i++){
    // z
    s = floor(coords[i][2] * box_inverse[2]);
    coords[i][2] -= s * box[2][2];
    coords[i][1] -= s * box[2][1];
    coords[i][0] -= s * box[2][0];
    // y
    s = floor(coords[i][1] * box_inverse[1]);
    coords[i][1] -= s * box[1][1];
    coords[i][0] -= s * box[1][0];
    // x
    s = floor(coords[i][0] * box_inverse[0]);
    coords[i][0] -= s * box[0][0];
  }
}
static void calc_distance_array_triclinic(coordinate* ref, int numref, coordinate* conf, int numconf, coordinate* box, double* distances)
{
  int i, j;
  double dx[3];
  float box_half[3], box_inverse[3];
  double rsq;
  
  box_half[0] = 0.5 * box[0][0];
  box_half[1] = 0.5 * box[1][1];
  box_half[2] = 0.5 * box[2][2];

  box_inverse[0] = 1.0 / box[0][0];
  box_inverse[1] = 1.0 / box[1][1];
  box_inverse[2] = 1.0 / box[2][2];
  // Move coords to inside box
  triclinic_pbc(ref, numref, box, box_inverse);
  triclinic_pbc(conf, numconf, box, box_inverse);

  for (i=0; i < numref; i++){
    for (j=0; j < numconf; j++){
      dx[0] = conf[j][0] - ref[i][0];
      dx[1] = conf[j][1] - ref[i][1];
      dx[2] = conf[j][2] - ref[i][2];
      minimum_image_triclinic(dx, box, box_half);
      rsq = (dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      *(distances + i*numconf + j) = sqrt(rsq);
    }
  }
}
static void calc_self_distance_array_triclinic(coordinate* ref, int numref, coordinate* box, double *distances, int distnum)
{
  int i, j, distpos;
  double dx[3];
  double rsq;
  float box_half[3], box_inverse[3];

  box_half[0] = 0.5 * box[0][0];
  box_half[1] = 0.5 * box[1][1];
  box_half[2] = 0.5 * box[2][2];

  box_inverse[0] = 1.0 / box[0][0];
  box_inverse[1] = 1.0 / box[1][1];
  box_inverse[2] = 1.0 / box[2][2];

  triclinic_pbc(ref, numref, box, box_inverse);

  distpos = 0;
  for (i=0; i < numref; i++){
    for (j=i+1; j < numref; j++){
      dx[0] = ref[j][0] - ref[i][0];
      dx[1] = ref[j][1] - ref[i][1];
      dx[2] = ref[j][2] - ref[i][2];
      minimum_image_triclinic(dx, box, box_half);
      rsq = (dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      *(distances + distpos) = sqrt(rsq);
      distpos += 1;
    }
  }
}

#endif
