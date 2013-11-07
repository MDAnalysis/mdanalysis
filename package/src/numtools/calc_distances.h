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

static void calc_bond_distance(coordinate* atom1, coordinate* atom2, int numatom, float* box, double* distances)
{
  int i;
  double dx[3];
  float inverse_box[3];
  double rsq;

  inverse_box[0] = 1.0/box[0];
  inverse_box[1] = 1.0/box[1];
  inverse_box[2] = 1.0/box[2];
 
  for (i=0; i<numatom; i++) {
    dx[0] = atom1[i][0] - atom2[i][0];
    dx[1] = atom1[i][1] - atom2[i][1];
    dx[2] = atom1[i][2] - atom2[i][2];
    // PBC time!
    minimum_image(dx, box, inverse_box);
    rsq = (dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]);
    *(distances+i) = sqrt(rsq);
  }
}

static void calc_bond_distance_noPBC(coordinate* atom1, coordinate* atom2, int numatom, double* distances)
{
  int i;
  double dx[3];
  double rsq;
 
  for (i=0; i<numatom; i++) {
    dx[0] = atom1[i][0] - atom2[i][0];
    dx[1] = atom1[i][1] - atom2[i][1];
    dx[2] = atom1[i][2] - atom2[i][2];
    rsq = (dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]);
    *(distances+i) = sqrt(rsq);
  }
}

static void calc_angle(coordinate* atom1, coordinate* atom2, coordinate* atom3, int numatom, double* angles)
{
  int i;
  double rij[3], rjk[3];
  double dotp, norm[2];

  for (i=0; i<numatom; i++) {
    rij[0] = atom1[i][0] - atom2[i][0];
    rij[1] = atom1[i][1] - atom2[i][1];
    rij[2] = atom1[i][2] - atom2[i][2];
    norm[0] = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);

    rjk[0] = atom3[i][0] - atom2[i][0];
    rjk[1] = atom3[i][1] - atom2[i][1];
    rjk[2] = atom3[i][2] - atom2[i][2];
    norm[1] = sqrt(rjk[0]*rjk[0] + rjk[1]*rjk[1] + rjk[2]*rjk[2]);

    dotp = rij[0] * rjk[0] + rij[1] * rjk[1] + rij[2] * rjk[2];

    *(angles+i) = acos(dotp/norm[0]/norm[1]);
  }
}

static void calc_torsion(coordinate* atom1, coordinate* atom2, coordinate* atom3, coordinate* atom4,
                         int numatom, double* angles)
{
  int i;
  double va[3], vb[3], vc[3];
  double rm[3], rn[3];
  double m, n, dotp;

  for (i=0; i<numatom; i++) {
    // connecting vectors between all 4 atoms: 1 -va-> 2 -vb-> 3 -vc-> 4
    va[0] = atom2[i][0] - atom1[i][0];
    va[1] = atom2[i][1] - atom1[i][1];
    va[2] = atom2[i][2] - atom1[i][2];

    vb[0] = atom3[i][0] - atom2[i][0];
    vb[1] = atom3[i][1] - atom2[i][1];
    vb[2] = atom3[i][2] - atom2[i][2];

    vc[0] = atom4[i][0] - atom3[i][0];
    vc[1] = atom4[i][1] - atom3[i][1];
    vc[2] = atom4[i][2] - atom3[i][2];

    // rm is normal vector to va vb
    // rn is normal vector is vb vc
    // m & n are norms to these vectors
    rm[0] = va[1]*vb[2] - va[2]*vb[1];
    rm[1] = va[0]*vb[2] - va[2]*vb[0];
    rm[2] = va[0]*vb[1] - va[1]*vb[0];
    m = sqrt(rm[0]*rm[0] + rm[1]*rm[1] + rm[2]*rm[2]);

    rn[0] = vb[1]*vc[2] - vb[2]*vc[1];
    rn[1] = vb[0]*vc[2] - vb[2]*vc[0];
    rn[2] = vb[0]*vc[1] - vb[1]*vc[0];
    n = sqrt(rn[0]*rn[0] + rn[1]*rn[1] + rn[2]*rn[2]);

    dotp = rm[0]*rn[0] + rm[1]*rn[1] + rm[2]*rn[2];

    *(angles + i) = acos(dotp/m/n);
  }

}
#endif
