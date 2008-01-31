
#ifndef __DISTANCES_H
#define __DISTANCES_H

#include <math.h>

typedef float coordinate[3];

static void minimum_image(double *x, float *box, float *box_half)
{
	if (fabs(x[0]) > box_half[0]) {
		if (x[0] < 0.0) x[0] += box[0];
		else x[0] -= box[0];
	}
	if (fabs(x[1]) > box_half[1]) {
		if (x[1] < 0.0) x[1] += box[1];
		else x[1] -= box[1];
	}
	if (fabs(x[2]) > box_half[2]) {
		if (x[2] < 0.0) x[2] += box[2];
		else x[2] -= box[2];
	}
}

static void calc_distance_array(coordinate* ref, int numref, coordinate* conf, int numconf, float* box, double* distances)
{
	int i, j;
	double dx[3];
	float box_half[3];
	double rsq;

	box_half[0] = box[0]/2;
	box_half[1] = box[1]/2;
	box_half[2] = box[2]/2;
	
	for (i=0; i < numref; i++) {
		for (j=0; j < numconf; j++) {
			dx[0] = conf[j][0]-ref[i][0];
			dx[1] = conf[j][1]-ref[i][1];
			dx[2] = conf[j][2]-ref[i][2];
			// Periodic boundaries
			minimum_image(dx, box, box_half);
			rsq = (dx[0]*dx[0])+(dx[1]*dx[1])+(dx[2]*dx[2]);
			*(distances+i*numconf+j) = sqrt(rsq);
		}
	}
}

#endif
