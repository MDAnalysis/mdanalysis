/*
ap.c --- C implementation of the Affinity Propagation clustering algorithm
Copyright (C) 2014 Wouter Boomsma, Matteo Tiberti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* Helper functions */

inline int trmIndex(int row, int col) {  // array index for triangular matrix
	return (row) > (col) ? ((row)+1)*(row)/2+(col) : ((col)+1)*(col)/2+(row);
}

inline int sqmIndex(int colsn, int row, int col) { // array index for square matrix
	return row*colsn + col;
}

inline double pwmax(double x, double y) { //pairwise min
	return x > y ? x : y;
}

inline double pwmin(double x, double y) { //pairwise max
	return x < y ? x : y;
}

double min(double * values, int length) { //array min
	double min = values[0];
	for (int i=1;i<length;i++) {
		if (values[i] < min) {
			min = values[i];
		}
	}
	return min;
}

double max(double * values, int length) { //array max
	double max = values[0];
	for (int i=1;i<length;i++) {
		if (values[i] > max) {
			max = values[i];
		}
	}
	return max;
}

void printarray(double* array, int lenarray) { //print an array, for debug purposes
	for (int i=0;i<lenarray;i++) {
		printf("%1.2f\n",array[i]);
	}
}

void printsqmatrix(double* array, int lenarray) { //print a square matrix, for debug purposes
	//printf("%s","\n");
	for (int i=0;i<lenarray;i++) {
		for (int j=0;j<lenarray;j++) {
			printf("%1.2f\t",array[sqmIndex(lenarray,i,j)]);
		}
		printf("%s","\n");
	}
}

void printtrmatrix(double* array, int lenarray) { //print a triangular matrix, for debug purposes
	//printf("%s","\n");
	for (int i=0; i<lenarray; i++) {
		for (int j=0; j<=i; j++) {
			//printf("%d\t",trmIndex(i,j));
			printf("%1.2f\t",array[trmIndex(i,j)]);
		}
	printf("%s","\n");
	}
}

int CAffinityPropagation(double *s, int n, double lambda, int max_iterations, int convergence, int noise, long* clusters) { // Affinity Propagation clustering algorithm

    /* n: number of elements
       s: similarity matrix
       lambda: damping parameter ([0.5;1.0[)
       max_iterations: maximum number of iterations
       convergence: convergence reached when centroids are stable for convergence iterations
       noise: apply noise to input similarities to eliminate redundancy */

	double *r =              (double *)  calloc( n*n , sizeof(double));  // N*N responsibilities matrix
	double *a  =             (double *)  calloc( n*n , sizeof(double));  // N*N availabilities matrix
    int *exemplars =          (int *)  malloc( n   * sizeof(int));        // N array of exemplars
    int *old_exemplars =      (int *)  malloc( n   * sizeof(int));        // N array of old exemplars, for convergence checking


	int i = 0;                        // index i over elements
	int k = 0;                        // index k over elements
	int j = 0;                        // generic index
    int idx = 0;                      // index for triangular matrix
    int sqm_idx = 0;                  // index for square matrix
    int currit = 0;                   // current iteration number
    int conv_count = 0;               // number of iterations with constant centroids so far
	double tmpsum = 0.0, maxsim = 0.0, this_tmpsum = 0.0;      // accumulators
	double tmp = 0.0;                  // temporary value
	double max1 = 0;
	double max2 = 0;
	int conv_reached = 0;        // convergence flag
	int has_cluster = 0;         // found clusters flag
	double lamprev = 1.0 - lambda;     // 1-lambda

    if (noise != 0) { // Add noise to data
        for (int i=0;i<n*(n+1)/2;i++) {
            s[i] = s[i] + (1e-16*s[i] )*(rand()/((double)RAND_MAX+1));
        }
     }

    for (int i=0;i<n*n;i++) {
        r[i] = 0.0;
        a[i] = 0.0;
    }

    for (i=0;i<n;i++) { // Initialize exemplars
		exemplars[i] = -1;
	}

    //printtrmatrix(s,n);
    //printf("\n");

    //printf("Preference %3.2f: running Affinity Propagation.\n",s[0]);    

	while (currit < max_iterations && conv_reached == 0) { // Start iterations

	// Update r

		for (int i=0;i<n;i++) {
			max1 = -DBL_MAX;
			max2 = -DBL_MAX;
			for (int j=0;j<n;j++) {
				sqm_idx = sqmIndex(n,i,j);
				idx = trmIndex(i,j);
				tmp = s[idx]+a[sqm_idx];
				if (tmp > max1) {
					max2 = max1;
					max1 = tmp;
				}
				else if (tmp > max2) {
					max2 = tmp;
				}
			}
			for (int k=0;k<n;k++) {
				idx = trmIndex(i,k);
				sqm_idx = sqmIndex(n,i,k);
				if (a[sqm_idx]+s[idx] == max1)
					r[sqm_idx] = lambda*r[sqm_idx] + lamprev*(s[idx] - max2);
				else
					r[sqm_idx] = lambda*r[sqm_idx] + lamprev*(s[idx] - max1);
			}
		}

		//printf("%s","r\n");
		//printsqmatrix(r, n);

	// Update a

		for (int k=0;k<n;k++) {
			tmpsum = 0.0;
			for (int j=0;j<n;j++) { //sum all the elements > 0 of column k
				tmp = r[sqmIndex(n,j,k)];
				if (j!=k) {
                    if (tmp > 0.0)
                        tmpsum = tmpsum + tmp; // if j!=k (r(j,k)): if > 0, sum it
                }
                else
                    tmpsum = tmpsum + tmp; // if j == k (r(k,k)): always sum it.
                    // we will have to remove the j==i (r(i,k)) case, for every i.
            }
            //printf("tempsum %1.2f\n",tmpsum);
			for (int i=0;i<n;i++) {
			    this_tmpsum = tmpsum;
				sqm_idx = sqmIndex(n,i,k);
				if (i != k) {
                    tmp = r[sqm_idx];
                    //printf("tmp %1.2f\n",tmp);
					if (tmp > 0.0)
						this_tmpsum = this_tmpsum - tmp; //subtract r(i,k)
                        //printf("tmpsum2 %1.2f\n",tmpsum);
					if (this_tmpsum < 0.0)
						a[sqm_idx] = lambda*a[sqm_idx] + lamprev*this_tmpsum;
					else
						a[sqm_idx] = lambda*a[sqm_idx];
				}
				else  // we don't need to remove the r(i,k) case, BUT wee need to remove to remove the r(k,k) case
					a[sqm_idx] = lambda*a[sqm_idx] + lamprev*(this_tmpsum - r[sqmIndex(n,k,k)]);
			}
		}

		//printf("%s","a\n");
		//printsqmatrix(a,n);

    //Check for convergence

        int* tmp_exemplars = old_exemplars; // current exemplars become old, before calculating the new ones. Pointer swap - fast and convenient
        old_exemplars = exemplars;
        exemplars = tmp_exemplars;

        has_cluster = 0;
        for (int i=0;i<n;i++) { // identify exemplars
            idx = sqmIndex(n,i,i);
            if (r[idx] + a[idx] > 0.0) {
                exemplars[i] = 1;
                has_cluster = 1;
            }
            else
                exemplars[i] = 0;
        }

        if (has_cluster != 0) {
            conv_count++;
            for (j=0;j<n;j++) { // check if exemplars have changed. If they have changed, or if no clusters have been identified, reset convergence counter.
                if (! exemplars[j] == old_exemplars[j]) {
                    conv_count = 0;
                    break;
                }
            }
        }
        else conv_count = 0;

        if (conv_count == convergence) conv_reached = 1; // check convergence

        currit++; // increment iteration number
    } // start a new iteration. If convergence or max_iterations reached

    if ( conv_reached == 1 ) {
        //printf("Preference %3.2f: Convergence reached at iteration %d!\n",currit); // print convergence info
        for (int i=0;i<n;i++) { // assign elements to clusters
            idx = sqmIndex(n,i,0);
            maxsim = r[idx]+a[idx];
		//printf("%3.1f, ",maxsim);
            for (k=1;k<n;k++) {
                idx = sqmIndex(n,i,k);
                tmpsum = r[idx]+a[idx];
			//Zprintf("%3.1f, ",tmpsum);
                if (tmpsum > maxsim) {
                    clusters[i] = k;
                    maxsim = tmpsum;
                }
            }
        }
    }
    else {
        for (int i=0;i<n;i++)
            clusters[i] = -1.0;
        //printf("\nPreference %3.2f: Convergence not reached in %d iterations.\n", s[0], currit+1);
    }
    //for (int i=0;i<n;i++) { if (exemplars[i] == 1) printf("%d\n",i); }

    //Free memory anyway
    free(r);
    free(a);
    free(exemplars);
    free(old_exemplars);

    return conv_reached == 1 ? currit : -currit;
}

