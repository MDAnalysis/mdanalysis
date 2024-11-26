/* -*- tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
  vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

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
#include <float.h>

/* Helper functions */

inline int trmIndex(int row, int col) {  // array index for triangular matrix
	return (row) > (col) ? ((row)+1)*(row)/2+(col) : ((col)+1)*(col)/2+(row);
}

inline int sqmIndex(int colsn, int row, int col) { // array index for square matrix
	return row*colsn + col;
}

#ifndef _WIN32
float min(float * values, int length) { //array min
	float min = values[0];
	for (int i=1;i<length;i++) {
		if (values[i] < min) {
			min = values[i];
		}
	}
	return min;
}

float max(float * values, int length) { //array max
	float max = values[0];
	for (int i=1;i<length;i++) {
		if (values[i] > max) {
			max = values[i];
		}
	}
	return max;
}
#endif

void printarray(float* array, int lenarray) { //print an array, for debug purposes
	for (int i=0;i<lenarray;i++) {
		printf("%1.2f\n",array[i]);
	}
}

void printsqmatrix(float* array, int lenarray) { //print a square matrix, for debug purposes
	//printf("%s","\n");
	for (int i=0;i<lenarray;i++) {
		for (int j=0;j<lenarray;j++) {
			printf("%1.2f\t",array[sqmIndex(lenarray,i,j)]);
		}
		printf("%s","\n");
	}
}

void printtrmatrix(float* array, int lenarray) { //print a triangular matrix, for debug purposes
	//printf("%s","\n");
	for (int i=0; i<lenarray; i++) {
		for (int j=0; j<=i; j++) {
			//printf("%d\t",trmIndex(i,j));
			printf("%1.2f\t",array[trmIndex(i,j)]);
		}
	printf("%s","\n");
	}
}

int CAffinityPropagation(float *s, int n, float lambda, int max_iterations, int convergence, int noise, long* clusters) { // Affinity Propagation clustering algorithm

    /* n: number of elements
       s: similarity matrix
       lambda: damping parameter ([0.5;1.0[)
       max_iterations: maximum number of iterations
       convergence: convergence reached when centroids are stable for convergence iterations
       noise: apply noise to input similarities to eliminate redundancy */

	float *r =              (float *)  calloc( n*n , sizeof(float));  // N*N responsibilities matrix
	float *a  =             (float *)  calloc( n*n , sizeof(float));  // N*N availabilities matrix
    int *exemplars =          (int *)  malloc( n   * sizeof(int));        // N array of exemplars
    int *old_exemplars =      (int *)  malloc( n   * sizeof(int));        // N array of old exemplars, for convergence checking


	int i = 0;                        // index i over elements
	int k = 0;                        // index k over elements
	int j = 0;                        // generic index
    int idx = 0;                      // index for triangular matrix
    int sqm_idx = 0;                  // index for square matrix
    int currit = 0;                   // current iteration number
    int conv_count = 0;               // number of iterations with constant centroids so far
	float tmpsum = 0.0, maxsim = 0.0, this_tmpsum = 0.0;      // accumulators
	float tmp = 0.0;                  // temporary value
	float max1 = 0;
	float max2 = 0;
	int conv_reached = 0;        // convergence flag
	int has_cluster = 0;         // found clusters flag
	float lamprev = 1.0 - lambda;     // 1-lambda
	int n_clusters = 0; 			// number of clusters


    if (noise != 0) { // Add noise to data
        for (int i=0;i<n*(n+1)/2;i++) {
            s[i] = s[i] + (1e-16*s[i] )*(rand()/((float)RAND_MAX+1));
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
			max1 = -FLT_MAX;
			max2 = -FLT_MAX;
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
                if (exemplars[j] != old_exemplars[j]) {
                    conv_count = 0;
                    break;
                }
            }
        }
        else conv_count = 0;

        if (conv_count == convergence) conv_reached = 1; // check convergence

        currit++; // increment iteration number
    } // start a new iteration. If convergence or max_iterations reached

    //printf("Preference %3.2f: Convergence reached at iteration %d!\n",currit); // print convergence info

    // find number of clusters
	for (int i=0;i<n;i++) {
		idx = sqmIndex(n,i,i);
		if (r[idx]+a[idx] > 0) {
			exemplars[n_clusters] = i;
			n_clusters++;
		}
	}

    for (int i=0;i<n;i++) { // assign elements to clusters
        idx = sqmIndex(n,i,0);
        maxsim = r[idx]+a[idx];
        k=0;
	//printf("%3.1f, ",maxsim);
        for (k=1;k<n;k++) {
            idx = sqmIndex(n,i,k);
            tmpsum = r[idx]+a[idx];
            if (tmpsum >= maxsim) {
                clusters[i] = k;
                maxsim = tmpsum;
            }
        }
    }

	for (int i=0;i<n_clusters;i++) {
		clusters[exemplars[i]] = exemplars[i];
	}





    //Free memory anyway
    free(r);
    free(a);
    free(exemplars);
    free(old_exemplars);

    return conv_reached == 1 ? currit : -currit;
}
