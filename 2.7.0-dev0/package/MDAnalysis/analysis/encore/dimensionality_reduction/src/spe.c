/* -*- tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
  vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

  MDAnalysis --- https://www.mdanalysis.org
  Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
  (see the file AUTHORS for the full list of names)

  Released under the GNU Public Licence, v2 or any higher version

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
#include <math.h>
#include <time.h>
#include <sys/types.h>

#ifdef _WIN32
    #include <io.h>
#else /* unix-like __unix__ || __APPLE__ */
    #include <unistd.h>
#endif

#define EPSILON 1e-8

typedef struct {
    int index;
    double value;
} IVWrapper;

inline int trmIndex(int row, int col) {  // array index for triangular matrix
	return (row) > (col) ? ((row)+1)*(row)/2+(col) : ((col)+1)*(col)/2+(row);
}

void printarray(double* arr, int len) {
    for (int i=0; i<len; i++) {
        printf("%.4f ",arr[i]); }
    printf("\n");
}

int cmp_ivwrapper(const void* ptra, const void* ptrb) {
    IVWrapper* a = (IVWrapper*) ptra;
    IVWrapper* b = (IVWrapper*) ptrb;
    if (a->value == b->value)
        return(0);
    else {
        if (a->value < b->value)
            return(-1);
        else
            return(1);
        }
}

// Euclidean distance
double ed(double* d_coords, int i, int j, int dim) {
    double d = 0.0;
    double delta = 0.0;
    int idxi = i*dim;
    int idxj = j*dim;
    int idxik = 0;
    int idxjk = 0;

    for (int k=0; k<dim; k++) {
        idxik = idxi+k;
        idxjk = idxj+k;
        delta = d_coords[idxik] - d_coords[idxjk];
        d = d + delta*delta;
    }
    return(sqrt(d));
}

double stress(double *s, double *d_coords, int dim, int elemsn) {
    double denom = 0.0;
    double numer = 0.0;
    double dab = 0.0;
    double delta = 0.0;
    int k = 0;

    for (int i=0; i<elemsn; i++) {
        for (int j=0; j<i; j++) {
            dab = ed(d_coords, i, j, dim);
            denom += s[k];
            delta = dab - s[k];
            numer += delta*delta / s[k];
            //printf("%d, dab %.3f, denom %.3f, numer %.3f, s[k]: %.3f\n", k, dab, denom, numer,s[k]);
            k++;
            //printf("a %d, b %d, dab %.3f, rab %.3f, delta %.3f\n",i,j,dab,s[k],dab-s[k]);
        }
        k++;
    }
    return( numer/denom );
}

int neighbours(double *s, int nelem, double cutoff, int* s_indeces, int* ioffsets, int* js) {
    int idx = 0;
    int global_counter = 0;
    int offset_counter = 0;
    ioffsets[0] = 0;
    for (int i=0;i<nelem;i++) {
        offset_counter = 0;
        for (int j=0;j<nelem;j++) {
            idx = trmIndex(i,j);
            if (s[idx] < cutoff) {
                s_indeces[global_counter] = idx; //XXX: replace numbers with addresses
                js[global_counter] = j;
                offset_counter++;
                global_counter++;
            }
        }
        ioffsets[i+1] = ioffsets[i] + offset_counter;
    }
    return global_counter;
}

int* nearest_neighbours(double *s, int nelem, int kn) {
    IVWrapper* ivpairs = (IVWrapper*) malloc((nelem-1)*sizeof(IVWrapper));
    int* neighbours = (int*) malloc(nelem*kn*sizeof(int));
    int totk = 0;
    int k = 0;
    for (int i=0;i<nelem;i++) {
        k = 0;
        for (int j=0;j<nelem;j++) {
            if (i!=j) {
                ivpairs[k].index = j;
                ivpairs[k].value = s[trmIndex(i,j)];
                k++;
            }
        }
        qsort(ivpairs, nelem-1, sizeof(IVWrapper), cmp_ivwrapper);
        for (int k=0;k<kn;k++) {
            neighbours[totk] = ivpairs[k].index;
            totk++;
        }
    }
    free(ivpairs);
    return neighbours;
}

double neighbours_stress(double *s, double *d_coords, int dim, int elemsn, double rco) {
    double denom = 0.0;
    double numer = 0.0;
    double dab = 0.0;
    double delta = 0.0;
    int k = 0;

    for (int i=0; i<elemsn; i++) {
        for (int j=0; j<i; j++) {
            dab = ed(d_coords, i, j, dim);
	    if (s[k] <= rco || dab < s[k]) {
	        denom += s[k];
                delta = dab - s[k];
                numer += delta*delta / s[k];
	    }
            //printf("%d, dab %.3f, denom %.3f, numer %.3f, s[k]: %.3f\n", k, dab, denom, numer,s[k]);
            k++;
            //printf("a %d, b %d, dab %.3f, rab %.3f, delta %.3f\n",i,j,dab,s[k],dab-s[k]);
        }
        k++;
    }
    return( numer/denom );
}


double CStochasticProximityEmbedding(
        double* s,
        double* d_coords,
        double rco,
        int nelem,
        int dim,
        double maxlam,
        double minlam,
        int ncycle,
        int nstep,
        int stressfreq) {

    int a = 0, b = 0, idxa = 0, idxb = 0, idxak = 0, idxbk = 0;
    double dab = 0.0, rab = 0.0;
    double lam = maxlam;
    double t = 0.0;
    double finalstress = 0.0;

    /* random init of d */
    srand(time(NULL)+getpid()*getpid());
    for (int i=0; i<nelem*dim; i++) {
        d_coords[i] = (double) rand() / (double) RAND_MAX;
    }

    /* start self organization */
    for (int i=0; i<ncycle; i++) {
        for (int j=0; j<nstep; j++) {

            a = rand() % nelem;
            while(1) {
                b = rand() % nelem;
                if (b != a) break;
            }

            dab = ed(d_coords, a, b, dim);
            rab = s[trmIndex(a, b)];

            if (rab <= rco || (rab > rco && dab < rab)) {
                idxa = a * dim;
                idxb = b * dim;
                t = lam * 0.5 * (rab - dab) / (dab + EPSILON);

                for (int k=0; k<dim; k++) {
                    idxak = idxa+k;
                    idxbk = idxb+k;
                    d_coords[idxak] = d_coords[idxak] + t*(d_coords[idxak] - d_coords[idxbk]);
                    d_coords[idxbk] = d_coords[idxbk] + t*(d_coords[idxbk] - d_coords[idxak]);
                }
            }
        }
        lam = lam - (maxlam - minlam) / (double)(ncycle - 1);
        //if (i % stressfreq == 0 && i != 0 && stressfreq > 0)
	    //printf("Cycle %d - Residual stress: %.3f, lambda %.3f\n", i, neighbours_stress(s, d_coords, dim, nelem, rco),lam);
    }
    finalstress = neighbours_stress(s, d_coords, dim, nelem, rco);
    //printf("Calculation finished. - Residual stress: %.3f\n", finalstress);
    return(finalstress);
}
