/*
spe.c --- C implementation of the Stochastic Proximity Embedding algorithm and variations
Copyright (C) 2014 Wouter Boomsma, Matteo Tiberti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <time.h>

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




double CkNNStochasticProximityEmbedding(
        double* s,
        double* d_coords,
        int kn,
        int nelem,
        int dim,
        double maxlam,
        double minlam,
        int ncycle,
        int nstep,
        int stressfreq) {

    int* tmp;
    int* neighbours;
    int a = 0, b = 0, idx = 0, idxa = 0, idxb = 0, idxak = 0, idxbk = 0, idxab = 0;
    double dab = 0.0, rab = 0.0;
    double lam = maxlam;
    double t = 0.0;
    double finalstress = 0.0;

    srand(time(NULL)+getpid()*getpid());
    //time_t tempo;
    //time(&tempo);
    //printf("%s: Finding neighbours...\n",ctime(&tempo));
    neighbours = nearest_neighbours(s, nelem, kn);
    //time(&tempo);
    //printf("%s: Done!\n", ctime(&tempo));
    for (int i=0; i<nelem*dim; i++) {
        d_coords[i] = (double)rand();
    }

    for (int i=0; i<ncycle; i++) {
        for (int j=0; j<nstep; j++) {
            a = rand() % nelem;

            for (int k=kn; k<(a+1)*kn; k++) {
                b = neighbours[k];

                dab = ed(d_coords, a, b, dim);
                rab = s[trmIndex(a,b)];

                if (dab > rab) {
                    idxa = a * dim;
                    idxb = b * dim;
                    t = lam * 0.5 * (rab - dab) / (dab + EPSILON);

                    for (int k = 0; k < dim; k++) {
                        idxak = idxa+k;
                        idxbk = idxb+k;
                        d_coords[idxak] = d_coords[idxak] + t*(d_coords[idxak] - d_coords[idxbk]);
                        d_coords[idxbk] = d_coords[idxbk] + t*(d_coords[idxbk] - d_coords[idxak]);
                    }
                }
            }
        }
        if (i % stressfreq == 0 && i != 0 && stressfreq > 0)
            printf("Cycle %d - Residual stress: %.3f, lambda %.3f\n", i, stress(s, d_coords, dim, nelem),lam);
        lam = lam - (maxlam - minlam) / (double)(ncycle - 1);
    }
    free(neighbours);
    return(stress(s, d_coords, dim, nelem));
}

double CkNeighboursStochasticProximityEmbedding(
        double* s,
        double* d_coords,
        double rco,
        int kn,
        int nelem,
        int dim,
        double maxlam,
        double minlam,
        int ncycle,
        int stressfreq) {

    int* tmp;
    int a = 0, b = 0, idx = 0, idxa = 0, idxb = 0, idxak = 0, idxbk = 0, idxab = 0;
    double dab = 0.0, rab = 0.0;
    double lam = maxlam;
    double t = 0.0;
    double finalstress = 0.0;
    int* s_indeces = (int*) malloc(nelem*nelem*sizeof(int));
    int* ioffsets  = (int*) malloc(nelem      *sizeof(int));
    int* js        = (int*) malloc(nelem*nelem*sizeof(int));
    int nlistlen = 0;

    srand(time(NULL)+getpid()*getpid());
    nlistlen = neighbours(s, nelem, rco, s_indeces, ioffsets, js);

    s_indeces = (int*) realloc(s_indeces, nlistlen*sizeof(int));
    ioffsets  = (int*) realloc(ioffsets, nelem*sizeof(int));
    js        = (int*) realloc(js, nlistlen*sizeof(int));

    for (int i=0; i<nelem*dim; i++) {
        d_coords[i] = (double) rand() / RAND_MAX;
    }
    for (int i=0; i<nelem+1; i++) {
        //printf("ioff %d\n",ioffsets[i]);
    }

    /* start self organization */
    for (int i=0; i<ncycle; i++) {
        for (int a=0; a<nelem; a++) {
            for (int k=0; k<kn; k++) {

                int test = ioffsets[1] - ioffsets[0];

                b = rand() % (ioffsets[a+1] - ioffsets[a]);
                idxab = s_indeces[ioffsets[a]+b];
                b = js[ioffsets[a]+b];

                dab = ed(d_coords, a, b, dim);
                rab = s[idxab];

                if (dab < rab) {
                    idxa = a * dim;
                    idxb = b * dim;
                    t = lam * 0.5 * (rab - dab) / (dab + EPSILON);

                    for (int k = 0; k < dim; k++) {
                        idxak = idxa+k;
                        idxbk = idxb+k;
                        d_coords[idxak] = d_coords[idxak] + t*(d_coords[idxak] - d_coords[idxbk]);
                        d_coords[idxbk] = d_coords[idxbk] + t*(d_coords[idxbk] - d_coords[idxak]);
                    }
                }
            }
        }
        lam = lam - (maxlam - minlam) / (double)(ncycle - 1);
        if (i % stressfreq == 0 && i != 0 && stressfreq > 0)
            printf("Cycle %d - Residual stress: %.3f, lambda %.3f\n", i, stress(s, d_coords, dim, nelem),lam);
    }
    finalstress = stress(s, d_coords, dim, nelem);
    printf("Calculation finished (%d dimensions). - Residual stress: %.3f\n", dim, finalstress);
    return(finalstress);
    /* cleanup */
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

    int a = 0, b = 0, idx = 0, idxa = 0, idxb = 0, idxak = 0, idxbk = 0;
    double dab = 0.0, rab = 0.0;
    double lam = maxlam;
    double t = 0.0;
    double finalstress = 0.0;

    srand(time(NULL)+getpid()*getpid());
    /* random init of d */


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
        if (i % stressfreq == 0 && i != 0 && stressfreq > 0)
	  printf("Cycle %d - Residual stress: %.3f, lambda %.3f\n", i, neighbours_stress(s, d_coords, dim, nelem, rco),lam);
    }
    finalstress = neighbours_stress(s, d_coords, dim, nelem, rco);
    printf("Calculation finished. - Residual stress: %.3f\n", finalstress);
    return(finalstress);
}

// int main() {
//     double indata[72] = {-0.022,-3.647,2.514,0.519,-2.968,1.340,2.029,-2.951,1.374,2.666,-3.478,2.283,-0.063,-1.527,1.282,0.188,-0.757,-0.045,-0.388,0.668,-0.004,-0.136,1.349,-1.308,-0.497,2.593,-1.608,-1.121,3.383,-0.780,-0.215,3.048,-2.788,0.584,-4.045,3.240,0.211,-3.538,0.443,0.336,-0.947,2.140,-1.156,-1.569,1.465,-0.252,-1.330,-0.887,1.276,-0.706,-0.255,0.080,1.230,0.834,-1.476,0.622,0.217,0.348,0.865,-2.071,-1.312,2.964,0.133,-1.365,4.326,-1.087,0.275,2.415,-3.423,-0.502,4.005,-3.000};
//     //double indata[72] = {9.9780,6.3530,12.5140, 10.5190,7.0320,11.3400, 12.0290,7.0490,11.3740, 12.6660,6.5220,12.2830, 9.9370,8.4730,11.2820, 10.1880,9.2430,9.9550, 9.6120,10.6680,9.9960, 9.8640,11.3490,8.6920, 9.5030,12.5930,8.3920, 8.8790,13.3830,9.2200, 9.7850,13.0480,7.2120, 10.5840,5.9550,13.2400, 10.2110,6.4620,10.4430, 10.3360,9.0530,12.1400, 8.8440,8.4310,11.4650, 9.7480,8.6700,9.1130, 11.2760,9.2940,9.7450, 10.0800,11.2300,10.8340, 8.5240,10.6220,10.2170, 10.3480,10.8650,7.9290, 8.6880,12.9640,10.1330, 8.6350,14.3260,8.9130, 10.2750,12.4150,6.5770, 9.4980,14.0050,7.0000};
//     for (int i=0; i<72; i=i+3) {
//         printf("%.4f  %.4f %.4f\n", indata[i], indata[i+1], indata[i+2]);
//     }


//     double s[24*25/2];
//     for (int i=0;i<24*25/2;i++) {
//     s[i] = 0.0;
//     }


//     for (int i=0;i<24;i++) {
//         for (int j=0;j<i;j++) {
//             s[trmIndex(i,j)] = ed(indata,i,j,3);
//             //printf("%.3f " , s[trmIndex(i,j)]);
//         }
//     }
//     //printf("\n");

//     for (int i=0;i<24*25/2;i++) {
//         //printf("d: %.3f\n", s[i]);
//     }

//     double* indatap = &indata[0];

//     double* d_coords = (double*) malloc(24*2 * sizeof(double));



//     // CStochasticProximityEmbedding(
//     //     s,
//     //     d_coords,
//     //     10.0,
//     //     24,
//     //     2,
//     //     2.0,
//     //     0.1,
//     //     50,
//     //     1000000,
//     //     5);

// // void CStochasticProximityEmbedding(
// //         double* s,
// //         double* d_coords,
// //         double rco,
// //         int nelem,
// //         int dim,
// //         double maxlam,
// //         double minlam,
// //         int ncycle,
// //         int nstep,
// //         int stressfreq) {




//     CkNeighboursStochasticProximityEmbedding(
//         s,
//         d_coords,
//         10.0,
//         15,
//         24,
//         2,
//         2.0,
//         0.1,
//         50,
//         5);
// /*
//         double* s,
//         double* d_coords,
//         double rco,
//         int kn,
//         int nelem,
//         int dim,
//         double maxlam,
//         double minlam,
//         int ncycle,
//         int stressfreq) {
// */

// CkNNStochasticProximityEmbedding(s, d_coords, 15, 24, 2, 2.0, 0.1, 50, 5);

// // double CkNNStochasticProximityEmbedding(  
// //         double* s,
// //         double* d_coords,
// //         int kn,
// //         int nelem,
// //         int dim,
// //         double maxlam,
// //         double minlam,
// //         int ncycle,
// //         int stressfreq) 



// for (int i=0; i<48; i=i+2) {
//     printf("%.4lf  %.4lf\n", d_coords[i], d_coords[i+1]);
// }

// }
