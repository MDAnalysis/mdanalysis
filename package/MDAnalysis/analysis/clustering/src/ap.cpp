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
#include <float.h>
#include <vector>
#include <algorithm>
#include <functional>

// inline int tri_i(int row_i, int col_i) {
//     return (row_i)>(col_i) ? ((row_i)+1)*(row_i)/2+(col_i) : ((col_i)+1)*(col_i)/2+(row_i);
// }

inline int get_idx(int row_i, int col_i, int n_col) {
	return row_i * n_col + col_i;
}

int CAffinityPropagation(float* sim, int length, float lambda,
                        int max_iter, int conv_threshold, int add_noise,
                        long* clusters) {
	
	float* res = new float[length*length];
	float* av = new float[length*length];
	float* simav = new float[length*length];
	

	if (add_noise > 0) {
		for (int i=0; i < length; i++) {
			for (int j = 0; j < i; j++) {

				int idx = get_idx(i, j, length);
				float val = sim[idx];
				float noise = 1e-16 * val * (rand() / (float)RAND_MAX+1);
				sim[idx] += noise;
				sim[get_idx(j, i, length)] += noise;
			}
		}
	} 

	// for (int i = 0; i < length; i++) {
	// 	for (int j = 0; j < i; j++) {
	// 		int trix = tri_i(i, j);
	// 		sim[idx] = similarity[trix];
	// 		sim[j][i] = similarity[trix];
	// 	}
	// }

	// initialize exemplars
	int exemplars[length];
	int old_exemplars[length];
	std::fill_n(exemplars, length, -1);

	// main loop
	int iter = 0;
	float lambda_comp = 1.0 - lambda;
	bool converged = false;
	while (iter < max_iter && !converged) {
		// update responsibilities
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < length; j++) {
				int idx = get_idx(i, j, length);
				simav[idx] = sim[idx] + av[idx];
			}
			std::sort(&simav[0], &simav[length]);
			float *unique = std::unique(&simav[0], &simav[length]);
			std::nth_element(&simav[0], &simav[2], unique, std::greater<float>());
			float max1 = simav[0];
			float max2 = simav[1];

			for (int j = 0; j < length; j++) {
				int idx = get_idx(i, j, length);
				if (av[idx] + sim[idx] == max1) {
					res[idx] = lambda * av[idx] + lambda_comp * (sim[idx] - max2);
				} else {
					res[idx] = lambda * av[idx] + lambda_comp * (sim[idx] - max1);
				}
			}
		}
		// update availabilities
		for (int j = 0; j < length; j++) {
			float rsum = 0.0;
			for (int i = 0; i < length; i++) {
				int idx = get_idx(i, j, length);
				float r = res[idx];
				if (i == j || r > 0.0) {
					rsum += r;
				}
			}
			for (int i = 0; i < length; i++) {
				int idx = get_idx(i, j, length);
				av[idx] *= lambda;
				if (i != j) {
					float r = res[idx];
					if (r > 0.0) {
						rsum -= r;
					}
					if (rsum >= 0.0) {
						av[idx] += lambda_comp * rsum;
					}
				} else {
					av[idx] = lambda_comp * (rsum - res[get_idx(j, j, length)]);
				}
			}
		}

		// convergence check
		int *tmp_exemplars = old_exemplars;
		int *old_exemplars = exemplars;
		int *exemplars = tmp_exemplars;

		bool has_cluster = false;
		for (int i = 0; i < length; i++) {
			int idx = get_idx(i, i, length);
			if (res[idx] + av[idx] > 0.0) {
				exemplars[i] = 1;
				has_cluster = true;
			} else {
				exemplars[i] = 0;
			}
		}

		int conv_count = 0;
		if (has_cluster) {
			conv_count++;
			for (int j = 0; j < length; j++) {
				if (exemplars[j] != old_exemplars[j]) {
					conv_count = 0;
					break;
				}
			}
		}

		if (conv_count == conv_threshold) {
			converged = true;
		}
		iter++;
	}

	// find number of clusters
	int n_clusters = 0;
	for (int i = 0; i<length; i++) {
		int idx = get_idx(i, i, length);
		if (res[idx] + av[idx] > 0) {
			exemplars[n_clusters] = i;
			n_clusters++;
		}
	}

	for (int i=0; i < length; i++) {
		float maxsim = res[i*length] + av[i*length];
		for (int j = 0; j < length; j++) {
			int idx = get_idx(i, j, length);
			float tmpsum = res[idx] + av[idx];
			if (tmpsum >= maxsim) {
				clusters[i] = j;
				maxsim = tmpsum;
			}
		}
	}

	for (int i = 0; i<n_clusters; i++) {
		clusters[exemplars[i]] = exemplars[i];
	}

	delete [] res;
	delete [] av;
	delete [] simav;

	return converged ? iter : -iter;
	
	}

