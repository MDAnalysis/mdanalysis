/*
spe.h --- header file for spe.c
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

int trmIndex(int, int);

double ed(double*, int, int, int);

double stress(double, double, int, int);

double neighbours_stress(double, double, int, int, double);

int neighbours(double, int, double, int*, int*, int*);

int cmp_ivwrapper(const void*, const void*);

int nearest_neighbours(double*, int, int);

double CkNNStochasticProximityEmbedding(
        double*,
        double*,
        int,
        int,
        int,
        double,
        double,
        int,
        int,
        int); 

double CkNeighboursStochasticProximityEmbedding(
	double*,
        double*,
        double,
        int,
        int,
        int,
        double,
        double,
        int,
        int);

double CStochasticProximityEmbedding(
        double*,
        double*,
        double,
        int,
        int,
        double,
        double,
        int,
        int,
        int);




