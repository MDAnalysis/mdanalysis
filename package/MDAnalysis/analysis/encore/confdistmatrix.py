# confdistmatrix.py --- Conformational distance matrix calculator
# Copyright (C) 2014 Wouter Boomsma, Matteo Tiberti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Distance Matrix calculation --- :mod:`MDAnalysis.analysis.ensemble.confdistmatrix`
=====================================================================


The module contains a base class to easily compute, using parallelization and shared memory, matrices of conformational distance between the structures stored in an Ensemble.
A class to compute an RMSD matrix in such a way is also available.

"""

from multiprocessing import Process, Array, cpu_count, Value, RawValue
try:
    from MDAnalysis.analysis.rms import rmsd
    from MDAnalysis.analysis.align import rotation_matrix
except:
    from MDAnalysis.analysis.align import rmsd, rotation_matrix # backwards compatibility for MDAnalysis < 0.10.0

from numpy import sum, average, transpose, dot, ones, asarray, mean, float64, object, bool, array, int
from ctypes import c_float
from cutils import *
from getpass import getuser
from socket import gethostname
from datetime import datetime
from utils import TriangularMatrix, trm_indeces, AnimatedProgressBar
from time import sleep

class ConformationalDistanceMatrixGenerator:
    '''
    Base class for conformational distance matrices generator between array of coordinates. Work for single matrix elements is performed by the private _simple_worker and _fitter_worker methods, which respectively do or don't perform pairwise alignment before calculating the distance metric. The class efficiently and automatically spans work over a prescribed number of cores, while keeping both input coordinates and the output matrix as shared memory. If logging level is low enough, a progress bar of the whole process is printed out. This class acts as a functor.
    '''
    
    def run(self, ensemble, ncores=None, pairwise_align=False, align_subset_coordinates=None, mass_weighted=True, metadata=True):
        '''
        Run the conformational distance matrix calculation.
        
        **Arguments:**
        
        `ensemble` : encore.Ensemble.Ensemble object
            Ensemble object for which the conformational distance matrix will be computed. 
        
        `pairwise_align` : bool
            Whether to perform pairwise alignment between conformations
            
        `align_subset_coordinates` : numpy.array or None
            Use these coordinates for superimposition instead of those from ensemble.superimposition_coordinates
        
        `mass_weighted` : bool
            Whether to perform mass-weighted superimposition and metric calculation
            
        `metadata` : bool
            Whether to build a metadata dataset for the calculated matrix
            
        `ncores` : int
            Number of cores to be used for parallel calculation
        
        **Returns:**
        
        `cond_dist_matrix` : encore.utils.TriangularMatrix object
            Conformational distance matrix in triangular representation.
        '''
        
        # Decide how many cores have to be used. Since the main process is stopped while the workers do their job, ncores workers will be spawned.
        if not ncores:
            ncores = cpu_count()
        if ncores < 1:
            ncores = 1

        # framesn: number of frames
        framesn = len(ensemble.coordinates)
                
        # Prepare metadata recarray
        if metadata: 
            metadata = array([(gethostname(), getuser(), str(datetime.now()), ensemble.topology_filename, framesn, pairwise_align, ensemble.superimposition_selection_string, mass_weighted)], 
                          dtype=[('host',object),('user',object),('date',object),('topology file',object),('number of frames',int),('pairwise superimposition',bool),('superimposition subset',object),('mass-weighted',bool)])

        # Prepare alignment subset coordinates as necessary
        subset_coords = None 
        if pairwise_align:
            subset_selection = ensemble.superimposition_selection
            if align_subset_coordinates == None:
                subset_coords = align_subset_coordinates
            else:
                subset_coords = ensemble.superimposition_coordinates

        # Prepare masses as necessary
        subset_masses = None
    
        if mass_weighted:
            masses = ensemble.atom_selection.masses
            if pairwise_align:
                subset_masses = subset_selection.masses
        else:
            masses = ones((ensemble.coordinates[0].shape[0]))
            if pairwise_align:
                subset_masses = ones((subset_coords[0].shape[0]))

        # matsize: number of elements of the triangular matrix, diagonal elements included.
        matsize = framesn*(framesn+1)/2

        # Calculate the number of matrix elements that each core has to calculate as equally as possible. 
        if ncores > matsize:
            ncores = matsize
        runs_per_worker = [ matsize / int(ncores) for x in range(ncores) ]
        unfair_work = matsize % ncores
        for i in range(unfair_work):
            runs_per_worker[i] += 1
            
        # Splice the matrix in ncores segments. Calculate the first and the last (i,j) 
        # matrix elements of the slices that will be assigned to each worker. Each of them will proceed in a column-then-row order
        # (e.g. 0,0 1,0 1,1 2,0 2,1 2,2 ... )
        i=0
        a=[0,0]
        b=[0,0]
        tasks_per_worker = []        
        for n in range(len(runs_per_worker)):
            while i*(i-1)/2 < sum(runs_per_worker[:n+1]):
                i += 1
            b = [ i-2, sum(runs_per_worker[0:n+1])-(i-2)*(i-1)/2-1 ]
            tasks_per_worker.append((tuple(a), tuple(b)))
            if b[0] == b[1]: 
                a[0] = b[0] + 1
                a[1] = 0
            else:
                a[0] = b[0]
                a[1] = b[1] + 1
        
        # Allocate for output matrix
        distmat = Array(c_float, matsize) 

	# Prepare progress bar stuff and run it
	pbar = AnimatedProgressBar(end=matsize, width=80)
	partial_counters = [RawValue('i',0) for i in range(ncores)]

        # Initialize workers. Simple worker doesn't perform fitting, fitter worker does.      
        if pairwise_align:
            workers = [Process(target=self._fitter_worker, args=(tasks_per_worker[i], ensemble.coordinates, subset_coords, masses, subset_masses, distmat, partial_counters[i])) for i in range(ncores)]
        else:
            workers = [Process(target=self._simple_worker, args=(tasks_per_worker[i], ensemble.coordinates, masses, distmat, pbar_counter)) for i in range(ncores)]

	workers += [Process(target=self._pbar_updater, args=(pbar, partial_counters, matsize))]

        # Start & join the workers
        for w in workers:
            w.start()
        for w in workers:
            w.join()
        
        # When the workers have finished, return a TriangularMatrix object
        return TriangularMatrix(distmat,metadata=metadata)

    def _simple_worker(self, tasks, coords, masses, rmsdmat, pbar_counter):
        '''Simple worker prototype; to be overriden in derived classes
        '''
        for i,j in trm_indeces(tasks[0],tasks[1]):
            pass
                
    def _fitter_worker(self, tasks, coords, subset_coords, masses, subset_masses, rmsdmat, pbar_counter): # Prototype fitter worker: pairwase align and calculate metric. To be ovverridden in heir classes
        '''Fitter worker prototype; to be overriden in derived classes
        '''

        if subset_coords == None:
            for i,j in trm_indeces(tasks[0],tasks[1]):
                coords[i] -= average(coords[i], axis=0, weights=masses)
                coords[j] -= average(coords[j], axis=0, weights=masses)
		pbar_counter.value += 1
                pass
        else:
            for i,j in trm_indeces(tasks[0],tasks[1]):
                com_i = average(coords[i], axis=0, weights=masses)
                translated_i = coords[i] - com_i
                subset1_coords = subset_coords[i] - com_i
                com_j = average(coords[j], axis=0, weights=masses)
                translated_j = coords[j] - com_j
                subset2_coords = subset_coords[j] - com_j
                rotamat = rotation_matrix(subset1_coords, subset2_coords, subset_masses)[0]
                rotated_i = transpose(dot(rotamat, transpose(translated_i)))
                pbar_counter.value += 1
                pass

    def _pbar_updater(self, pbar,  pbar_counters, max_val, update_interval=0.2):
        '''Method that updates and prints the progress bar, upon polling progress status from workers.
            
        **Attributes:**
            
        `pbar` : encore.utils.AnimatedProgressBar object
            Progress bar object
        
        `pbar_counters` : list of multiprocessing.RawValue
            List of counters. Each worker is given a counter, which is updated at every cycle. In this way the _pbar_updater process can asynchronously fetch progress reports.
        
	`max_val` : int
            Total number of matrix elements to be calculated
        
	`update_interval` : float
            Number of seconds between progress bar updates
            
            '''
        
        
        
        val = 0 
        while val < max_val:
            val = 0
            for c in pbar_counters:
                val += c.value
            pbar.update(val)
            pbar.show_progress()
            sleep(update_interval)

    __call__ = run
	
class RMSDMatrixGenerator(ConformationalDistanceMatrixGenerator):
    '''
        RMSD Matrix calculator. Simple workers doesn't perform fitting, while fitter worker does.
    '''
    def _simple_worker(self, tasks, coords, masses, rmsdmat, pbar_counter):
        '''
        Simple RMSD Matrix calculator.
            
        **Arguments:**
        
        `tasks` : iterator of int of length 2
            Given a triangular matrix, this worker will calculate RMSD values from element tasks[0] to tasks[1]. Since the matrix is triangular, the trm_indeces matrix automatically calculates the corrisponding i,j matrix indices. The matrix is written as an array in a row-major order (see the TriangularMatrix class for details).
        
	`coords` : numpy.array
            Array of the ensemble coordinates
        
	`masses` : numpy.array
            Array of atomic masses, having the same order as the coordinates array
        
	`rmsdmat` : encore.utils.TriangularMatrix
            Memory-shared triangular matrix object
        
	`pbar_counter` : multiprocessing.RawValue
            Thread-safe shared value. This counter is updated at every cycle and used to evaluate the progress of each worker.
            '''
        for i,j in trm_indeces(tasks[0],tasks[1]):
            #masses = asarray(masses)/mean(masses)
            summasses = sum(masses)
            rmsdmat[(i+1)*i/2+j] = PureRMSD(coords[i].astype(float64), coords[j].astype(float64), coords[j].shape[0], masses, summasses)
            pbar_counter.value += 1

    def _fitter_worker(self, tasks, coords, subset_coords, masses, subset_masses, rmsdmat, pbar_counter):
        '''
            Fitter RMSD Matrix calculator: performs least-square fitting between each pair of structures before calculating the RMSD.
            
            **Arguments:**
            
            `tasks` : iterator of int of length 2
            Given a triangular matrix written in a row-major order, this worker will calculate RMSD values from element tasks[0] to tasks[1]. Since the matrix is triangular. the trm_indeces function automatically calculates the corrisponding i,j matrix indeces. (see the see encore.utils.TriangularMatrix for details).

            `coords` : numpy.array
                Array of the ensemble coordinates

            `subset_coords` : numpy.array or None
                Array of the coordinates used for fitting

            `masses` : numpy.array or None
                Array of atomic masses, having the same order as the coordinates array. If None, coords will be used instead.

            `subset_masses` : numpy.array
                Array of atomic masses, having the same order as the subset_coords array

            `rmsdmat` : encore.utils.TriangularMatrix
                Memory-shared triangular matrix object

            `pbar_counter` : multiprocessing.RawValue
                Thread-safe shared value. This counter is updated at every cycle and used to evaluate the progress of each worker.
            '''

        if subset_coords == None:
            for i,j in trm_indeces(tasks[0],tasks[1]):
                coords[i] -= average(coords[i], axis=0, weights=masses)
                coords[j] -= average(coords[j], axis=0, weights=masses)
                weights = asarray(masses)/mean(masses)
                rmsdmat[(i+1)*i/2+j] = rmsd(coords[i],coords[j],weights=weights)
                pbar_counter.value += 1
        else:
            for i,j in trm_indeces(tasks[0],tasks[1]):
                summasses = sum(masses)
                subset_weights = asarray(subset_masses)/mean(subset_masses)
                com_i = average(subset_coords[i], axis=0, weights=subset_masses)
                translated_i = coords[i] - com_i
                subset1_coords = subset_coords[i] - com_i
                com_j = average(subset_coords[j], axis=0, weights=subset_masses)
                translated_j = coords[j] - com_j
                subset2_coords = subset_coords[j] - com_j
                rotamat = rotation_matrix(subset1_coords, subset2_coords, subset_weights)[0]
                rotated_i = transpose(dot(rotamat, transpose(translated_i)))
                rmsdmat[(i+1)*i/2+j] = PureRMSD(rotated_i.astype(float64), translated_j.astype(float64), coords[j].shape[0], masses, summasses)            
                pbar_counter.value += 1

class MinusRMSDMatrixGenerator(ConformationalDistanceMatrixGenerator):
    '''
        -RMSD Matrix calculator. See encore.confdistmatrix.RMSDMatrixGenerator for details.
    '''

    def _simple_worker(self, tasks, coords, masses, rmsdmat, pbar_counter):
        '''
            Simple RMSD Matrix calculator. See encore.confdistmatrix.RMSDMatrixGenerator._simple_worker for details.
        '''
        for i,j in trm_indeces(tasks[0],tasks[1]):
            #masses = asarray(masses)/mean(masses)
            summasses = sum(masses)
            rmsdmat[(i+1)*i/2+j] = MinusRMSD(coords[i].astype(float64), coords[j].astype(float64), coords[j].shape[0], masses, summasses)            
            pbar_counter.value += 1

    def _fitter_worker(self, tasks, coords, subset_coords, masses, subset_masses, rmsdmat, pbar_counter):
        '''
        Fitter RMSD Matrix calculator. See encore.confdistmatrix.RMSDMatrixGenerator._fitter_worker for details.
        '''

        if subset_coords == None:
            for i,j in trm_indeces(tasks[0],tasks[1]):
                coords[i] -= average(coords[i], axis=0, weights=masses)
                coords[j] -= average(coords[j], axis=0, weights=masses)
                weights = asarray(masses)/mean(masses)
                rmsdmat[(i+1)*i/2+j] = - rmsd(coords[i],coords[j],weights=weights)
                pbar_counter.value += 1
        else:
            for i,j in trm_indeces(tasks[0],tasks[1]):
                #masses = asarray(masses)/mean(masses)
                summasses = sum(masses)
                com_i = average(subset_coords[i], axis=0, weights=subset_masses)
                translated_i = coords[i] - com_i
                subset1_coords = subset_coords[i] - com_i
                com_j = average(subset_coords[j], axis=0, weights=subset_masses)
                translated_j = coords[j] - com_j
                subset2_coords = subset_coords[j] - com_j
                rotamat = rotation_matrix(subset1_coords, subset2_coords, subset_masses)[0]
                rotated_i = transpose(dot(rotamat, transpose(translated_i)))
                rmsdmat[(i+1)*i/2+j] = MinusRMSD(rotated_i.astype(float64), translated_j.astype(float64), coords[j].shape[0], masses, summasses)   
                pbar_counter.value += 1

