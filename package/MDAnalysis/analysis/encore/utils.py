# utils.py --- Utility functions for encore package
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


from multiprocessing.sharedctypes import SynchronizedArray
from multiprocessing import Process, Manager
from numpy import  savez, load, zeros, array, float64, sqrt, atleast_2d, reshape, newaxis, zeros, dot, sum, exp
import numpy as np
from scipy.stats import gaussian_kde
import sys 
import time
import optparse
import copy


class TriangularMatrix:
    """Triangular matrix class. This class is designed to provide a memory-efficient representation of a triangular matrix that still behaves as a square symmetric one. The class wraps a numpy.array object, in which data are memorized in row-major order. It also has few additional facilities to conveniently load/write a matrix from/to file. It can be accessed using the [] and () operators, similarly to a normal numpy array.

**Attributes:**

`size` : int
	Size of the matrix (number of rows or number of columns)

`metadata` : dict
	Metadata for the matrix (date of creation, name of author ...)
"""    
    
    
    def __init__(self, size, metadata=None, loadfile=None):
        """Class constructor.
        
	**Attributes:**

	`size` : int or multiprocessing.SyncrhonizeArray
        Size of the matrix (number of rows or columns). If an array is provided instead, the size of the triangular matrix will be calculated and the array copied as the matrix elements. Otherwise, the matrix is just initialized to zero.

        `metadata` : dict or None
        Metadata dictionary. Used to generate the metadata attribute.

        `loadfile` : str or None
        Load the matrix from this file. All the attributes and data will be determined by the matrix file itself (i.e. metadata will be ignored); size has to be provided though.
        """
        self.metadata = metadata
        self.size = size
        if loadfile:
            self.loadz(loadfile)
            return
        if type(size) == int:
            self.size = size
            self._elements = zeros((size+1)*size/2, dtype=float64)
            return
        if type(size) == SynchronizedArray:
            self._elements = array(size.get_obj(), dtype=float64)
            self.size = int((sqrt(1+8*len(size))-1)/2)
            return
        else:
            raise TypeError
        
    def __call__(self,x,y):
        if x < y:
            x, y = y, x
        return self._elements[x*(x+1)/2+y]
        
    def __getitem__(self,args):
        x, y = args
        if x < y:
            x, y = y, x
        return self._elements[x*(x+1)/2+y]
        
    def __setitem__(self,args,val):
        x, y = args
        if x < y:
            x, y = y, x
        self._elements[x*(x+1)/2+y] = val
                                
    def savez(self,fname):
        """Save matrix in the npz compressed numpy format. Save metadata and data as well.

        **Arguments**:
        
        `fname` : str
        Name of the file to be saved.

        """
        savez(fname, elements=self._elements,metadata=self.metadata)
        
    def loadz(self,fname):
        """Load matrix from the npz compressed numpy format. 

        **Arguments**:

        `fname` : str
        Name of the file to be loaded.
        """
        loaded = load(fname)
        if loaded['metadata'] != None:
            if loaded['metadata']['number of frames'] != self.size:
                raise TypeError
            self.metadata = loaded['metadata']
        else:
            if self.size != len(loaded['elements']):
                raise TypeError
        self._elements = loaded['elements']

    def trm_print(self, justification=10):
        """ 
        Print the triangular matrix as triangular
        """
        for i in xrange(0,self.size):
            for j in xrange(i+1):
                print "%.3f".ljust(justification) % self.__getitem__((i,j)),
            print ""



class ParallelCalculation:
    """
    Generic parallel calculation class. Can use arbitrary functions,
    arguments to functions and kwargs to functions. 

    **Attributes:**
	`ncores` : int 
		Number of cores to be used for parallel calculation
	
	`function` : callable object
		Function to be run in parallel.

	`args` : list of tuples
		Each tuple contains the arguments that will be passed to function(). This means that a call to function() is performed for each tuple. function is called as function(*args, **kwargs). Runs are distributed on the requested numbers of cores.

	`kwargs` : list of dicts
		Each tuple contains the named arguments that will be passed to function, similarly as described for the args attribute.

	`nruns` : int
		Number of runs to be performed. Must be equal to len(args) and len(kwargs).
    """
    def __init__(self, ncores, function, args=[], kwargs=None):
	""" Class constructor.

	**Arguments:**
	
	`ncores` : int 
		Number of cores to be used for parallel calculation

	`function` : object that supports __call__, as functions
		function to be run in parallel.

	`args` : list of tuples
		Arguments for function; see the ParallelCalculation class description.

	`kwargs` : list of dicts or None
		kwargs for function; see the ParallelCalculation class description.
	"""
	
        # args[0] should be a list of args, one for each run
        self.ncores = ncores
        self.function = function

        # Arguments should be present
        self.args = args
            
        # If kwargs are not present, use empty dicts
        if kwargs:
            self.kwargs = kwargs
        else:
            self.kwargs = [ {} for i in self.args ]            

        self.nruns  = len(args)

    def worker(self, q, results):
        """
        Generic worker. Will run function with the prescribed args and kwargs.

	**Arguments:**

	`q` : multiprocessing.Manager.Queue object
		work queue, from which the worker fetches arguments and messages

	`results` : multiprocessing.Manager.Queue object
		results queue, where results are put after each calculation is finished

        """
        while True:
            i = q.get()
            if i == 'STOP':
                return
            results.put( (i,self.function(*self.args[i],**self.kwargs[i]) ) )

    def run(self):
        """
        Run parallel calculation.

	**Returns:**

	`results` : tuple of ordered tuples (int, object)
		int is the number of the calculation corresponding to a certain argument in the args list, and object is the result of corresponding calculation. For instance, in (3, output), output is the return of function(\*args[3], \*\*kwargs[3]).  
        """
        manager = Manager()
        q = manager.Queue()
        results = manager.Queue()

        workers = [ Process(target=self.worker, args=(q,results)) for i in range(self.ncores) ]
        
        for w in workers:
            w.start()
        
        for i in range(self.nruns):
            q.put(i)
        for w in workers:
            q.put('STOP')
                
        for w in workers:
            w.join()
        
        results_list = []
        
        results.put('STOP')
        for i in iter(results.get, 'STOP'):
            results_list.append(i)
        
        return tuple(sorted(results_list, key=lambda x: x[0]))

class ProgressBar(object):
    """Handle and draw a progress barr.  From https://github.com/ikame/progressbar
    """
    def __init__(self, start=0, end=10, width=12, fill='=', blank='.', format='[%(fill)s>%(blank)s] %(progress)s%%', incremental=True):
        super(ProgressBar, self).__init__()

        self.start = start
        self.end = end
        self.width = width
        self.fill = fill
        self.blank = blank
        self.format = format
        self.incremental = incremental
        self.step = 100 / float(width) #fix
        self.reset()

    def __add__(self, increment):
        increment = self._get_progress(increment)
        if 100 > self.progress + increment:
            self.progress += increment
        else:
            self.progress = 100
        return self

    def __str__(self):
        progressed = int(self.progress / self.step) #fix
        fill = progressed * self.fill
        blank = (self.width - progressed) * self.blank
        return self.format % {'fill': fill, 'blank': blank, 'progress': int(self.progress)}

    __repr__ = __str__

    def _get_progress(self, increment):
        return float(increment * 100) / self.end

    def reset(self):
        """Resets the current progress to the start point"""
        self.progress = self._get_progress(self.start)
        return self
    def update(self, progress):
	"""Update the progress value instead of incrementing it"""
	this_progress = self._get_progress(progress)
	if this_progress < 100:
	    self.progress = this_progress
        else:
            self.progress = 100
	
class AnimatedProgressBar(ProgressBar):
    """Extends ProgressBar to allow you to use it straighforward on a script.
    Accepts an extra keyword argument named `stdout` (by default use sys.stdout).
    The progress status may be send to any file-object. 
    """
    def __init__(self, *args, **kwargs):
        super(AnimatedProgressBar, self).__init__(*args, **kwargs)
        self.stdout = kwargs.get('stdout', sys.stdout)

    def show_progress(self):
        if hasattr(self.stdout, 'isatty') and self.stdout.isatty():
            self.stdout.write('\r')
        else:
            self.stdout.write('\n')
        self.stdout.write(str(self))
        self.stdout.flush()


def trm_indeces(a,b):
    """
    Generate (i,j) indeces of a triangular matrix, between elements a and b. The matrix size is automatically determined from the number of elements.
    For instance: trm_indexes((0,0),(2,1)) yields (0,0) (1,0) (1,1) (2,0) (2,1).

    **Arguments:**

    `a` : (int i, int j) tuple
	starting matrix element. 

    `b` : (int i, int j) tuple
	final matrix element.

    """
    i, j = a
    while i < b[0]:
        if i == j:
            yield (i,j)
            j = 0
            i += 1
        else:
            yield (i,j)
            j += 1
    while j <= b[1]:
        yield (i,j)
        j+=1

def trm_indeces_nodiag(n):
    """generate (i,j) indeces of a triangular matrix of n rows (or columns), without diagonal (e.g. no elements (0,0),(1,1),...,(n,n))

    **Arguments:**

    `n` : int 
	Matrix size
"""

    for i in xrange(1,n):
        for j in xrange(i):
            yield (i,j)