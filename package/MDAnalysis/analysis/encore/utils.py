# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from multiprocessing.sharedctypes import SynchronizedArray
from multiprocessing import Process, Manager
import numpy as np
import sys
import MDAnalysis as mda
from ...coordinates.memory import MemoryReader


class TriangularMatrix(object):
    """Triangular matrix class. This class is designed to provide a
    memory-efficient representation of a triangular matrix that still behaves
    as a square symmetric one. The class wraps a numpy.array object,
    in which data are memorized in row-major order. It also has few additional
    facilities to conveniently load/write a matrix from/to file. It can be
    accessed using the [] and () operators, similarly to a normal numpy array.

    """

    def __init__(self, size, metadata=None, loadfile=None):
        """Class constructor.

        Parameters
        ----------

        size : int or multiprocessing.SyncrhonizeArray
            Size of the matrix (number of rows or columns). If an
            array is provided instead, the size of the triangular matrix
            will be calculated and the array copied as the matrix
            elements. Otherwise, the matrix is just initialized to zero.

        metadata : dict or None
            Metadata dictionary. Used to generate the metadata attribute.

        loadfile : str or None
            Load the matrix from this file. All the attributes and data will
            be determined by the matrix file itself (i.e. metadata will be
            ignored); size has to be provided though.

        """
        self.metadata = metadata
        self.size = size
        if loadfile:
            self.loadz(loadfile)
            return
        if type(size) == int:
            self.size = size
            self._elements = np.zeros((size + 1) * size / 2, dtype=np.float64)
            return
        if type(size) == SynchronizedArray:
            self._elements = np.array(size.get_obj(), dtype=np.float64)
            self.size = int((np.sqrt(1 + 8 * len(size)) - 1) / 2)
            return
        else:
            raise TypeError

    def __getitem__(self, args):
        x, y = args
        if x < y:
            x, y = y, x
        return self._elements[x * (x + 1) / 2 + y]

    def __setitem__(self, args, val):
        x, y = args
        if x < y:
            x, y = y, x
        self._elements[x * (x + 1) / 2 + y] = val

    def as_array(self):
        """Return standard numpy array equivalent"""
        a = np.zeros((self.size, self.size))
        a[np.tril_indices(self.size)] = self._elements
        a[np.triu_indices(self.size)] = a.T[np.triu_indices(self.size)]
        return a

    def savez(self, fname):
        """Save matrix in the npz compressed numpy format. Save metadata and
        data as well.

        Parameters
        ----------

        fname : str
            Name of the file to be saved.
        """
        np.savez(fname, elements=self._elements, metadata=self.metadata)

    def loadz(self, fname):
        """Load matrix from the npz compressed numpy format.

        Parameters
        ----------

        fname : str
            Name of the file to be loaded.
        """
        loaded = np.load(fname)

        if loaded['metadata'].shape != ():
            if loaded['metadata']['number of frames'] != self.size:
                raise TypeError
            self.metadata = loaded['metadata']
        else:
            if self.size*(self.size-1)/2+self.size != len(loaded['elements']):
                raise TypeError
        self._elements = loaded['elements']

    def __add__(self, scalar):
        """Add scalar to matrix elements.

        Parameters
        ----------

        scalar : float
            Scalar to be added.
        """
        newMatrix = self.__class__(self.size)
        newMatrix._elements = self._elements + scalar;
        return newMatrix

    def __iadd__(self, scalar):
        """Add scalar to matrix elements.

        Parameters
        ----------

        scalar : float
            Scalar to be added.
        """
        self._elements += scalar
        return self


    def __mul__(self, scalar):
        """Multiply with scalar.

        Parameters
        ----------

        scalar : float
            Scalar to multiply with.
        """
        newMatrix = self.__class__(self.size)
        newMatrix._elements = self._elements * scalar;
        return newMatrix

    def __imul__(self, scalar):
        """Multiply with scalar.

        Parameters
        ----------

        scalar : float
            Scalar to multiply with.
        """
        self._elements *= scalar
        return self




    __rmul__ = __mul__

    def __str__(self):
        return str(self.as_array())


class ParallelCalculation(object):
    """
    Generic parallel calculation class. Can use arbitrary functions,
    arguments to functions and kwargs to functions.

    Attributes
    ----------

    ncores : int
            Number of cores to be used for parallel calculation

    function : callable object
            Function to be run in parallel.

    args : list of tuples
            Each tuple contains the arguments that will be passed to
            function(). This means that a call to function() is performed for
            each tuple. function is called as function(\*args, \*\*kwargs). Runs
            are distributed on the requested numbers of cores.

    kwargs : list of dicts
            Each tuple contains the named arguments that will be passed to
            function, similarly as described for the args attribute.

    nruns : int
            Number of runs to be performed. Must be equal to len(args) and
            len(kwargs).
    """

    def __init__(self, ncores, function, args=None, kwargs=None):
        """ Class constructor.

        Parameters
        ----------

        ncores : int
            Number of cores to be used for parallel calculation

        function : object that supports __call__, as functions
            function to be run in parallel.

        args : list of tuples
            Arguments for function; see the ParallelCalculation class
            description.

        kwargs : list of dicts or None
            kwargs for function; see the ParallelCalculation
            class description.
        """

        # args[i] should be a list of args, one for each run
        self.ncores = ncores
        self.functions = function
        if not hasattr(self.functions, '__iter__'):
            self.functions = [self.functions]*len(args)
        if len(self.functions) != len(args):
            self.functions = self.functions[:]*(len(args)/len(self.functions))

        # Arguments should be present
        if args is None:
            args = []
        self.args = args

        # If kwargs are not present, use empty dicts
        if kwargs:
            self.kwargs = kwargs
        else:
            self.kwargs = [{} for i in self.args]

        self.nruns = len(args)

    def worker(self, q, results):
        """
        Generic worker. Will run function with the prescribed args and kwargs.

        Parameters
        ----------

        q : multiprocessing.Manager.Queue object
                work queue, from which the worker fetches arguments and
                messages

        results : multiprocessing.Manager.Queue object
                results queue, where results are put after each calculation is
                finished

        """
        while True:
            i = q.get()
            if i == 'STOP':
                return

            results.put((i, self.functions[i](*self.args[i], **self.kwargs[i])))

    def run(self):
        """
        Run parallel calculation.

        Returns
        -------

        results : tuple of ordered tuples (int, object)
                int is the number of the calculation corresponding to a
                certain argument in the args list, and object is the result of
                corresponding calculation. For instance, in (3, output), output
                is the return of function(\*args[3], \*\*kwargs[3]).
        """
        results_list = []
        if self.ncores == 1:
            for i in range(self.nruns):
                results_list.append((i, self.functions[i](*self.args[i],
                                                          **self.kwargs[i])))
        else:
            manager = Manager()
            q = manager.Queue()
            results = manager.Queue()

            workers = [Process(target=self.worker, args=(q, results)) for i in
                       range(self.ncores)]

            for i in range(self.nruns):
                q.put(i)
            for w in workers:
                q.put('STOP')

            for w in workers:
                w.start()

            for w in workers:
                w.join()

            results.put('STOP')
            for i in iter(results.get, 'STOP'):
                results_list.append(i)

        return tuple(sorted(results_list, key=lambda x: x[0]))


class ProgressBar(object):
    """Handle and draw a progress barr.
    From https://github.com/ikame/progressbar
    """

    def __init__(self, start=0, end=10, width=12, fill='=', blank='.',
                 format='[%(fill)s>%(blank)s] %(progress)s%%',
                 incremental=True):
        super(ProgressBar, self).__init__()

        self.start = start
        self.end = end
        self.width = width
        self.fill = fill
        self.blank = blank
        self.format = format
        self.incremental = incremental
        self.step = 100 / float(width)  # fix
        self.reset()

    def __add__(self, increment):
        increment = self._get_progress(increment)
        if 100 > self.progress + increment:
            self.progress += increment
        else:
            self.progress = 100
        return self

    def __str__(self):
        progressed = int(self.progress / self.step)  # fix
        fill = progressed * self.fill
        blank = (self.width - progressed) * self.blank
        return self.format % {'fill': fill, 'blank': blank,
                              'progress': int(self.progress)}

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
    Accepts an extra keyword argument named `stdout`
    (by default use sys.stdout).
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


def trm_indeces(a, b):
    """
    Generate (i,j) indeces of a triangular matrix, between elements a and b.
    The matrix size is automatically determined from the number of elements.
    For instance: trm_indeces((0,0),(2,1)) yields (0,0) (1,0) (1,1) (2,0)
    (2,1).

    Parameters
    ----------

    a : (int i, int j) tuple
        starting matrix element.

    b : (int i, int j) tuple
        final matrix element.
    """
    i, j = a
    while i < b[0]:
        if i == j:
            yield (i, j)
            j = 0
            i += 1
        else:
            yield (i, j)
            j += 1
    while j <= b[1]:
        yield (i, j)
        j += 1


def trm_indices_nodiag(n):
    """generate (i,j) indeces of a triangular matrix of n rows (or columns),
    without diagonal (e.g. no elements (0,0),(1,1),...,(n,n))

    Parameters
    ----------

    n : int
        Matrix size
"""

    for i in xrange(1, n):
        for j in xrange(i):
            yield (i, j)


def trm_indices_diag(n):
    """generate (i,j) indeces of a triangular matrix of n rows (or columns),
    with diagonal

    Parameters
    ----------

    n : int
        Matrix size
"""

    for i in xrange(0, n):
        for j in xrange(i+1):
            yield (i, j)


def merge_universes(universes):
    """
    Merge list of universes into one

    Parameters
    ----------
    universes : list of Universe objects


    Returns
    ----------
    Universe object
    """

    for universe in universes:
        universe.transfer_to_memory()

    return mda.Universe(
        universes[0].filename,
        np.concatenate(tuple([e.trajectory.timeseries(format='fac') for e in universes]),
                       axis=0),
        format=MemoryReader)
