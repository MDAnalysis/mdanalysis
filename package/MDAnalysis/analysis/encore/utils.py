# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import numbers
from multiprocessing.sharedctypes import SynchronizedArray
from multiprocessing import Process, Manager
from joblib import cpu_count
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

        size : int / array_like
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
        if isinstance(metadata, dict):
            self.metadata = np.array(metadata.items(), dtype=object)
        else:
            self.metadata = metadata

        self.size = size
        if loadfile:
            self.loadz(loadfile)
        elif isinstance(size, numbers.Integral):
            self.size = size
            self._elements = np.zeros((size + 1) * size // 2, dtype=np.float64)
        elif isinstance(size, SynchronizedArray):
            self._elements = np.array(size.get_obj(), dtype=np.float64)
            self.size = int((np.sqrt(1 + 8 * len(size)) - 1) / 2)
        elif isinstance(size, np.ndarray):
            self._elements = size
            self.size = int((np.sqrt(1 + 8 * len(size)) - 1) / 2)
        else:
            raise TypeError

    def __getitem__(self, args):
        x, y = args
        if x < y:
            x, y = y, x
        return self._elements[x * (x + 1) // 2 + y]

    def __setitem__(self, args, val):
        x, y = args
        if x < y:
            x, y = y, x
        self._elements[x * (x + 1) // 2 + y] = val

    def as_array(self):
        """Return standard numpy array equivalent"""
        # pylint: disable=unsubscriptable-object
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
        loaded = np.load(fname, allow_pickle=True)

        if loaded["metadata"].shape != ():
            if loaded["metadata"]["number of frames"] != self.size:
                raise TypeError
            self.metadata = loaded["metadata"]
        else:
            if self.size * (self.size - 1) / 2 + self.size != len(
                loaded["elements"]
            ):
                raise TypeError
        self._elements = loaded["elements"]

    def __add__(self, scalar):
        """Add scalar to matrix elements.

        Parameters
        ----------

        scalar : float
            Scalar to be added.
        """
        newMatrix = self.__class__(self.size)
        newMatrix._elements = self._elements + scalar
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
        newMatrix._elements = self._elements * scalar
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
    r"""
    Generic parallel calculation class. Can use arbitrary functions,
    arguments to functions and kwargs to functions.

    Attributes
    ----------
    n_jobs : int
        Number of cores to be used for parallel calculation. If -1 use all
        available cores.
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

    def __init__(self, n_jobs, function, args=None, kwargs=None):
        """
        Parameters
        ----------
        n_jobs : int
            Number of cores to be used for parallel calculation. If -1 use all
            available cores.
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
        self.n_jobs = n_jobs
        if self.n_jobs == -1:
            self.n_jobs = cpu_count()

        self.functions = function
        if not hasattr(self.functions, "__iter__"):
            self.functions = [self.functions] * len(args)
        if len(self.functions) != len(args):
            self.functions = self.functions[:] * (
                len(args) // len(self.functions)
            )

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
            if i == "STOP":
                return

            results.put(
                (i, self.functions[i](*self.args[i], **self.kwargs[i]))
            )

    def run(self):
        r"""
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
        if self.n_jobs == 1:
            for i in range(self.nruns):
                results_list.append(
                    (i, self.functions[i](*self.args[i], **self.kwargs[i]))
                )
        else:
            manager = Manager()
            q = manager.Queue()
            results = manager.Queue()

            workers = [
                Process(target=self.worker, args=(q, results))
                for i in range(self.n_jobs)
            ]

            for i in range(self.nruns):
                q.put(i)
            for w in workers:
                q.put("STOP")

            for w in workers:
                w.start()

            for w in workers:
                w.join()

            results.put("STOP")
            for i in iter(results.get, "STOP"):
                results_list.append(i)

        return tuple(sorted(results_list, key=lambda x: x[0]))


def trm_indices(a, b):
    """
    Generate (i,j) indeces of a triangular matrix, between elements a and b.
    The matrix size is automatically determined from the number of elements.
    For instance: trm_indices((0,0),(2,1)) yields (0,0) (1,0) (1,1) (2,0)
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

    for i in range(1, n):
        for j in range(i):
            yield (i, j)


def trm_indices_diag(n):
    """generate (i,j) indeces of a triangular matrix of n rows (or columns),
    with diagonal

    Parameters
    ----------

    n : int
        Matrix size
    """

    for i in range(0, n):
        for j in range(i + 1):
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
        np.concatenate(
            tuple([e.trajectory.timeseries(order="fac") for e in universes]),
            axis=0,
        ),
        format=MemoryReader,
    )
