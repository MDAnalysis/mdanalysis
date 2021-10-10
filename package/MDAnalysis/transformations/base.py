# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
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

"""\
Transformations Base Class --- :mod:`MDAnalysis.transformations.base`
=====================================================================

.. autoclass:: TransformationBase

"""
from threadpoolctl import threadpool_limits


class TransformationBase(object):
    """Base class for defining on-the-fly transformations

    The class is designed as a template for creating on-the-fly
    Transformation classes. This class will

    1) set up a context manager framework on limiting the threads
    per call, which may be the multi-thread OpenBlas backend of NumPy.
    This backend may kill the performance when subscribing hyperthread
    or oversubscribing the threads when used together with other parallel
    engines e.g. Dask.
    (PR `#2950 <https://github.com/MDAnalysis/mdanalysis/pull/2950>`_)

    Define ``max_threads=1`` when that is the case.

    2) set up a boolean attribute `parallelizable` for checking if the
    transformation can be applied in a **split-apply-combine** parallelism.
    For example, the :class:`~MDAnalysis.transformations.positionaveraging.PositionAverager`
    is history-dependent and can not be used in parallel analysis natively.
    (Issue `#2996 <https://github.com/MDAnalysis/mdanalysis/issues/2996>`_)

    To define a new Transformation, :class:`TransformationBase`
    has to be subclassed.
    ``max_threads`` will be set to ``None`` by default,
    i.e. does not do anything and any settings in the environment such as
    the environment variable :envvar:`OMP_NUM_THREADS`
    (see the `OpenMP specification for OMP_NUM_THREADS <https://www.openmp.org/spec-html/5.0/openmpse50.html>`_)
    are used.
    ``parallelizable`` will be set to ``True`` by default.
    You may need to double check if it can be used in parallel analysis;
    if not, override the value to ``False``.
    Note this attribute is not checked anywhere in MDAnalysis yet.
    Developers of the parallel analysis have to check it in their own code.

    .. code-block:: python

       class NewTransformation(TransformationBase):
           def __init__(self, ag, parameter,
                        max_threads=1, parallelizable=True):
               super().__init__(max_threads=max_threads,
                                 parallelizable=parallelizable)
               self.ag = ag
               self._param = parameter

           def _transform(self, ts):
               #  REQUIRED
               ts.positions = some_function(ts, self.ag, self._param)
               return ts

    Afterwards the new transformation can be run like this.

    .. code-block:: python

       new_transformation = NewTransformation(ag, param)
       u.trajectory.add_transformations(new_transformation)


    .. versionadded:: 2.0.0
       Add the base class for all transformations to limit threads and
       check if it can be used in parallel analysis.
    """

    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        max_threads: int, optional
           The maximum thread number can be used.
           Default is ``None``, which means the default or the external setting.
        parallelizable: bool, optional
           A check for if this can be used in split-apply-combine parallel
           analysis approach.
           Default is ``True``.
        """
        self.max_threads = kwargs.pop('max_threads', None)
        self.parallelizable = kwargs.pop('parallelizable', True)

    def __call__(self, ts):
        """The function that makes transformation can be called as a function

        The thread limit works as a context manager with given `max_threads`
        wrapping the real :func:`_transform` function
        """
        with threadpool_limits(self.max_threads):
            return self._transform(ts)

    def _transform(self, ts):
        """Transform the given `Timestep`

        It deals with the transformation of a single `Timestep`.
        """
        raise NotImplementedError("Only implemented in child classes")
