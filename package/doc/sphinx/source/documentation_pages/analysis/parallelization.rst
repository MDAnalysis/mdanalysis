.. -*- coding: utf-8 -*-

=================
Parallel analysis
=================

Starting from v2.8.0, MDAnalysis has introduced process-based parallelization for
virtually all analysis classes. This page will briefly explain what has changed
in the :class:`MDanalysis.analysis.base.AnalysisBase` protocol, how it affects
users and developers, when you should use parallelization (almost always!), and
when you should abstain from doing so (rarely).

All this work was done as part of Google Summer of Code 2023 by @marinegor
and MDAnalysis GSoC mentors.


How to use parallelization
==========================

In order to use parallelization in a built-in analysis class ``SomeClass``, simply check which
backends are available, and then just enable them by providing
``backend='multiprocessing'`` and ``n_workers=...`` to ``SomeClass.run(...)`` method:

.. code-block:: python

    u = mda.Universe(...)
    my_run = SomeClass(trajectory)
    assert SomeClass.get_supported_backends() == ('serial', 'multiprocessing', 'dask')

    my_run.run(backend='multiprocessing', n_workers=12)

For some classes, such as :class:`MDAnalysis.analysis.rms.RMSF`,
split-apply-combine parallelization isn't possible, and running them will be
impossible with any but ``serial`` backend.

Also, backends are getting added to the new classes slowly -- in 2.8.0 only
:class:`MDAnalysis.analysis.rms.RMSD` got new backends. But since adding
parallelization to a class is very simple, it won't take much time until it's
introduced to all trivial classes.


How does parallelization work
=============================

The main idea behind its current implementation is that a trajectory analysis is
often trivially parallelizable, meaning you can analyze all frames
independently, and then merge them in a single object. This approach is also
known as "split-apply-combine", and isn't new to MDAnalysis users, since it was
first introduced in [pmda](https://github.com/mdanalysis/pmda). Version 2.8.0 of
MDAnalysis brings this approach to the main library.


split-apply-combine
-------------------

The following scheme explains the current ``AnalysisBase.run()`` protocol
(user-implemented methods are highlighted in orange):

.. figure:: /images/AnalysisBase_parallel.png


In short, after checking input parameters and configuring the backend,
``AnalysisBase`` splits all the frames into *computation groups* (equally sized
sequential groups of frames to be processed by each worker). All groups then get
**split** between workers of a backend configured early, the main instance gets
serialized and distributed between workers, and then
:meth:`MDAnalysis.analysis.AnalysisBase._compute()` method gets called for all
frames of a computation group. Within this method, a user-implemented
``_single_frame`` method gets **applied** to each frame in a computation group.
After that, the main instance gets an object that will **combine** all the
objects from remote workers, and all instances get *merged* with an instance of
:class:`MDAnalysis.analysis.results.ResultsGroup`. Then, a normal
user-implemented ``_compute`` method is called.

Since one of the goals of the project was to **not** break any existing code,
the protocol mimics the single-process parallelization where possible. Meaning,
user-implemented methods such as ``_prepare``, ``_single_frame`` and
``_conclude`` won't need to know they are operating on an instance within the
main python process, or on a remote instance, since the executed code is the
same in both cases.


New methods in ``AnalysisBase``
-------------------------------

From a user point of new, there are no new (non-trivial) methods. If
you want to write your own analysis class, you still have to implement only
``_prepare``, ``_single_frame`` and ``_conclude``. However, from a developer point of
view, there are a few new methods:

#. :meth:`MDAnalysis.analysis.base.AnalysisBase._define_run_frames`
#. :meth:`MDAnalysis.analysis.base.AnalysisBase._prepare_sliced_trajectory`
#. :meth:`MDAnalysis.analysis.base.AnalysisBase._configure_backend`
#. :meth:`MDAnalysis.analysis.base.AnalysisBase._setup_computation_groups`
#. :meth:`MDAnalysis.analysis.base.AnalysisBase._compute`
#. :meth:`MDAnalysis.analysis.base.AnalysisBase._get_aggregator`

The first two methods share the functionality of ``_setup_frames``.
``_define_run_frames`` is run once during analysis, as it checks that input
parameters ``start, stop, step`` or ``frames`` are consistent with the given
trajectory and prepares the ``slicer`` object that defines the iteration pattern
through the trajectory. ``_prepare_sliced_trajectory`` assigns
``self._sliced_trajectory`` attribute, and also number of frames in it, and
``self.frames`` and ``self.times`` arrays. In case the computation will be later
split between other processes, this method will be called again on each of the
computation groups.

Method ``_configure_backend`` performs basic health checks for a given analysis
class -- namely, it compares a given backend (if it's a ``str`` instance, such as
``'multiprocessing'``) with the list of builtin backends (and also backends
implemented for a given analysis subclass), and configures a
:class:`MDAnalysis.analysis.backends.BackendBase` instance accordingly. If the user
decides to provide a custom backend (any subclass of
:class:`MDAnalysis.analysis.backends.BackendBase`, or anything with ``apply``
method), it ensures that number of workers wasn't specified twice (during backend
initialization and in ``run()`` arguments).

After a backend is configured, ``_setup_computation_groups`` splits the frames
prepared earlier in ``self._prepare_sliced_trajectory`` into a number of groups,
by default equal to the number of workers. 

In the ``_compute`` method, frames get initialized again with
``_prepare_sliced_trajectory``, and attributes necessary for a specific
analysis get initialized with ``_prepare`` method. Then the function iterates over
``self._sliced_trajectory``, assigning ``self._frame_index`` and ``self._ts`` as
frame index (within a computation group) and timestamp, and also setting
respective ``self.frames`` and ``self.times`` array values.

After ``_compute`` has finished, the main analysis instance calls
``_get_aggretator`` method, which merges the ``self.results`` attributes from
other processes into a single :class:`MDAnalysis.analysis.results.Results`
instance, making it look for the subsequent ``_conclude`` method like the run
was performed in a serial fashion, without parallelization.


New classes: ``ResultsGroup`` and ``BackendBase``
=================================================

``ResultsGroup``
----------------

:class:`MDAnalysis.analysis.results.ResultsGroup` extends the functionality of
the :class:`MDAnalysis.analysis.results.Results` class. Since the ``Results``
class is basically a dictionary that also keeps track of assigned attributes, it
is possible to iterate over all these attributes later. ``ResultsGroup`` does
exactly that: given a list of the ``Results`` objects with the same attributes,
it applies a specific aggregation function to every attribute, and stores it as
a same attribute of the returned object:

.. code-block:: python

    from MDAnalysis.analysis.results import ResultsGroup, Results
    group = ResultsGroup(lookup={'mass': ResultsGroup.float_mean})
    obj1 = Results(mass=1)
    obj2 = Results(mass=3)
    assert group.merge([obj1, obj2]) == Results(mass=2.0)


``BackendBase``
---------------

:class:`MDAnalysis.analysis.backends.BackendBase` holds all backend attributes,
and also implements an :meth:`MDAnalysis.analysis.backends.BackendBase.apply`
method, applying a given function to a list of its parameters, but in a parallel
fashion. Although in ``AnalysisBase`` it is used to apply a ``_compute``
function, in principle it can be used to any arbitrary function and arguments,
given they're serializable.


When to use parallelization? (Known limitations)
================================================

For now, the syntax for running parallel analysis is explicit, meaning by
default the ``serial`` version will be run, and the parallelization won't be
enabled by default. Although we expect the parallelization to be useful in most
cases, there are some known caveats from the inital benchmarks.

Fast ``_single_frame`` compared to reading from disk
--------------------------------------------------

In all cases, parallelization will not be useful only when frames are being
processed faster than being read from the disk, otherwise reading is the
bottleneck here. Hence, you'll benefit from parallelization only if you have
relatively much compute per frame, or a fast drive, as illustrated below:

.. figure:: /images/parallelization_time.png

In other words, if you have *fast* analysis (say,
:class:`MDAnalysis.analysis.rms.RMSD`) **and** a slow HDD drive, you are likely
to not get any benefits from parallelization. Otherwise, you should be fine.

Serialization issues
--------------------

For built-in analysis classes, the default serialization with both
``multiprocessing`` and ``dask`` is known to work. If you're using some custom
analysis class that e.g. stores a non-serializable object in one of its
attributes, you might get a serialization error (``PicklingError`` if you're
using a ``multiprocessing`` backend). If you want to get around that, we suggest
trying ``backend='dask'`` (it uses ``dask`` serialization engine instead of
``pickle``).

Out of memory issues
--------------------

If you have large memory footprint of each worker, you can run into
out-of-memory errors (i.e. your server freezes when executing a run). In this
case we suggest decreasing the number of workers from all available CPUs (that
you can get with ``multiprocessing.cpu_count()``) to a smaller number.

Progress bar is missing
-----------------------

It is yet not possible to get a progress bar running with any parallel backend.
If you want an ETA of your analysis, we suggest running it in ``serial`` mode
for the first 10-100 frames with ``verbose=True``, and then running it with
multiple workers. Processing time scales almost linearly, so you can get your
ETA by dividing ``serial`` ETA by the number of workers.


Adding parallelization to your own analysis class
=================================================

If you want to add parallelization to your own analysis class, first make sure
your algorithm allows you to do that, i.e. you can process each frame independently.
Then it's rather simple -- let's look at the actual code that added
parallelization to the :class:`MDAnalysis.analysis.rms.RMSD`:

.. code-block:: python

    from MDAnalysis.analysis.base import AnalysisBase
    from MDAnalysis.analysis.results import ResultsGroup

    class RMSD(BackendBase):
        @classmethod
        def get_supported_backends(cls):
            return ('serial', 'multiprocessing', 'dask',)

        @classmethod
        def is_parallelizable(self):
            return True
        
        def _get_aggregator(self):
            return ResultsGroup(lookup={'rmsd': ResultsGroup.ndarray_vstack})


That's it! First two methods are boilerplate -- ``get_supported_backends``
returns a tuple with built-in backends that will work for your class (if there
are no serialization issues, it should be all three), and ``is_parallelizable``
is ``True`` (which is set to ``False`` in ``AnalysisBase``, hence we have to
re-define it), and ``_get_aggregator`` will be used as described earlier. Note
that :mod:`MDAnalysis.analysis.results` also provides few convenient functions
(defined as class methods of ``ResultsGroup``) for results aggregation:

#. :meth:`MDAnalysis.analysis.results.ResultsGroup.flatten_sequence`
#. :meth:`MDAnalysis.analysis.results.ResultsGroup.ndarray_sum`
#. :meth:`MDAnalysis.analysis.results.ResultsGroup.ndarray_mean`
#. :meth:`MDAnalysis.analysis.results.ResultsGroup.float_mean`
#. :meth:`MDAnalysis.analysis.results.ResultsGroup.ndarray_hstack`
#. :meth:`MDAnalysis.analysis.results.ResultsGroup.ndarray_vstack`


So you'll likely find appropriate functions for basic aggregation there.

Writing custom backends
=======================

In order to write your custom backend (e.g. using ``dask.distributed``), inherit
the :class:`MDAnalysis.analysis.backends.BackendBase` and (re)-implement
``__init__`` and ``apply`` methods. Optionally, you can implement methods for
validation of correct backend initialization -- ``_get_checks`` and
``get_warnings``:

.. code-block:: python

    from MDAnalysis.analysis.backends import BackendBase
    class ThreadsBackend(BackendBase):
        def __init__(self, n_workers: int, starting_message: str = "Useless backend"):
            self.n_workers = n_workers
            self.starting_message = starting_message
            self._validate()

        def _get_warnings(self):
            return {True: 'warning: this backend is useless'}
        
        def _get_checks(self):
            return {isinstance(self.n_workers, int), 'error: self.n_workers is not an integer'}

        def apply(self, func, computations):
            from multiprocessing.dummy import Pool

            with Pool(processes=self.n_workers) as pool:
                print(self.starting_message)
                results = pool.map(func, computations)
            return results
    

In order to use a custom backend, you must explicitly state that you're about to use an unsupported_backend:

.. code-block:: python

    from MDAnalysis.analysis.rms import RMSD
    R = RMSD(...) # setup the run
    n_workers = 2
    backend = ThreadsBackend(n_workers=n_workers)
    R.run(backend=backend, unsupported_backend=True)