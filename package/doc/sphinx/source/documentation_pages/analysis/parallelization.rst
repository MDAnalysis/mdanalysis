.. -*- coding: utf-8 -*-

=================
Parallel analysis
=================

Starting v.2.8.0, MDAnalysis has introduced process-based parallelization for
virtually all analysis classes. This page will briefly explain what has changed
in the :class:`mdanalysis.analysis.base.AnalysisBase` protocol, how it affects
users and developers, when you should use parallelization (almost always!), and
when you can abstain from doing so (rarely).

All this work was done within Google Summer of Code summer project by @marinegor
and MDAnalysis core developers team.


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
parallelization to a class is very simple, it won't take much time untill it's
introduced to all trivial classes.


How does parallelization work
=============================

The main idea behind its current version is that a trajectory analysis is almost
always trivially parallelizable, meaning you can analyze all frames
independently, and then merge them in a single object. This approach is also
known as "split-apply-combine", and isn't new to MDAnalysis users, since it was
first introduced in [pmda](https://github.com/mdanalysis/pmda). Version 2.8.0 of
MDAnalysis brings this approach to the main library.


split-apply-combine
-------------------

The following scheme explains the current ``AnalysisBase.run()`` protocol
(user-implemented methods are highlighted in orange):

.. figure:: /images/AnalysisBase_parallel.png


In short, after checking input parameters and configuring backend,
``AnalysisBase`` splits all the frames into *computation_groups* (equally sized
sequencial groups of frames to be processed by each worker). All groups then get
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

From a user point of new, there are no new (non-trivial) methods. Meaning, if
you want to write your own analysis class, you still have to implement only
`_prepare`, `_single_frame` and `_conclude`. However, from a developer point of
view, there are few:

#. :meth:`MDAnalysis.analysis.base.AnalysisBase._define_run_frames`
#. :meth:`MDAnalysis.analysis.base.AnalysisBase._prepare_sliced_trajectory`
#. :meth:`MDAnalysis.analysis.base.AnalysisBase._configure_backend`
#. :meth:`MDAnalysis.analysis.base.AnalysisBase._setup_computation_groups`
#. :meth:`MDAnalysis.analysis.base.AnalysisBase._compute`
#. :meth:`MDAnalysis.analysis.base.AnalysisBase._get_aggregator`

First two methods share the functionality of ``_setup_frames``.
``_define_run_frames`` is run once during analysis, as it checks that input
parameters ``start, stop, step`` or ``frames`` are consistent with the given
trajectory and prepares the ``slicer`` object that defines the iteration pattern
through trajectory. ``_prepare_sliced_trajectory`` assigns
``self._sliced_trajectory`` attribute, and also number of frames in it, and
``self.frames`` and ``self.times`` arrays. In case the computation will be later
split between other processes, this method will be called again on each of the
computation groups.

Method ``_configure_backend`` performs basic health checks for a given analysis
class -- namely, it compares given backend (if it's a ``str`` instance, such as
``'multiprocessing'``) with the list of builtin backends (and also backends
implemented for a given analysis subclass), and configures a
:class:`MDAnalysis.analysis.backends.BackendBase` instance accordingly. If user
decides to provide a custom backend (any subclass of
:class:`MDAnalysis.analysis.backends.BackendBase`, or anything with ``apply``
method), ensures that number of workers wasn't specified twice (during backend
initialization and in ``run()`` arguments).

After backend is configured, ``_setup_computation_groups`` splits the frames
prepared earlier in ``self._prepare_sliced_trajectory`` into number of groups,
by default equal to the number of workers. 

In the ``_compute`` method, frames get initialized again with
``_prepare_sliced_trajectory``, and also attributes necessary for a specific
analysis get initialized with `_prepare` method. Then the function iterates over
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


When to use parallelization *aka* known limitations
===================================================

Fast `_single_frame` compared to reading from disk
--------------------------------------------------

Seriralization issues
---------------------

Out of memory issues
--------------------

Progress bar is missing
-----------------------



Adding parallelization to your own analysis class
=================================================

Which methods need to be implemented, and how to make sure of that -- see RMSD as example.

