=====================
MDAnalysis benchmarks
=====================

Benchmarking MDAnalysis with Airspeed Velocity.

Usage with MDAnalysis
---------------------

Airspeed Velocity builds clean conda environments to
benchmark the performance of MDAnalysis at different
time points in its history (for different git commit
hashes).

To build / install Airspeed Velocity it should
suffice to clone the `git repo`_, building the develop
branch with::

    python setup.py install --user

`Airspeed Velocity commands`_ are described in detail in their
documentation. A common usage example for evaluating the
performance of a feature branch / pull request would be::

    asv continuous --bench GROReadBench d76c9348 e0bc303 -e

In the above, ``GROReadBench`` is handled as a regular
expression for the specific benchmark test to run between
the provided git commit hashes.

To evaluate a benchmark test over the course of project
history one would commonly use ``asv run``. For example,
to probe performance for trajectory readers at 20 commit
hashes evenly spread over a few years of the project one
might run::

     asv run -e -s 20 ddb57592..e0bc3034 --bench TrajReaderBench -j

It is also possible to specify ``ALL`` to space the performance
tests over the entire lifetime of the project, but exercise
caution as very early commits may represent a state of the
project where many features are not available and / or
files are not in the expected locations. Using ``--merges`` is also
frequently advisable as merge commits are more likely to build
and run successfully.

The ``asv run`` command will store detailed benchmark data locally
as ``JSON`` files, which can be converted into interactive website
data and hosted locally with::

    asv publish
    asv preview

.. _git repo: https://github.com/airspeed-velocity/asv
.. _Airspeed Velocity commands: http://asv.readthedocs.io/en/latest/commands.html

Writing benchmarks
------------------

The Airspeed Velocity `documentation for writing benchmarks`_ is a
suitable reference. As a quick summary of guidelines:

- wrap imports from MDAnalysis in the test modules because older
  commit hashes may not have the name imported for various reasons::

     try:
        from MDAnalysis.coordinates.DCD import DCDReader
     except ImportError:
        pass

  The benchmarks themselves will automatically handle the missing
  features if the above is done.

- leave the timing code to ASV -- don't implement your own

- the benchmarks are written as Python modules in the `benchmark`
  directory (and subdirectories thereof). There are no formal
  naming requirements for these modules. Benchmarks are generally
  written as functions of the form `time_feature()` with the `time`
  prefix, sometimes within a broader class object. A `setup` method
  or attribute may be used to perform operations that are required
  for the test suite, but should not be included in the performance
  timing.

- parametrized benchmarks can be quite powerful for testing several
  inputs to a test; note that all possible combinations will be tested,
  and it is very useful to label the parameters with names as these
  will be nicely summarized in the output::

       params = (['XTC', 'TRR', 'DCD'],
                 [10, 20, 30])
       param_names = ['traj_format', 'num_atoms']

- memory usage can also be profiled with test prefixes including `mem`
  and `peakmem`

.. _documentation for writing benchmarks: http://asv.readthedocs.io/en/latest/writing_benchmarks.html

Advanced Notes
--------------

- the depedencies installed in the clean conda benchmarking environments,
  and indeed the decision to use conda over virtualenv, can be controlled
  in the ``asv.conf.json`` file, as can which versions of Python are probed

- the above file also controls which branch of MDAnalysis is used for a
  first-pass check of the benchmarks that are written--regardless of where you
  run the benchmarks from, the current ``JSON`` file indicates that ASV
  will check your benchmark code against the latest commit hash on develop
  branch first, before running the actual benchmarks for the specified commit
  hashes
