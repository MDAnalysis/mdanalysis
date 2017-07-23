The unit tests and the test data are bundled together in the package **MDAnalyisTests-_release_**. In order to run the tests, this package must be installed in addition to MDAnalysis.

The tests also rely on the `pytest` and `numpy` packages, and require both to run.

# Users
Install [MDAnalysisTests](MDAnalysisTests) via
```
pip install --upgrade MDAnalysisTests
```
or download the tar file, unpack, and run `python setup.py install`.

Run the tests by invoking
```
python -c 'from MDAnalysis.tests import test; test()'
```

All tests should pass (i.e. no **FAIL**, **ERROR**); *SKIPPED* or *XFAIL* are ok. For anything that fails or gives an error [ask on the user mailing list](http://users.mdanalysis.org) or [raise an issue](/MDAnalysis/mdanalysis/issues).

# Developers #

All tests should pass (i.e. no **FAIL**, **ERROR**); *SKIPPED* or *XFAIL* are ok. For anything that fails or gives an error **fix your code** (or *raise an issue*).

**Do not push code that fails to the development branch.** Instead, push it to a feature or issue branch. It will run automatically through the unit tests by travis-ci and it will be available for comment and discussion by other developers.

## Recommended ##
Use the tests from the [git source repository](Source), which are located in the [testsuite/MDAnalysisTests](https://github.com/MDAnalysis/mdanalysis/tree/develop/testsuite) directory:
```
cd testsuite/MDAnalysisTests
pytest  --disable-pytest-warnings --numprocesses 4
```
(Try increasing the number of processes; with 24 processes on 12 cores (+hyperthreading) this took ~40 seconds; in serial it takes ~30 min).

To run specific tests just specify the path to the test file:
```
pytest test_analysis.py
```

Specific test classes inside test files, and even specific test methods, can also be specified:
```
# Test the entire TestContactMatrix class
pytest test_analysis.py::TestContactMatrix
# Test only test_sparse in the TestContactMatrix class
pytest test_analysis.py::TestContactMatrix::test_sparse
```

## Alternatives
You can install `MDAnalysisTests` and then run the tests anywhere. Extra functionality, afforded by our nose plugins, is added only if the tests are run through the `mda_nosetests` script, or by directly invoking `MDAnalysis.tests.test()` (which is what the 3-line `mda_nosetests` script does under the hood).

You can run all tests in serial from the python interpreter like this:
```python
from MDAnalysis.tests import test
test()
```

Fore more details see below.

## Examples ##
Examples for output in various modes. For debugging, verbose mode is more useful as one can identify failing tests while they are running.

### Serial testing ###
For example, a successful test might look like the following
```
>>> MDAnalysis.tests.test()
......S...S............................................................................................................................................................K.KK...........................................................................................................................................................................K...............................................................................................................................................................................................................K..............................................................K............................................................................................................................................................................................................................................................................................................................................................................K......K....................K.....................K...................K....................K...........................K...................K...................K.............................................................K...................K........................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................
============ 5497 passed, 31 skipped, 1 xfailed in 1056.77 seconds =============
```

### Coverage ###
We test code coverage of the unit tests with the  [coverage](http://nedbatchelder.com/code/modules/rees-coverage.html) plugin of nose. Currently, this is done automatically as part of the Travis job on Python 2.7 and is viewable on **[coveralls](https://coveralls.io/r/MDAnalysis/mdanalysis?branch=develop)**.

If you want to generate a coverage report manually you can run
```
cd testsuite
pytest --numprocesses 4 --cov=MDAnalysis --cov-report html
```

# Details #
We are borrowing some of NumPy's testing frame work; thus, numpy **must** be installed for the tests to run at all. The tests require at least numpy 1.5.

## Running tests from within python ##
Run all the tests with
```
   import MDAnalysis.tests
   MDAnalysis.tests.test()
```

## Running tests from the command line ##
Instead of running tests from within python, one can also run them via the `mda_nosetests` script that ships with MDAnalysisTests. With version 0.11.0 the test subsystem was overhauled to allow the incorporation of customized nose plugins. In order for them to work tests must be invoked via our own wrapper function `MDAnalysis.tests.test()`. This is what the `mda_nosetests` script does for you. Alternatively, you can call `MDAnalysis.tests.test()` from the interpreter. `mda_nosetests` strives to be compatible and interchangeable with nose's [nosetests](http://somethingaboutorange.com/mrl/projects/nose/0.11.2/usage.html) script, with [added functionality](#compatibility-with-nosetests).

Go into the tests directory (or the package root)
```
cd /testsuite/MDAnalysisTests
```
and invoke ./mda_nosetests directly to run **all tests** on two processors in parallel ("`%`" is the shell prompt and should _not_ be typed):
```
% pytest --numprocesses 4
```

When you have written a **new unit test** it is helpful to check that it passes without running the entire suite. For example, in order to test everything in, say, `test\_selections.py` run
```
% pytest test_selections
..............
----------------------------------------------------------------------
Ran 14 tests in 3.421s

OK
```
One can also test individual test classes. For instance, after working on the XYZReader one can check just the TestCompressedXYZReader tests with
```
% pytest test_coordinates::TestCompressedXYZReader
....
----------------------------------------------------------------------
Ran 4 tests in 0.486s

OK
```
where we are testing the class `TestCompressedXYZReader` which can be found in the module (file) `test\_coordinates.py`.


### Running tests with setuptools ###

Setuptools can also use [nose](http://somethingaboutorange.com/mrl/projects/nose) directly (and it takes care of having all the libraries in place):
```
python setup.py nosetests
```

If you have the [coverage](http://nedbatchelder.com/code/modules/rees-coverage.html) package installed, you can also check code coverage of the tests:
```
python setup.py nosetests --with-coverage --cover-package=MDAnalysis --cover-erase --cover-tests
```

# TODO: Fix further

## Compatibility with numpy's test wrapper ##

Up to MDAnalysisTests version 0.11.0 numpy's wrapper was used to run tests when invoked through `MDAnalysis.tests.test()` (but not through `nosetests`). This is now replaced by our own wrapper.

Main differences are that numpy-specific arguments to `MDAnalysis.tests.test()` are now either emulated or simply not implemented. Here's a list of the ones most commonly used with MDAnalysisTests:

* `label`: this allows the selection of tests based on whether or not they were decorated with the `@dec.slow` decorator. Current behavior recognizes only `label='fast'`. Any other label (or its absence) defaults to running all tests.
* `verbose`: this argument is no longer accepted. Pass one of `-v`, `--verbose`, or `--verbosity=n` in `argv`;
* `extra_argv`: this argument allows an extra list of flags to be passed to nose. It is still accepted, but deprecated in favor of nose's identical `argv` argument.

Any other numpy-specific arguments will not be accepted and will cause the test run to fail.

Under numpy the behavior was not to silence any test stdout output. The behavior now is to silence it; this can be reversed with the `-s` flag.

## Data ##
The simulation data used in tests are all released under the same license as MDAnalysis or are in the Public Domain (such as PDBs from the Protein Databank). An incomplete list of sources:
* from Beckstein et al. (2009) (`adk.psf`,`adk_dims.dcd`)
  * _adk\_dims_      Trajectory of a macromolecular transition of the enzyme adenylate kinase between a closed and an open conformation. The simulation was run in [CHARMM](http://www.charmm.org) c35a1.
* unpublished simulations (O. Beckstein)
  * _adk\_oplsaa_    Ten frames from the first 1 ns of a equilibrium trajectory of AdK in water with Na+ counter ions. The OPLS/AA forcefield is used with the TIP4P water model. The simulation was run with [Gromacs](http://www.gromacs.org) 4.0.2.
* contributions from developers and users
* Protein Databank

## Purposely failing test ##

It is sometime convenient to have a failing test. This is helpful, for instance, when testing a plugin for nosetests, or to make sure a script catches failure as expected. The `test_failure` test in `test_failure.py` can be made to fail by setting the `MDA_FAILURE_TEST` environment variable.

### References ###

  * O. Beckstein, E.J. Denning, J.R. Perilla and T.B. Woolf, Zipping and Unzipping of Adenylate Kinase: Atomistic Insights into the Ensemble of Open-Closed Transitions. J Mol Biol 394 (2009), 160--176, doi:[10.1016/j.jmb.2009.09.009](http://dx.doi.org/10.1016/j.jmb.2009.09.009)


# Writing test cases #

The tests are in a separate package, together with any data files required for running the tests (see [Issue 87](http://issues.mdanalysis.org/87) for details). Whenever you _add a new feature_ to the code you _should also add a test case_ (ideally, in the same git commit so that the code and the test case are treated as one unit).

The unit tests use the [unittest module](http://docs.python.org/library/unittest.html) together with [nose](http://somethingaboutorange.com/mrl/projects/nose/0.11.3/index.html). See the examples in the [MDAnalysisTests](https://github.com/MDAnalysis/mdanalysis/tree/develop/testsuite/MDAnalysisTests) package.

The [SciPy testing guidelines](http://projects.scipy.org/numpy/wiki/TestingGuidelines#id11) are a good howto for writing test cases.

Conventions for MDAnalysis
  * Relative import statements are now banned from unit testing modules (see [Issue #189](/MDAnalysis/mdanalysis/issues/189) for details)
  * using `os.chdir()` is banned because it can break the tests in really weird ways (see [Issue #556](https://github.com/MDAnalysis/mdanalysis/issues/556)): use `with tempdir.in_tempdir()` or something similar
  * Test input data is stored in  [MDAnalysisTests/data](https://github.com/MDAnalysis/mdanalysis/tree/develop/testsuite/MDAnalysisTests/data).
    * Keep files small if possible; for trajectories 10 frames or less are sufficient.
    * Add the file name of test data files to [MDAnalysisTests/datafiles.py](https://github.com/MDAnalysis/mdanalysis/blob/develop/testsuite/MDAnalysisTests/datafiles.py) (see the code for details).
    * Add the file(s) or a glob pattern to the `package_data` in [setup.py](https://github.com/MDAnalysis/mdanalysis/blob/develop/testsuite/setup.py); otherwise the file will not be included in the python package.
    * If you use data from a published paper then add a reference to _this wiki page_ and the doc string in [MDAnalysisTests/\_\_init\_\_.py](https://github.com/MDAnalysis/mdanalysis/blob/develop/testsuite/MDAnalysisTests/__init__.py).
  * Tests are currently organized by top-level module. Each file containing tests must start with `test_` by convention (this is how nose/unittest works). Tests itself also have to follow the appropriate naming conventions. See the docs above or the source.
  * Tests that take longer than 3 seconds to run should be marked `@slow` (see e.g. the XTC tests in [MDAnalysisTests/test\_coordinates.py](https://github.com/MDAnalysis/mdanalysis/blob/develop/testsuite/MDAnalysisTests/test_coordinates.py). They will only be run if `labels="full"` is given as an argument to the `test()` function.
  * Add a test for
    * new functionality
    * fixed issues (typically named `test_IssueXX` or referencing the issue in the doc string (to avoid regression)
    * anything you think worthwhile – the more the better!


# Changes with releases #

The way we organized the unit tests changed between releases. The procedure for the current release is detailed at the very top of the page. The following list is for historical reference and in case you ever want to go back to a previous release.

  1. since **0.11.0**: the testing subsystem was overhauled to allow the use of plugins external to nose. We also no longer use numpy's `test()` wrapper. `mda_nosetests` is now the preferred way to run the tests from the command-line in a mostly backward-compatible way with the usage of `nosetests`. Most numpy-specific arguments to `test()` are now deprecated in favor of nose flags.
  1. since **0.7.5**: tests _and_ data are together in package **MDAnalysisTests**. See [Issue 87](http://issues.mdanalysis.org/87) for details.
  1. release **0.7.4**: tests are in **MDAnalysis** and data is in **MDAnalysisTestData** (for MDAnalysis == 0.7.4). To install [MDAnalysisTestData](MDAnalysisTestData) download the `MDAnalysisTestData-0.7.4.tar.gz` from the [Download](http://code.google.com/p/mdanalysis/downloads/list) section or try ```easy_install http://mdanalysis.googlecode.com/files/MDAnalysisTestData-0.7.4.tar.gz```
  1. release **0.6.1** to **0.7.3**: tests and data were included with **MDAnalysis**
  1. release **0.4** to **0.6.0**: no tests included