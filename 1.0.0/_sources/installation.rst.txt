.. -*- coding: utf-8 -*-

====================
Installation
====================

The latest versions of **MDAnalysis** can be installed using `conda` or `pip`. 
Currently, the conda releases only support serial calculations.
If you plan to use the parallel OpenMP algorithms, you need to 
install MDAnalysis with pip and have a working OpenMP installation.

MDAnalysis has a separate :ref:`test suite <mdanalysistests>` **MDAnalysisTests** that is required to run the test cases and examples. 
The test files change less frequently, take up around 90 MB of space, 
and are not needed for daily use of MDAnalysis. However, they are often used in examples,
including many in this User Guide. If you are not interested in developing 
MDAnalysis or using the example files, you most likely don't need the tests. If you want to 
run examples in the User Guide, install the tests. 
The tests are distributed separately from the main package. 

conda
=====
To install the latest stable version of MDAnalysis via ``conda``, use the following command. This installs all dependencies needed for full analysis functionality (excluding external programs such as `HOLE`_):

.. code-block:: bash

    conda install -c conda-forge mdanalysis

To upgrade:

.. code-block:: bash

    conda update mdanalysis

To install the tests:

.. code-block:: bash

    conda install -c conda-forge MDAnalysisTests

pip
=====
The following command will install or upgrade the latest stable version of MDAnalysis via ``pip``, with core dependencies. This means that some packages required by specific analysis modules will not be installed.

.. code-block:: bash

    pip install --upgrade MDAnalysis

If you need to install a fully-featured MDAnalysis, add the ``analysis`` tag. As with ``conda``, this will not install external programs such as `HOLE`_.

.. code-block:: bash

    pip install --upgrade MDAnalysis[analysis]

To install/upgrade tests:

.. code-block:: bash

    pip install --upgrade MDAnalysisTests

Development versions
====================
To install development versions of MDAnalysis, you can compile it from source.

.. code-block:: bash

    git clone https://github.com/MDAnalysis/mdanalysis
    cd mdanalysis
    pip install -e .

In order to install from source, you will need ``numpy`` and ``cython``. See :ref:`create-virtual-environment` for instructions on how to create a full development environment.

Testing
-------

The tests rely on the `pytest` and `numpy` packages, which must also be installed. Run tests with: 

.. code-block:: bash

    pytest --disable-pytest-warnings --pyargs MDAnalysisTests

All tests should pass (i.e. no FAIL, ERROR); SKIPPED or XFAIL are ok. If anything fails or gives an error, 
`ask on the user mailing list <http://users.mdanalysis.org/>`_ or `raise an issue <https://github.com/MDAnalysis/mdanalysis/issues>`_.

Testing MDAnalysis can take a while, as there are quite a few tests. 
The plugin `pytest-xdist <https://github.com/pytest-dev/pytest-xdist>`_ can be used to run tests in parallel.

.. code-block:: bash

    pip install pytest-xdist
    pytest --disable-pytest-warnings --pyargs MDAnalysisTests --numprocesses 4


Additional datasets
===================

:ref:`MDAnalysisData is an additional package <mdanalysisdata>` with datasets that can be used in example tutorials. You can install it with ``conda`` or ``pip``:

.. code-block:: bash

    # conda
    conda install -c conda-forge mdanalysisdata
    # pip
    pip install --upgrade MDAnalysisData

This installation does not download all the datasets; instead, the datasets are cached when they are first downloaded using a Python command. 


.. _`HOLE`: http://www.holeprogram.org
