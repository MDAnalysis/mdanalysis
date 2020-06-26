.. -*- coding: utf-8 -*-
.. _datasets:

============
Example data
============

MDAnalysis offers a collection of example data files and datasets to run tests and tutorials. These are split into two packages:

    * MDAnalysisTests: primarily for unit tests of the code
    * MDAnalysisData: datasets for workshops and tutorials

.. _mdanalysistests:

MDAnalysisTests
===============

While this is installed as a separate package, you should import files like so::

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PSF, DCD

    u = mda.Universe(PSF, DCD)

A complete list of files and descriptions is in the ``mdanalysis/testsuite/MDAnalysisTests/datafiles.py`` file. The actual files are stored in the ``mdanalysis/testsuite/MDAnalysisTests/data/`` directory. 

.. _mdanalysisdata:

MDAnalysisData
==============

The `MDAnalysisData <https://www.mdanalysis.org/MDAnalysisData/>`_ package is an interface to download, cache, and access certain datasets hosted on external repositories (e.g. figshare_, zenodo_, DataDryad_). Data is not downloaded upon installation, so the package itself is small; but the directory where the datasets are cached can grow significantly.

You can access datasets like so::

    import MDAnalysis as mda
    from MDAnalysisData import datasets
    adk = datasets.fetch_adk_equilibrium()
    u = mda.Universe(adk.topology, adk.trajectory)

    # to see the description of the dataset
    print(adk.DESCR)

The cached files are stored by default in the ``~/MDAnalysis_data`` directory. This can be changed by setting the environment variable ``MDANALYSIS_DATA``. You can change this for a Terminal session with the below command:

.. code-block:: bash

    export MDANALYSIS_DATA=/my/chosen/path/MDAnalysis_data

Add it to your ``.bashrc`` for a permanent change.

.. code-block:: bash

    echo 'export MDANALYSIS_DATA=/my/chosen/path/MDAnalysis_data' >> ~/.bashrc

In Python, you can check the location of your caching directory::

    MDAnalysisData.base.get_data_home()

And clear the directory with::

    MDAnalysisData.base.clear_data_home()

`A list of datasets <https://www.mdanalysis.org/MDAnalysisData/usage.html>`_ can be found at the MDAnalysisData documentation.


.. _figshare: https://figshare.com/
.. _zenodo: https://zenodo.org/
.. _DataDryad: https://www.datadryad.org/