
.. _working-with-user-guide:

==============================
Contributing to the user guide
==============================

MDAnalysis maintains two kinds of documentation: 

    #. `This user guide <https://www.mdanalysis.org/UserGuide/>`__: a map of how MDAnalysis works, combined with tutorial-like overviews of specific topics (such as the analyses)
    
    #. `The documentation generated from the code itself <https://www.mdanalysis.org/docs/>`__. Largely built from code docstrings, these are meant to provide a clear explanation of the usage of individual classes and functions. They often include technical or historical information such as in which version the function was added, or deprecation notices.

This guide is about how to contribute to the user guide. If you are looking to add to documentation of the main code base, please see :ref:`working-with-mdanalysis-docs`.

The user guide makes use of a number of Sphinx extensions to ensure that the code examples are always up-to-date. These include `nbsphinx <https://nbsphinx.readthedocs.io>`_ and the `ipython directive <http://matplotlib.org/sampledoc/ipython_directive.html>`__.

The ``ipython`` directive lets you put code in the documentation which will be run
during the doc build. For example:

    ::

        .. ipython:: python

            x = 2
            x**3

will be rendered as:

    .. ipython::

        In [1]: x = 2

        In [2]: x**3
        Out[2]: 8

Many code examples in the docs are run during the
doc build. This approach means that code examples will always be up to date,
but it does make the doc building a bit more complex.

Here is an overview of the development workflow for the user guide, as expanded on throughout the rest of the page.

    #. :ref:`Fork the MDAnalysis repository <forking-user-guide>` from the mdanalysis account into your own account
    #. :ref:`Set up an isolated virtual environment <create-virtual-environment-user-guide>` for your documentation
    #. :ref:`Create a new branch off the master branch <create-code-branch>`
    #. Add your new documentation.
    #. :ref:`Build and view the documentation <build-user-guide>`.
    #. :ref:`Commit and push your changes, and open a pull request <add-changes-user-guide>`.


.. _forking-user-guide:

Forking and cloning the User Guide
==================================

Go to the `MDAnalysis project page <https://github.com/MDAnalysis/UserGuide>`_ and hit the :guilabel:`Fork` button. You will
want to clone your fork to your machine:

    .. code-block:: bash

        git clone https://github.com/your-user-name/UserGuide.git
        cd UserGuide
        git remote add upstream https://github.com/MDAnalysis/UserGuide

This creates the directory `UserGuide` and connects your repository to
the upstream (main project) MDAnalysis repository.


.. _create-virtual-environment-user-guide:

Creating a development environment
==================================

:ref:`Create a new virtual environment <create-virtual-environment>` for the user guide. Install the required dependencies, and activate the ``nglview`` extension. We use ``nglview`` for visualizing molecules in Jupyter notebook tutorials.

If using conda:

    .. code-block:: bash

        cd UserGuide/
        conda env create python=3.6 -f environment.yml --quiet
        conda activate mda-user-guide
        jupyter-nbextension enable nglview --py --sys-prefix


If using pip:

    .. code-block:: bash

        cd UserGuide/
        pip install -r requirements.txt
        jupyter-nbextension enable nglview --py --sys-prefix

.. _build-user-guide:

Building the user guide
=======================

Navigate to the ``doc/`` directory and run ``make html``:

    .. code-block:: bash

        cd doc/
        make html

The HTML output will be in ``doc/build/``, which you can open in your browser of choice. The homepage is ``doc/build/index.html``.

If rebuilding the documentation becomes tedious after a while, install the :ref:`sphinx-autobuild <autobuild-sphinx>` extension. 

Saving state in Jupyter notebooks
=================================

One of the neat things about ``nglview`` is the ability to interact with molecules via the viewer. This ability can be preserved for the HTML pages generated from Jupyer notebooks by ``nbsphinx``, if you save the notebook widget state after execution.

.. _add-changes-user-guide:

Adding changes to the UserGuide
===============================

As with the code, :ref:`commit and push <adding-code-to-mda>` your code to GitHub. Then :ref:`create a pull request <create-a-pull-request>`. The only test run for the User Guide is: that your file compile into HTML documents without errors. As usual, a developer will review your PR and merge the code into the User Guide when it looks good.

It is often difficult to review Jupyter notebooks on GitHub, especially if you embed widgets and images. One way to make it easier on the developers who review your changes is to build the changes on your forked repository and link the relevant sections in your pull request. To do this, create a ``gh-pages`` branch and merge your new branch into it. 

.. code-block:: bash

    # the first time
    git checkout -b gh-pages
    git merge origin/my-new-branch

Fix any merge conflicts that arise. Then edit ``UserGuide/doc/source/conf.py`` and change the URL of the site, which is set to ``site_url = "https://www.mdanalysis.org/UserGuide"``. Change it to your personal site, e.g. ::

    site_url = "https://www.my_user_name.github.io/UserGuide"


Now you can build your pages with the ``make github`` macro in the ``UserGuide/doc/`` directory, which builds the files and copies them to the top level of your directory.

.. code-block:: bash

    make github

You should be able to open one of these new HTML files (e.g. ``UserGuide/index.html``) in a browser and navigate your new documentation. Check that your changes look right. If they are, push to your `gh-pages` branch from the ``UserGuide/`` directory.

.. code-block:: bash

    git add .
    git commit -m 'built my-new-branch'
    git push -f origin gh-pages

On GitHub, navigate to your fork of the repository and go to **Settings**. In the **GitHub Pages** section, select the "gh-pages branch" from the **Source** dropdown. Check that your website is published at the given URL.

.. image:: images/gh-pages-settings.png

For each time you add changes to another branch later, just merge into gh-pages and rebuild. 

.. code-block:: bash

    git checkout gh-pages
    git merge origin/my_branch
    cd doc/
    make github

.. _autobuild-sphinx:

Automatically building documentation
====================================

Constantly rebuilding documentation can become tedious when you have many changes to make. Use `sphinx-autobuild <https://pypi.org/project/sphinx-autobuild>`_ to rebuild documentation every time you make changes to any document, including Jupyter notebooks. Install ``sphinx-autobuild``:

    .. code-block:: bash

        pip install sphinx-autobuild

Then, run the following command in the ``doc/`` directory:

    .. code-block:: bash

        python -m sphinx_autobuild source build

This will start a local webserver at http://localhost:8000/, which will refresh every time you save changes to a file in the documentation. This is helpful for both the user guide (first navigate to ``UserGuide/doc``) and the main repository documentation (navigate to ``package/doc/sphinx``).




