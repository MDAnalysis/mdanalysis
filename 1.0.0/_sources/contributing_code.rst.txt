.. _working-with-mdanalysis-code:

=================================
Contributing to the main codebase
=================================

If you would like to contribute, start by searching through the `issues <https://github.com/MDAnalysis/mdanalysis/issues>`_ and :ref:`pull requests <https://github.com/MDAnalysis/mdanalysis/pulls>` to see whether someone else has raised a similar idea or question.

If you don't see your idea or problem listed, do one of the following:

    * If your contribution is **minor**, such as a typo fix, go ahead and fix it by following the guide below and :ref:`open a pull request <adding-code-to-mda>`.

    * If your contribution is **major**, such as a bug fix or a new feature, start by opening an issue first. That way, other people can weigh in on the discussion before you do any work. If you also create a pull request, you should link it to the issue by including the issue number in the pull request’s description.

Here is an overview of the development workflow for code or inline code documentation, as expanded on throughout the rest of the page.

    #. :ref:`Fork the MDAnalysis repository <forking-code-repo>` from the mdanalysis account into your own account
    #. :ref:`Set up an isolated virtual environment <create-virtual-environment>` for code development
    #. :ref:`Build development versions <build-mdanalysis-develop>` of MDAnalysis and MDAnalysisTests on your computer into the virtual environment
    #. :ref:`Create a new branch off the develop branch <create-code-branch>`
    #. :ref:`Add your new feature or bug fix <writing new code>` or :ref:`add your new documentation <guidelines-for-docstrings>`
    #. :ref:`Add and run tests <testing>` (if adding to the code)
    #. :ref:`Build and view the documentation <building-code-documentation>` (if adding to the docs)
    #. :ref:`Commit and push your changes, and open a pull request. <adding-code-to-mda>`


Working with the code
=====================

.. _forking-code-repo:

-------
Forking
-------

You will need your `own fork <https://help.github.com/en/github/getting-started-with-github/fork-a-repo>`_ to work on the code. Go to the `MDAnalysis project page <https://github.com/MDAnalysis/mdanalysis>`_ and hit the :guilabel:`Fork` button. You will
want to `clone your fork <https://help.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository>`_ to your machine:

.. code-block:: bash

    git clone https://github.com/your-user-name/mdanalysis.git
    cd mdanalysis
    git remote add upstream https://github.com/MDAnalysis/mdanalysis

This creates the directory `mdanalysis` and connects your repository to
the upstream (main project) MDAnalysis repository.

.. _create-virtual-environment:

----------------------------------
Creating a development environment
----------------------------------

To change code and test changes, you'll need to build both **MDAnalysis** and **MDAnalysisTests** 
from source. This requires a Python environment. We highly recommend that you use 
virtual environments. This allows you to have multiple experimental development versions 
of MDAnalysis that do not interfere with each other, or your own stable version. 
Since MDAnalysis is split into the actual package and a test suite, you need to install 
both modules in development mode.

You can do this either with :ref:`conda <dev-with-conda>` or :ref:`pip <dev-with-pip>`.

.. _dev-with-conda:

With conda
----------

Install either `Anaconda <https://www.anaconda.com/download/>`_ 
or `miniconda <https://conda.io/miniconda.html>`_.
Make sure your conda is up to date:

    .. code-block:: bash

        conda update conda

Create a new environment with ``conda create``. This will allow you to change code in 
an isolated environment without touching your base Python installation, and without 
touching existing environments that may have stable versions of MDAnalysis. :

    .. code-block:: bash

        conda create --name mdanalysis-dev

Activate the environment to build MDAnalysis into it:

    .. code-block:: bash

        conda activate mdanalysis-dev

To view your environments:

    .. code-block:: bash

        conda info -e

To list the packages installed in your current environment:

    .. code-block:: bash

        conda list

To return to your root environment:

    .. code-block:: bash

        conda deactivate

See the full conda docs `here <http://conda.pydata.org/docs>`__.

.. _dev-with-pip:

With pip and virtualenv
-----------------------

Like conda, virtual environments managed with `virtualenv <https://virtualenv.pypa.io/en/latest/>`_ allow you to use different versions of Python and Python packages for your different project. Unlike conda, virtualenv is not a general-purpose package manager. Instead, it leverages what is available on your system, and lets you install Python packages using pip.

To use virtual environments you have to install the virtualenv package first. This can be done with either pip or the package manager of your system:

    .. code-block:: bash

        pip install virtualenv
        # or on ubuntu
        sudo apt install virtualenv
        # or on fedora
        sudo dnf install python-virtualenv

Virtual environments can be created for each project directory.

    .. code-block:: bash

        cd my-project/
        virtualenv my-project-env

This will create a new folder ``my-project-env``. This folder contains the virtual environment and all packages you have installed in it. To activate it in the current terminal run:

    .. code-block:: bash

        source myproject-env/bin/activate

Now you can install packages via pip without affecting your global environment. The packages that you install when the environment is activated will be available in terminal sessions that have the environment activated. You can deactivate the virtual environment by running:

    .. code-block:: bash

        deactivate

The `virtualenvwrapper package <https://virtualenvwrapper.readthedocs.io/en/latest/>`_ makes virtual environments easier to use. It provides some very useful features:

    - it organises the virtual environment into a single user-defined directory, so they are not scattered throughout the file system;
    - it defines commands for the easy creation, deletion, and copying of virtual environments;
    - it defines a command to activate a virtual environment using its name;
    - all commands defined by ``virtualenvwrapper`` have tab-completion for virtual environment names.

You first need to install ``virtualenvwrapper`` *outside* of a virtual environment:

    .. code-block:: bash

        pip install virtualenvwrapper
        # or on ubuntu
        sudo apt install virtualenvwrapper
        # or on fedora
        sudo dnf install python-virtualenvwrapper

Then, you need to load it into your terminal session. Add the following lines in ``~/.bashrc``. They will be executed every time you open a new terminal session:

    .. code-block:: bash

        # Decide where to store the virtual environments
        export WORKON_HOME=~/Envs
        # Make sure the directory exists
        mkdir -p ${WORKON_HOME}
        # Load virtualenvwrapper
        source /usr/local/bin/virtualenvwrapper.sh

Open a new terminal or run ``source ~/.bashrc`` to update your session. You can now create a virtual environment with:

    .. code-block:: bash

        mkvirtualenv my-project

Regardless of your current working directory, the environment is created in ``~/Envs/`` and it is now loaded in our terminal session.

You can load your virtual environments by running ``workon my-project``, and exit them by running ``deactivate``.

Virtual environments, especially with ``virtualenvwrapper``, can do much more. For example, you can create virtual environments with different python interpreters with the ``-p`` flag. The Hitchhiker's Guide to Python has a good `tutorial <https://docs.python-guide.org/dev/virtualenvs/>`_ that gives a more in-depth explanation of virtual environments. The `virtualenvwrapper documentation <https://virtualenvwrapper.readthedocs.io/en/latest/>`_ is also a good resource to read.

On a Mac
--------

One more step is often required on macOS, because of the default number of files that a process can open simultaneously is quite low (256). To increase the number of files that can be accessed, run the following command:

    .. code-block:: bash

        ulimit -n 4096

This sets the number of files to 4096. However, this command only applies to your currently open terminal session. To keep this high limit, add the above line to your ``~/.profile``.



.. _build-mdanalysis-develop:

-------------------
Building MDAnalysis
-------------------

Make sure that you have :ref:`cloned the repository <forking-code-repo>`  
and activated your virtual environment. First we need to install dependencies:

    .. code-block:: bash

        # if using conda
        conda install -c biobuilds -c conda-forge \
            pip cython numpy mmtf-python mock six biopython \
            networkx cython matplotlib scipy griddataformats \
            hypothesis gsd codecov "seaborn>=0.7.0,<=0.9" \
            clustalw=2.1 netcdf4 scikit-learn "joblib>=0.12"\
            psutil pytest

        # if using conda with python 3.7 or 3.8, also run
        conda install -c conda-forge parmed

        # if using conda with other versions of python, also run
        pip install parmed

    .. code-block:: bash

        # if using pip and virtualenv
        pip install cython numpy mmtf-python mock six biopython \
            networkx cython matplotlib scipy griddataformats \
            hypothesis gsd codecov "seaborn>=0.7.0,<=0.9" \
            netcdf4 scikit-learn "joblib>=0.12" parmed psutil pytest

Ensure that you have a working C/C++ compiler (e.g. gcc or clang). You will also need Python ≥ 3.4. We will now install MDAnalysis. 

    .. code-block:: bash

        # go to the mdanalysis source directory
        cd mdanalysis/

        # Build and install the MDAnalysis package
        cd package/
        pip install -e .
        
        # Build and install the test suite
        cd ../testsuite/
        pip install -e .

At this point you should be able to import MDAnalysis from your locally built version. If you are running the development version, this is visible from the version number ending in "-dev0". For example:

    .. code-block:: bash

        $ python  # start an interpreter
        >>> import MDAnalysis as mda
        >>> mda.__version__
        '0.20.2-dev0'

If your version number does not end in "-dev0", you may be on the ``master`` branch. In your ``mdanalysis/`` directory, switch to the ``develop`` branch:

    .. code-block:: bash

        $ git checkout develop
        Switched to branch 'develop'


.. _branches-in-mdanalysis:

----------------------
Branches in MDAnalysis
----------------------

There are two important branches in MDAnalysis:

    - ``master``: for production-ready code
    - ``develop``: for development code

The ``master`` branch is only for stable, production-ready code. Development code should *never* be committed to this branch. Typically, code is only committed by the release manager, when a release is ready.

The ``develop`` branch can be considered an "integration" branch for including your code into the next release. Only working, tested code should be committed to this branch. Code contributions ("features") should branch off ``develop`` rather than ``master``.


.. _create-code-branch:

Creating a branch
-----------------

The develop branch should only contain approved, tested code, so create a
feature branch for making your changes. For example, to create a branch called 
``shiny-new-feature`` from ``develop``:

    .. code-block:: bash

        git checkout -b shiny-new-feature develop

This changes your working directory to the ``shiny-new-feature`` branch.  Keep any
changes in this branch specific to one bug or feature so it is clear
what the branch brings to MDAnalysis. You can have many branches with different names
and switch in between them using the ``git checkout my-branch-name`` command.

There are several special branch names that you should not use for your feature branches:

    - ``master``
    - ``develop``
    - ``release-*``


``release`` branches are used to :ref:`prepare a new production release <preparing-release>` and should be handled by the release manager only.

.. _writing-new-code:

----------------
Writing new code
----------------

Code formatting in Python
-------------------------

MDAnalysis is a project with a long history and many contributors; it hasn't used a consistent coding style. Since version 0.11.0, we are trying to update all the code to conform with `PEP8`_. Our strategy is to update the style every time we touch an old function and thus switch to `PEP8`_ continuously.

**Important requirements (from PEP8):**
    - keep line length to **79 characters or less**; break long lines sensibly
    - indent with **spaces** and use **4 spaces per level**
    - naming:

        - classes: `CapitalClasses` (i.e. capitalized nouns without spaces)
        - methods and functions: `underscore_methods` (lower case, with underscores for spaces)

We recommend that you use a Python Integrated Development Environment (IDE) (`PyCharm`_ and others) or external tools like `flake8`_ for code linting. For integration of external tools with emacs and vim, check out `elpy`_ (emacs) and `python-mode`_ (vim).

To apply the code formatting in an automated way, you can also use code formatters. External tools include `autopep8`_ and `yapf`_. Most IDEs either have their own code formatter or will work with one of the above through plugins.


.. _`PEP8`: https://www.python.org/dev/peps/pep-0008/
.. _`flake8`: http://flake8.readthedocs.org/en/latest/
.. _`PyCharm`: https://www.jetbrains.com/pycharm/
.. _`elpy`: https://github.com/jorgenschaefer/elpy
.. _`python-mode`: https://github.com/klen/python-mode
.. _`autopep8`: https://github.com/hhatto/autopep8
.. _`yapf`: https://github.com/google/yapf


Modules and dependencies
------------------------

MDAnalysis strives to keep dependencies small and lightweight. Code outside the :mod:`MDAnalysis.analysis` and :mod:`MDAnalysis.visualization` modules should only rely on the :ref:`core dependencies <core-module-dependencies>`, which are always installed. Analysis and visualization modules can use any :ref:`any package, but the package is treated as optional <optional-modules>`.

Imports in the code should follow the :ref:`general-rules-for-importing`.

.. seealso::

    See :ref:`module-imports` for more information.


Developing in Cython
--------------------

The ``setup.py`` script first looks for the `.c` files included in the standard MDAnalysis distribution. These are not in the GitHub repository, so ``setup.py`` will use Cython to compile extensions. `.pyx` source files are used instead of `.c` files. From there, `.pyx` files are converted to `.c` files if they are newer than the already present `.c` files or if the ``--force`` flag is set (i.e. ``python setup.py build --force``). End users (or developers) should not trigger the `.pyx` to `.c` conversion, since `.c` files delivered with source packages are always up-to-date. However, developers who work on the `.pyx` files will automatically trigger the conversion since `.c` files will then be outdated. 

Place all source files for compiled shared object files into the same directory as the final shared object file.

`.pyx` files and cython-generated `.c` files should be in the same directory as the `.so` files. External dependent C/C++/Fortran libraries should be in dedicated ``src/`` and ``include/`` folders. See the following tree as an example:

    ::

        MDAnalysis 
            |--lib
            |   |-- _distances.so
            |   |-- distances.pyx
            |   |-- distances.c
            |-- coordinates
                |-- _dcdmodule.so
                |-- src
                    |-- dcd.c
                |-- include
                    |-- dcd.h

.. _test-code:

-----------------
Testing your code
-----------------

MDAnalysis takes testing seriously. All code added to MDAnalysis should have tests to ensure that it works as expected; we aim for 90% coverage. See :ref:`testing` for more on :ref:`writing <write-new-tests>`, :ref:`running <run-test-suite>`, and interpreting tests.


---------------------
Documenting your code
---------------------

Changes to the code should be reflected in the ongoing ``CHANGELOG``. Add an entry here to document your fix, enhancement, or change. In addition, add your name to the author list. If you are addressing an issue, make sure to include the issue number.


.. _adding-code-to-mda:

------------------------------
Adding your code to MDAnalysis
------------------------------

Committing your code
--------------------

When you are happy with a set of changes and :ref:`all the tests pass <test-code>`, it is time to commit. All changes in one revision should have a common theme. If you implemented two rather different things (say, one bug fix and one new feature), then split them into two commits with different messages.

Once you’ve made changes to files in your local repository, you can see them by typing:

    .. code-block:: bash

        git status

Tell git to track files by typing:

    .. code-block::

        git add path/to/file-to-be-added.py

Doing ``git status`` again should give something like:

    .. code-block::

        # On branch shiny-new-feature
        #
        #       modified:   /relative/path/to/file-you-added.py
        #

Then commit with:

    .. code-block:: bash

        git commit -m

This opens up a message editor. 

*Always* add a descriptive comment for your commit message (feel free to be verbose!):

    - use a short (<50 characters) subject line that summarizes the change
    - leave a blank line
    - optionally, add additional more verbose descriptions; paragraphs or bullet lists (with ``-`` or ``*``) are good
    - manually break lines at 80 characters
    - manually indent bullet lists

.. seealso::

    See `Tim Pope's A Note About Git Commit Messages <http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html>`_ for a rationale for these rules.


Pushing your code to GitHub
---------------------------

When you want your changes to appear publicly on your GitHub page, push your forked feature branch’s commits:

    .. code-block:: bash

        git push origin shiny-new-feature

Here `origin` is the default name given to your remote repository on GitHub. You can see the remote repositories:

    .. code-block:: bash

        git remote -v

If you added the upstream repository as described above you will see something like:

    .. code-block:: bash

        origin	git@github.com:your-username/mdanalysis.git (fetch)
        origin	git@github.com:your-username/mdanalysis.git (push)
        upstream	git@github.com:MDAnalysis/mdanalysis.git (fetch)
        upstream	git@github.com:MDAnalysis/mdanalysis.git (push)

Now your code is on GitHub, but it is not yet a part of the MDAnalysis project. For that to happen, a pull request needs to be submitted on GitHub. 

.. _rebase-code:

Rebasing your code
------------------

Often the upstream MDAnalysis develop branch will be updated while you are working on your own code.
You will then need to update your own branch with the new code to avoid merge conflicts.
You need to first retrieve it and then `rebase <https://www.atlassian.com/git/tutorials/rewriting-history/git-rebase>`_
your branch so that your changes apply to the new code:

    .. code-block:: bash

        git fetch upstream
        git rebase upstream/develop

This will replay your commits on top of the latest development code from MDAnalysis.  If this
leads to merge conflicts, you must resolve these before submitting your pull
request.  If you have uncommitted changes, you will need to ``git stash`` them
prior to updating.  This will effectively store your changes and they can be
reapplied after updating with ``git stash apply``. 

Once rebased, push your changes:

    .. code-block:: bash

        git push -f origin shiny-new-feature

and `create a pull request <https://github.com/MDAnalysis/mdanalysis/pulls>`_.

.. _create-a-pull-request:

Creating a pull request
-----------------------

The typical approach to adding your code to MDAnalysis is to make a `pull request <https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests>`_ on GitHub. Please make sure that your contribution :ref:`passes all tests <test-code>`. If there are test failures, you will need to address them before we can review your contribution and eventually merge them. If you have problems with making the tests pass, please ask for help! (You can do this in the comments of the pull request). 

    #. Navigate to your repository on GitHub
    #. Click on the :guilabel:`Pull Request` button
    #. You can then click on :guilabel:`Commits` and :guilabel:`Files Changed` to make sure everything looks okay one last time
    #. Write a description of your changes and follow the PR checklist

        - check that docs are updated
        - check that tests run
        - check that you've updated CHANGELOG
        - reference the issue that you address, if any

    #. Click :guilabel:`Send Pull Request`.

Your pull request is then sent to the repository maintainers. After this, the following happens:

    #. A :ref:`suite of tests are run on your code <continuous-integration>` with the tools :ref:`travis`, :ref:`appveyor` and :ref:`codecov`. If they fail, please fix your pull request by pushing updates to it.
    #. Developers will ask questions and comment in the pull request. You may be asked to make changes. 
    #. When everything looks good, a core developer will merge your code into the ``develop`` branch of MDAnalysis. Your code will be in the next release.

If you need to make changes to your code, you can do so on your local repository as you did before. Committing and pushing the changes will  update your pull request and restart the automated tests.

.. _working-with-mdanalysis-docs:

Working with the code documentation
===================================

MDAnalysis maintains two kinds of documentation: 

    #. `This user guide <https://www.mdanalysis.org/UserGuide/>`__: a map of how MDAnalysis works, combined with tutorial-like overviews of specific topics (such as the analyses)
    
    #. `The documentation generated from the code itself <https://www.mdanalysis.org/docs/>`__. Largely built from code docstrings, these are meant to provide a clear explanation of the usage of individual classes and functions. They often include technical or historical information such as in which version the function was added, or deprecation notices.

This guide is for the documentation generated from the code. If you are looking to contribute to the user guide, please see :ref:`working-with-user-guide`.

MDAnalysis has a lot of documentation in the Python doc strings. The docstrings follow the `Numpy Docstring Standard <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`__, which is used widely
in the Scientific Python community. They are nice to read as normal text and are converted by sphinx to normal ReST through `napoleon <http://sphinxcontrib-napoleon.readthedocs.org/en/latest/index.html>`__.

This standard specifies the format of
the different sections of the docstring. See `this document
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
for a detailed explanation, or look at some of the existing functions to
extend it in a similar manner.

Note that each page of the  `online documentation <https://www.mdanalysis.org/docs/>`_ has a link to the *Source* of the page. You can look at it in order to find out how a particular page has been written in reST and copy the approach for your own documentation.

.. _building-code-documentation:

--------------------------
Building the documentation
--------------------------

The online documentation is generated from the pages in ``mdanalysis/package/doc/sphinx/source/documentation_pages``. The documentation for the current release are hosted at www.mdanalysis.org/docs, while the development version is at www.mdanalysis.org/mdanalysis/. 

In order to build the documentation, you must first :ref:`clone the main MDAnalysis repo <forking-code-repo>`. :ref:`Set up a virtual environment <create-virtual-environment>` in the same way as you would for the code (you can use the same environment as you do for the code). You will need to install several packages for the docs.

    .. code-block:: bash

        pip install sphinx sphinx-sitemap sphinx_rtd_theme

In addition, build the development version of MDAnalysis (if you haven't done this already):

    .. code-block:: bash

        pip install -e .

Then, generate the docs with:

    .. code-block:: bash

        python setup.py build_sphinx -E

This generates and updates the files in ``doc/html``. If the above command fails with an ``ImportError``, run

    .. code-block:: bash

        python setup.py build_ext --inplace

and retry.

You will then be able to open the home page, ``doc/html/index.html``, and look through the docs. In particular, have a look at any pages that you tinkered with. It is typical to go through multiple cycles of fix, rebuild the docs, check and fix again.

If rebuilding the documentation becomes tedious after a while, install the :ref:`sphinx-autobuild <autobuild-sphinx>` extension. 

-------------------------
Where to write docstrings
-------------------------

When writing Python code, you should always add a docstring to each public (visible to users):

    * module
    * function
    * class
    * method
 
\When you add a new module, you should include a docstring with a short sentence describing what the module does, and/or a long document including examples and references. 

.. _guidelines-for-docstrings:

---------------------------------
Guidelines for writing docstrings
---------------------------------

A typical function docstring looks like the following:

    ::

        def func(arg1, arg2):
            """Summary line.

            Extended description of function.

            Parameters
            ----------
            arg1 : int
                Description of `arg1`
            arg2 : str
                Description of `arg2`


            Returns
            -------
            bool
                Description of return value

            """
            return True

.. seealso::

    The `napoleon documentation <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`_ has further breakdowns of docstrings at the module, function, class, method, variable, and other levels.

* When writing reST markup, make sure that there are **at least two blank lines above** the reST after a numpy heading. Otherwise, the Sphinx/napoleon parser does not render correctly.

    .. code-block:: RST

        some more docs bla bla

        Notes
        -----
        THE NEXT TWO BLANK LINES ARE IMPORTANT.


        .. versionadded:: 0.16.0
  
* Do not use "Example" or "Examples" as a normal section heading (e.g. in module level docs): *only* use it as a `NumPy doc Section <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`__. It will not be rendered properly, and will mess up sectioning.


* When writing multiple common names in one line, Sphinx sometimes tries to reference the first name. In that case, you have to split the names across multiple lines. See below for an example:

    .. code-block:: RST

        Parameters
        ----------
        n_atoms, n_residues : int
            numbers of atoms/residues

* We are using MathJax with sphinx so you can write LaTeX code in math tags. 

    In blocks, the code below

        .. code-block:: rst

            #<SPACE if there is text above equation>
            .. math::
                e^{i\pi} = -1

    renders like so:

        .. math::
            e^{i\pi} = -1
    

    Math directives can also be used inline.

        .. code-block:: rst

            We make use of the identity :math:`e^{i\pi} = -1` to show...

    Note that you should *always* make doc strings with math code **raw** python strings **by prefixing them with the letter "r"**, or else you will get problems with backslashes in unexpected places.

        ::

            def rotate(self, R):
                r"""Apply a rotation matrix *R* to the selection's coordinates.

                :math:`\mathsf{R}` is a 3x3 orthogonal matrix that transforms a vector
                :math:`\mathbf{x} \rightarrow \mathbf{x}'`:

                .. math::

                \mathbf{x}' = \mathsf{R}\mathbf{x}
                """

    .. seealso::
    
        See `Stackoverflow: Mathjax expression in sphinx python not rendering correctly <http://stackoverflow.com/questions/16468397/mathjax-expression-in-sphinx-python-not-rendering-correclty">`_ for further discussion.


-------------------
Documenting changes
-------------------

.. _versionadded: https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#directive-versionadded
.. _versionchanged: https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#directive-versionchanged
.. _deprecated: https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#directive-deprecated

We use reST constructs to annotate *additions*, *changes*, and *deprecations* to the code so that users can quickly learn from the documentation in which version of MDAnalysis the feature is available.

A **newly added module/class/method/attribute/function** gets a `versionadded`_  directive entry in its primary doc section, as below.

.. code-block:: rst

   .. versionadded:: X.Y.Z

For parameters and attributes, we typically mention the new entity in a `versionchanged`_ section of the function or class (although a `versionadded`_ would also be acceptable).

**Changes** are indicated with a `versionchanged`_ directive

.. code-block:: rst

   .. versionchanged:: X.Y.Z
      Description of the change. Can contain multiple descriptions.
      Don't assume that you get nice line breaks or formatting, write your text in
      full sentences that can be read as a paragraph.

**Deprecations** (features that are not any longer recommended for use and that will be removed in future releases) are indicated by the `deprecated`_ directive:

.. code-block:: rst

   .. deprecated:: X.Y.Z
      Describe (1) alternatives (what should users rather use) and 
      (2) in which future release the feature will be removed.

When a feature is removed, we remove the deprecation notice and add a `versionchanged`_ to the docs of the enclosing scope. For example, when a parameter of a function is removed, we update the docs of the function. Function/class removal are indicated in the module docs. When we remove a whole module, we typically indicate it in the top-level reST docs that contain the TOC tree that originally included the module.



--------------------------------------
Writing docs for abstract base classes
--------------------------------------

MDAnalysis contains a number of abstract base classes, such as :class:`~MDAnalysis.analysis.base.AnalysisBase`. Developers who define new base classes, or modify existing ones, should follow these rules:

    - The *class docstring* needs to contain a list of methods that can be overwritten by inheritance from the base class. Distinguish and document methods as required or optional.
    - The class docstring should contain a minimal example for how to derive this class. This demonstrates best practices, documents ideas and intentions behind the specific choices in the API, helps to promote a unified code base, and is useful for developers as a concise summary of the API.
    - A more detailed description of methods should come in the *method docstring*, with a note specifying if the method is required or optional to overwrite.

See the documentation of :class:`MDAnalysis.analysis.base.AnalysisBase` for an example of this documentation.

---------------------------------------
Adding your documentation to MDAnalysis
---------------------------------------

As with any contribution to an MDAnalysis repository, :ref:`commit and push <adding-code-to-mda>` your documentation contributions to GitHub. If *any fixes in the restructured text* are needed, *put them in their own commit* (and do not include any generated files under `docs/html`). Try to keep all reST fixes in the one commit. ``git add FILE`` and ``git commit --amend`` is your friend when piling more and more small reST fixes onto a single "fixed reST" commit.

We recommend :ref:`building the docs locally first <building-code-documentation>` to preview your changes. Then, :ref:`create a pull request <create-a-pull-request>`. All the tests in the MDAnalysis test suite will run, but only one checks that the documents compile correctly.

---------------------------------------
Viewing the documentation interactively
---------------------------------------

In the Python interpreter one can simply say:

    ::

        import MDAnalysis
        help(MDAnalysis)
        help(MDAnalysis.Universe)

In ``ipython`` one can use the question mark operator:

    .. ipython::
        :verbatim:

        In [1]: MDAnalysis.Universe?
