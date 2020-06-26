.. -*- coding: utf-8 -*-
.. _preparing-release:

===================
Preparing a release
===================

Rules for a release branch:

    - May branch from: ``develop``
    - Must be merged into: ``master`` (and ``develop`` if needed)
    - Naming convention: ``release-*`` where ``*`` is a version number

Release policy and release numbering
====================================

We use a **MAJOR.MINOR.PATCH** scheme to label releases. We adhere to the idea of `semantic versioning <http://semver.org/>`_ (semantic versioning was introduced with release 0.9, see `Issue 200`_): Given a version number **MAJOR.MINOR.PATCH**, we increment the:

  * **MAJOR** version when we make **incompatible API changes**,
  * **MINOR** version when we **add functionality** in a backwards-compatible manner, and
  * **PATCH** version when we make backwards-compatible **bug fixes**.

However, as long as the **MAJOR** number is **0** (i.e. the API has not stabilized), even **MINOR** increases *may* introduce incompatible API changes. As soon as we have a 1.0.0 release, the public API can only be changed in a backward-incompatible manner with an increase in MAJOR version.

Additional labels for pre-release and build metadata are available as extensions to the MAJOR.MINOR.PATCH format.

The `CHANGELOG <https://github.com/MDAnalysis/mdanalysis/blob/develop/package/CHANGELOG>`_ lists important changes for each release.

*MAJOR*, *MINOR* and *PATCH* numbers are integers that increase monotonically.

The **release number** is set in `setup.py <https://github.com/MDAnalysis/mdanalysis/blob/develop/package/setup.py>`_ *and* in ``MDAnalysis.__version__`` (MDAnalysis/version.py), e.g. ::

    RELEASE = '0.7.5'


**While the code is in development** (i.e. whenever we are not preparing a release!) the release number gets the suffix **-dev0**, e.g. ::

    RELEASE = '0.7.6-dev0'

so that people using the :ref:`develop branch <branches-in-mdanalysis>` from the source repository can immediately see that it is not a final release. For example, "0.7.6-dev0" is the state *before* the 0.7.6 release.

.. _`Issue 200`: https://github.com/MDAnalysis/mdanalysis/issues/200

Typical workflow for preparing a release
========================================

#. Declare feature freeze on ``develop`` via the `developer mailing list`_

#. Create a release branch from ``develop``::

    git checkout -b release-0.7.6 develop

#. Finalise the ``CHANGELOG`` with the release number and date. Summarize important changes and add all authors that contributed to this release.

#. Make sure the version number is right::

    ./maintainer/change_release.sh 0.7.6

#. Check that the documentation is up-to-date and tests pass. Check that any new Cython code has compiled.

#. Commit the finalized state::

    git commit -m "release 0.7.6 ready"

#. Build a source distribution tarballs under ``package/dist/MDAnalysis-MAJOR-MINOR-PATCH.tar.gz`` and ``testsuite/dist/MDAnalysisTests-MAJOR-MINOR-PATCH.tar.gz``:

    .. code-block:: bash

        # MDAnalysis
        cd package/
        python setup.py sdist

        # MDAnalysisTests
        cd ../testsuite/
        python setup.py sdist

#. Test the distribution in a ``tmp`` directory.

    #. Unpack and try to build it:

        .. code-block:: bash

            mkdir tmp && cd tmp
            tar -zxvf ../dist/MDAnalysis-0.7.5.tar.gz
            cd MDAnalysis-0.7.5
            python setup.py build --build-lib=.

    
    #. Run the tests again:

        .. code-block::

            python
            >>> import MDAnalysis.tests
            >>> MDAnalysis.tests.test(label='full', extra_argv=['--exe'])

        
        The above should work at least on Linux and Mac OS X. If it fails then go back and fix things and *do not release*.


#. If everything works, merge the branch into master and tag the release::

    git checkout master
    git merge --no-ff release-0.7.6
    git tag -m 'release 0.7.5 of MDAnalysis and MDAnalysisTests' release-0.7.5
    git push --tags origin master

#. Merge the branch back into ``develop`` (this is not required if the only change was the version number)::

    git checkout develop
    git merge --no-ff release-0.7.6
    ./maintainer/change_release.sh 0.7.7-devel
    git commit -a -m "version number changed to 0.7.7-devel"

#. Build and deploy the docs manually. (You may need to first ``pip install sphinx==2.2.0 sphinx_sitemap sphinx_rtd_theme``)::

    cd package/
    python setup.py build_sphinx
    cd ..

    # You need a OAUTH token that gives commit access to the MDAnalysis/docs repo
    export GH_TOKEN=<secret>

    ./maintainer/deploy_master_docs.sh

#. Update the release on the Python package index (Pypi)

    #. Upload the package to Pypi. You need to have run ``python setup.py register`` previously.

        .. code-block:: bash

            twine upload -r pypi dist/MDAnalysis-0.16.2.tar.gz 

    #. Upload the docs to Pypi

    #. Make the new tar ball a *featured* release so that it shows up on the front page (and *unfeature* any older releases).

    #. Provide a short description (a condensed version of the ``CHANGELOG``)

#. Update the release on Anaconda

    conda packages are built on conda-forge.

    #. Create a pull request from https://github.com/MDAnalysis/mdanalysis-feedstock for https://github.com/conda-forge/mdanalysis-feedstock
    #. Create a pull request from https://github.com/MDAnalysis/mdanalysistests-feedstock to https://github.com/conda-forge/mdanalysistests-feedstock

#. Create a ReleaseXYZ wiki page, modelled after e.g. `Release062 <https://github.com/MDAnalysis/mdanalysis/wiki/Release062>`_ (using the ``CHANGELOG`` as a reference). Add it to the `Release Notes <https://github.com/MDAnalysis/mdanalysis/wiki/Release-Notes>`_.


#. Delete the release branch::

    git branch -d release-0.7.6

.. _`developer mailing list`: https://groups.google.com/forum/#!forum/mdnalysis-devel
