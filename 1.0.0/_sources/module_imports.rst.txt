.. -*- coding: utf-8 -*-
.. _module-imports:

============================
Module imports in MDAnalysis
============================

We are striving to keep module dependencies small and lightweight (i.e., easily installable with ``pip``).

.. _general-rules-for-importing:

General rules for importing
===========================

    - Imports should all happen at the start of a module (not inside classes or functions).  
    - Modules must be imported in the following order:

        - `future <https://docs.python.org/2/library/__future__.html>`_ (``from __future__ import absolute_import, print_function, division``)
        - Compatibility imports (e.g. ``six``)
        - global imports (installed packages)
        - local imports (MDAnalysis modules)
    - use **absolute imports** in the library (i.e. relative imports must be explicitly indicated)

For example::

    from __future__ import absolute_import
    from six.moves import range

    import numpy as np

    import .core
    import ..units

--------------------------------------------
Module imports in :mod:`MDAnalysis.analysis`
--------------------------------------------

#. In :mod:`MDAnalysis.analysis`, all imports must be at the top level (as in the General Rule) — see `Issue 666`_ for more information.
#. :ref:`Optional modules <optional-modules>` can be imported
#. No analysis module is imported automatically at the :mod:`MDAnalysis.analysis` level to avoid breaking the installation when optional dependencies are missing.

.. _`Issue 666`: https://github.com/MDAnalysis/mdanalysis/issues/666

.. _module-imports-in-tests:

--------------------------------
Module imports in the test suite
--------------------------------

    - Use the module import order in :ref:`general-rules-for-importing`, but import :mod:`MDAnalysis` modules before :mod:`MDAnalysisTests` imports
    - Do not use *relative imports* (e.g. ``import .datafiles``) in the test suite. This breaks running the tests from inside the test directory (see `Issue 189`_ for more information)
    - Never import the :mod:`MDAnalysis` module from the ``__init__.py`` of :mod:`MDAnalysisTests` or from any of its plugins (it's ok to import from test files). Importing :mod:`MDAnalysis` from the test setup code will cause severe coverage underestimation.

.. _`Issue 189`: https://github.com/MDAnalysis/mdanalysis/issues/189

Module dependencies in the code
===============================

.. _core-module-dependencies:

--------------------------------
List of core module dependencies
--------------------------------

Any module from the standard library can be used, as well as the following nonstandard libraries:

   * :mod:`six`
   * :mod:`numpy`
   * :mod:`biopython`
   * :mod:`gridDataFormats`
   * :mod:`mmtf-python`
   * :mod:`scipy`
   * :mod:`matplotlib`

because these packages are always installed.

If you must depend on a new external package, first discuss its use on the `developer mailing list`_ or as part of the issue/pull request.

.. _`developer mailing list`: https://groups.google.com/forum/#!forum/mdnalysis-devel


.. _core-modules:

---------------------
Modules in the "core"
---------------------

The core of MDAnalysis contains all packages that are not in :mod:`MDAnalysis.analysis` or :mod:`MDAnalysis.visualization`. Only packages in the :ref:`core-module-dependencies` can be imported.

.. _optional-modules:

----------------------------------------------------------------------------------
Optional modules in :mod:`MDAnalysis.analysis` and :mod:`MDAnalysis.visualization`
----------------------------------------------------------------------------------

Modules under :mod:`MDAnalysis.analysis` are considered independent from the core package. Each analysis module can have its own set of dependencies. We strive to keep them small, but module authors are, in principle, free to import what they need. When analysis modules import packages outside of :ref:`core-module-dependencies`, the dependencies are considered **optional** (and should be listed in ``setup.py`` under *analysis*). (See also `Issue 1159`_ for more discussion.)

A user who does not have a specific optional package installed must still be able to import everything else in MDAnalysis. An analysis module *may* simply raise an :code:`ImportError` if a package is missing. However, it is recommended that the module should print and log an *error message* notifying the user that a specific additional package needs to be installed to use it.

If a large portion of the code in the module does not depend on a specific optional module then you should:

    - guard the import at the top level with a :code:`try/except`
    - print and log a *warning*
    - only raise an :code:`ImportError` in the specific function or method that would depend on the missing module.

.. _`Issue 1159`: https://github.com/MDAnalysis/mdanalysis/issues/1159
