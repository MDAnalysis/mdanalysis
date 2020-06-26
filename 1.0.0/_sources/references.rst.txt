.. -*- coding: utf-8 -*-
.. note: make sure that no lines accidentaly start with a single character
..       followed by a period: reST interprets it as an enumerated list and
..       messes up the formatting

.. The references are accessible globally; you can cite these papers anywhere
.. in the docs.

.. _references:

==========
References
==========

MDAnalysis and the included algorithms are scientific software that
are described in academic publications. **Please cite these papers when you use
MDAnalysis in published work.**

It is possible to :ref:`automatically generate a list of references
<citations-using-duecredit>` for any program that uses
MDAnalysis. This list (in common reference manager formats) contains
the citations associated with the specific algorithms and libraries
that were used in the program.




Citations using Duecredit
=========================

Citations can be automatically generated using duecredit_, depending on the
packages used. Duecredit is easy to install via ``pip``. Simply type:

.. code-block:: bash

   pip install duecredit

duecredit_ will remain an optional dependency, i.e. any code using
MDAnalysis will work correctly even without duecredit installed.

A list of citations for ``yourscript.py`` can be obtained using simple
commands.

.. code-block:: bash

   cd /path/to/yourmodule
   python -m duecredit yourscript.py

or set the environment variable :envvar:`DUECREDIT_ENABLE`

.. code-block:: bash

   DUECREDIT-ENABLE=yes python yourscript.py

Once the citations have been extracted (to a hidden file in the
current directory), you can use the :program:`duecredit` program to
export them to different formats. For example, one can display them in
BibTeX format, using:

.. code-block:: bash
 
   duecredit summary --format=bibtex 


**Please cite your use of MDAnalysis and the packages and algorithms
that it uses. Thanks!**


.. _duecredit: https://github.com/duecredit/duecredit

.. bibliography:: references.bib