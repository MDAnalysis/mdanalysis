===================================================
 HOLE analysis --- :mod:`MDAnalysis.analysis.hole2`
===================================================

:Author: Lily Wang
:Year: 2020
:Copyright: GNU Public License v3

.. versionadded:: 1.0.0

This module provides an updated interface for the HOLE_ suite of tools
:footcite:p:`Smart1993,Smart1996` to analyse an ion channel pore or transporter
pathway :footcite:p:`Stelzl2014` as a function of time or arbitrary order
parameters. It replaces :mod:`MDAnalysis.analysis.hole`.

HOLE_ must be installed separately and can be obtained in binary form
from http://www.holeprogram.org/ or as source from
https://github.com/osmart/hole2. (HOLE is open source and available
under the Apache v2.0 license.)


.. deprecated:: 2.8.0
  This module is deprecated in favour of the mdakit
  `mdahole2 <https://www.mdanalysis.org/mdahole2/>`_ and
  will be removed in MDAnalysis 3.0.0.


See Also
--------
:mod:`mdahole2.analysis.hole`


Module
------

.. automodule:: MDAnalysis.analysis.hole2
		


Utility functions and templates
-------------------------------

.. automodule:: MDAnalysis.analysis.hole2.utils
    :members:


.. automodule:: MDAnalysis.analysis.hole2.templates
    :members:
