.. -*- coding: utf-8 -*-
.. _CONFIG-format:

========================================
CONFIG (DL_Poly Config)
========================================

.. include:: classes/CONFIG.txt

.. _HISTORY-format:

========================================
HISTORY (DL_Poly Config)
========================================

.. include:: classes/HISTORY.txt

MDAnalysis can read information both from DL Poly_ config and DL Poly history files. Although DL Poly input file units can be flexible, output files appear to have the following units:

    - Time: ps
    - Length: Angstrom
    - Mass: amu (Dalton)
    - Velocity: Angstrom/ps
    - Force: Angstrom Dalton / ps :sup:`2`

MDAnalysis currently does not convert these into the native kJ/(mol A) force units when reading files in. See `Issue 2393`_ for discussion on units.


.. _Poly: http://www.stfc.ac.uk/SCD/research/app/ccg/software/DL_POLY/44516.aspx

.. _`Issue 2393`: https://github.com/MDAnalysis/mdanalysis/issues/2393