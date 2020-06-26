.. -*- coding: utf-8 -*-
.. _GRO-format:

=======================================
GRO (GROMACS structure file)
=======================================

.. include:: classes/GRO.txt

GRO files provide topology, coordinate, and sometimes velocity information. 

Reading in
==========

Prior to MDAnalysis version 0.21.0 and GROMACS 2019.5, MDAnalysis failed to parse GRO files with box sizes where an axis length was longer than 10 characters. 

.. important::

    A Universe created with a GRO file and a Universe created with a corresponding TPR file will have *different* :ref:`atom and residue numbering <TPR-format>`, due to how a TPR file is parsed.

Writing out
===========

AtomGroups can be written out to a GRO file. However, this format does not support multi-frame trajectories.