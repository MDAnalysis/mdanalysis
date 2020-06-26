.. -*- coding: utf-8 -*-
.. _ITP-format:

=====================================
ITP (GROMACS portable topology files)
=====================================

.. include:: classes/ITP.txt

A ITP_ file is a portable topology file. 

.. important::

    Unlike :ref:`TPR <TPR-format>` files, atom ``ids`` and residues ``resids`` in ITP files are indexed from 1. This means that a TPR file created from your ITP files will have *different* numbering in MDAnalysis than the ITP file.

.. _ITP: http://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html#molecule-itp-file