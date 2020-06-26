.. -*- coding: utf-8 -*-
.. _LAMMPSDUMP-format:

=========================================
LAMMPSDUMP (LAMMPS ascii dump file)
=========================================

.. include:: classes/LAMMPSDUMP.txt

Reading in
==========

MDAnalysis expects ascii dump files to be written with the default
`LAMMPS dump format`_ of 'atom'.It will automatically convert positions from their scaled/fractional representation to their real values.

.. important::

    Lennard-Jones units are not implemented. See :ref:`units`
    for other recognized values and the documentation for the LAMMPS
    `units command`_.

.. _units command: http://lammps.sandia.gov/doc/units.html
.. _`LAMMPS dump format`: http://lammps.sandia.gov/doc/dump.html