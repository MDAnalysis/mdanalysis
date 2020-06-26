.. -*- coding: utf-8 -*-
.. _DCD-format:

========================================
DCD (CHARMM, NAMD, or LAMMPS trajectory)
========================================

.. include:: classes/DCD.txt

DCD is used by NAMD, CHARMM and LAMMPS as the default trajectory format. 



Reading in
==========

**Unitcell dimensions**

Generally, DCD trajectories produced by any code can be read (with the :class:`~MDAnalysis.coordinates.DCD.DCDReader`) although there can be issues with the unitcell dimensions (simulation box). 
Currently, MDAnalysis tries to guess the correct **format for the unitcell representation** but it can be wrong. 
**Check the unitcell dimensions**, especially for triclinic unitcells (see `Issue 187`_). 

MDAnalysis always uses ``(*A*, *B*, *C*, *alpha*, *beta*, *gamma*)`` to
represent the unit cell. Lengths *A*, *B*, *C* are in the MDAnalysis length
unit (Å), and angles are in degrees.

The ordering of the angles in the unitcell is the same as in recent
versions of VMD's DCDplugin_ (2013), namely the `X-PLOR DCD format`_: The
original unitcell is read as ``[A, gamma, B, beta, alpha, C]`` from the DCD
file. If any of these values are < 0 or if any of the angles are > 180
degrees then it is assumed it is a new-style CHARMM unitcell (at least
since c36b2) in which box vectors were recorded.

.. important:: 

    Check your unit cell dimensions carefully, especially when using triclinic boxes. Old CHARMM trajectories might give wrong unitcell values.

**Units**

The DCD file format is not well defined. In particular, NAMD and CHARMM use it differently. 
DCD trajectories produced by CHARMM and NAMD( >2.5) record time in AKMA units. 
If other units have been recorded (e.g., ps) then employ the configurable :ref:`LAMMPS DCD format <LAMMPS-format>` and set the time unit as an optional argument. 
You can find a list of units used in the DCD formats on the MDAnalysis `wiki`_.

.. _`X-PLOR DCD format`: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
.. _Issue 187: https://github.com/MDAnalysis/mdanalysis/issues/187
.. _DCDplugin: http://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/dcdplugin_8c-source.html#l00947
.. _wiki: https://github.com/MDAnalysis/mdanalysis/wiki/FileFormats#dcd


Writing out
===========

The writer follows recent NAMD/VMD convention for the unitcell (box lengths in Å and angle-cosines, ``[A, cos(gamma), B, cos(beta), cos(alpha), C]``). It writes positions in Å and time in AKMA time units.

Reading and writing these trajectories
within MDAnalysis will work seamlessly. However, if you process those trajectories
with other tools, you need to watch out that time and unitcell dimensions
are correctly interpreted.