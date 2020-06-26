.. -*- coding: utf-8 -*-
.. _units:

===================
Units and constants
===================

The units of MDAnalysis trajectories are the **Å** (**ångström**) for
**length** and **ps** (**picosecond**) for **time**. Regardless of how the 
original MD format stored the trajectory data, MDAnalysis converts it to 
MDAnalysis units when reading the data in, and converts back when writing the 
data out. Analysis classes generally also use these default units. Exceptions 
to the default units are always noted in the documentation; for example, mass
densities can be given in :math:`g/cm^3`. 

Other base units are listed in the table :ref:`table-baseunits`.

.. _table-baseunits:

.. Table:: Base units in MDAnalysis

   =========== ============== ===============================================
   Quantity    Unit            SI units
   =========== ============== ===============================================
   length       Å              :math:`10^{-10}` m
   time         ps             :math:`10^{-12}` s
   energy       kJ/mol         :math:`1.66053892103219 \times 10^{-21}` J
   charge       :math:`e`      :math:`1.602176565 \times 10^{-19}` As
   force        kJ/(mol·Å)     :math:`1.66053892103219 \times 10^{-11}` J/m
   speed        Å/ps           :math:`100` m/s
   mass         u              :math:`1.66053906660(50) \times 10^{-27}` kg
   angle        degrees        :math:`\frac{\pi}{180}` rad
   =========== ============== ===============================================


Unit conversion
===============

Quantities can be converted from units with :func:`~MDAnalysis.units.convert`.
:func:`~MDAnalysis.units.convert` simply multiplies the initial quantity with a 
precomputed conversion factor, as obtained from 
:func:`~MDAnalysis.units.get_conversion_factor`. 

The computed conversion factors for each quantity type are stored in :mod:`MDAnalysis.units` and shown below.

Constants
=========

.. include:: generated/units_table.txt
