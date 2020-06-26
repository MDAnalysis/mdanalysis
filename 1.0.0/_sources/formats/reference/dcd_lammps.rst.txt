.. -*- coding: utf-8 -*-
.. _LAMMPS-format:

========================================
DCD (Flexible LAMMPS trajectory)
========================================

.. include:: classes/LAMMPS.txt

LAMMPS can `write DCD`_ trajectories but unlike a `CHARMM trajectory`_
(which is often called a DCD, even though CHARMM itself calls them
"trj") the time unit is not fixed to be the AKMA_ time unit but can depend on
settings in LAMMPS. The most common case for biomolecular simulations
appears to be that the time step is recorded in femtoseconds (command
`units real`_ in the input file) and lengths in ångströms. Other cases
are unit-less Lennard-Jones time units.

This presents a problem for MDAnalysis, because it cannot autodetect
the unit from the file. By default, we assume that the unit for
length is the ångström and the unit for time is the femtosecond. If this is
not true, then the user *should supply the appropriate units* in the
keywords ``timeunit`` and/or ``lengthunit`` to :class:`~MDAnalysis.coordinates.LAMMPS.DCDWriter` and
:class:`~MDAnalysis.core.universe.Universe` (which then calls
:class:`~MDAnalysis.coordinates.LAMMPS.DCDReader`).

.. _LAMMPS: http://lammps.sandia.gov/
.. _write DCD: http://lammps.sandia.gov/doc/dump.html
.. _CHARMM trajectory: http://www.charmm.org/documentation/c36b1/dynamc.html#%20Trajectory
.. _AKMA: http://www.charmm.org/documentation/c36b1/usage.html#%20AKMA
.. _units real: http://lammps.sandia.gov/doc/units.html
.. _units command: http://lammps.sandia.gov/doc/units.html
.. _`Issue 64`: https://github.com/MDAnalysis/mdanalysis/issues/64
.. _`Issue 84`: https://github.com/MDAnalysis/mdanalysis/issues/84
.. _`LAMMPS dump format`: http://lammps.sandia.gov/doc/dump.html