.. -*- coding: utf-8 -*-
.. _TRJ-format:

===========================================
TRJ, MDCRD, CRDBOX (AMBER ASCII trajectory)
===========================================

.. include:: classes/TRJ.txt

MDAnalysis supports reading of `AMBER ASCII trajectories`_ ("traj") and :ref:`binary trajectories <NCDF-format>` ("netcdf").

.. _`AMBER ASCII trajectories`: https://ambermd.org/FileFormats.php#trajectory

.. important::

   In the AMBER community, these trajectories are often saved with the
   suffix '.crd'. This extension conflicts with the CHARMM CRD
   format and MDAnalysis *will not correctly autodetect AMBER ".crd"
   trajectories*. Instead, explicitly provide the ``format="TRJ"``
   argument to :class:`~MDAnalysis.core.universe.Universe`::

     u = MDAnalysis.Universe("top.prmtop", "traj.crd", format="TRJ")

Reading in
==========

Units are assumed to be the following default AMBER units:

    * length: Angstrom
    * time: ps

.. rubric:: Limitations

* Periodic boxes are only stored as box lengths A, B, C in an AMBER
trajectory; the reader always assumes that these are orthorhombic
boxes.
* The trajectory does not contain time information so we simply set
the time step to 1 ps :ref:`(or the user could provide it with the dt argument) <universe-kwargs>`
* Trajectories with fewer than 4 atoms probably fail to be read (BUG).
* If the trajectory contains exactly *one* atom then it is always
assumed to be non-periodic (for technical reasons).
* Velocities are currently *not supported* as ASCII trajectories.
