# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2020 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""HOLE Analysis --- :mod:`MDAnalysis.analysis.hole2`
=====================================================================================

:Author: Lily Wang
:Year: 2020
:Copyright: GNU Public License v3

.. versionadded:: 1.0.0

This module contains the tools to interface with HOLE_
:cite:p:`Smart1993,Smart1996` to analyse an ion channel pore or transporter
pathway :cite:p:`Stelzl2014`.

Using HOLE on a PDB file
------------------------

Use the :func:``hole`` function to run `HOLE`_ on a single PDB file. For example,
the code below runs the `HOLE`_ program installed at '~/hole2/exe/hole' ::

    from MDAnalysis.tests.datafiles import PDB_HOLE
    from MDAnalysis.analysis import hole2
    profiles = hole2.hole(PDB_HOLE, executable='~/hole2/exe/hole')
    # to create a VMD surface of the pore
    hole2.create_vmd_surface(filename='hole.vmd')

``profiles`` is a dictionary of HOLE profiles, indexed by the frame number. If only
a PDB file is passed to the function, there will only be one profile at frame 0.
You can visualise the pore by loading your PDB file into VMD, and in
Extensions > Tk Console, type::

    source hole.vmd

You can also pass a DCD trajectory with the same atoms in the same order as
your PDB file with the ``dcd`` keyword argument. In that case, ``profiles`` will
contain multiple HOLE profiles, indexed by frame.

The HOLE program will create some output files:

    * an output file (default name: hole.out)
    * an sphpdb file (default name: hole.sph)
    * a file of van der Waals' radii
      (if not specified with ``vdwradii_file``. Default name: simple2.rad)
    * a symlink of your PDB or DCD files (if the original name is too long)
    * the input text (if you specify ``infile``)

By default (`keep_files=True`), these files are kept. If you would like to
delete the files after the function is wrong, set `keep_files=False`. Keep in
mind that if you delete the sphpdb file, you cannot then create a VMD surface.


Using HOLE on a trajectory
--------------------------

You can also run HOLE on a trajectory through the :class:`HoleAnalysis` class.
This behaves similarly to the ``hole`` function, although arguments such as ``cpoint``
and ``cvect`` become runtime arguments for the :meth:`~HoleAnalysis.run` function.

The class can be set-up and run like a normal MDAnalysis analysis class::

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import MULTIPDB_HOLE
    from MDAnalysis.analysis import hole2

    u = mda.Universe(MULTIPDB_HOLE)

    ha = hole2.HoleAnalysis(u, executable='~/hole2/exe/hole') as h2:
    ha.run()
    ha.create_vmd_surface(filename='hole.vmd')

The VMD surface created by the class updates the pore for each frame of the trajectory.
Use it as normal by loading your trajectory in VMD and sourcing the file in the Tk Console.

You can access the actual profiles generated in the ``results`` attribute::

    print(ha.results.profiles)

Again, HOLE writes out files for each frame. If you would like to delete these files
after the analysis, you can call :meth:`~HoleAnalysis.delete_temporary_files`::

    ha.delete_temporary_files()

Alternatively, you can use HoleAnalysis as a context manager that deletes temporary
files when you are finished with the context manager::

    with hole2.HoleAnalysis(u, executable='~/hole2/exe/hole') as h2:
        h2.run()
        h2.create_vmd_surface()


Using HOLE with VMD
-------------------

The :program:`sos_triangle` program that is part of HOLE_ can write an input
file for VMD_ to display a triangulated surface of the pore found by
:program:`hole`. This functionality is available with the
:meth:`HoleAnalysis.create_vmd_surface` method
[#create_vmd_surface_function]_. For an input trajectory MDAnalysis writes a
*trajectory* of pore surfaces that can be animated in VMD together with the
frames from the trajectory.


Analyzing a full trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~

To analyze a full trajectory and write pore surfaces for all frames to file
:file:`hole_surface.vmd`, use ::

    import MDAnalysis as mda
    from MDAnalysis.analysis import hole2

    # load example trajectory MULTIPDB_HOLE
    from MDAnalysis.tests.datafiles import MULTIPDB_HOLE

    u = mda.Universe(MULTIPDB_HOLE)

    with hole2.HoleAnalysis(u, executable='~/hole2/exe/hole') as h2:
        h2.run()
        h2.create_vmd_surface(filename="hole_surface.vmd")

In VMD, load your trajectory and then in the tcl console
(e.g.. :menuselection:`Extensions --> Tk Console`) load the surface
trajectory:

.. code-block:: tcl

   source hole_surface.vmd

If you only want to *subsample the trajectory* and only show the surface at
specific frames then you can either load the trajectory with the same
subsampling into VMD or create a subsampled trajectory.


Creating subsampled HOLE surface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For example, if we want to start displaying at frame 1 (i.e., skip frame 0), stop at frame 7, and
only show every other frame (step 2) then the HOLE analysis will be ::

    with hole2.HoleAnalysis(u, executable='~/hole2/exe/hole') as h2:
        h2.run(start=1, stop=9, step=2)
        h2.create_vmd_surface(filename="hole_surface_subsampled.vmd")

The commands produce the file ``hole_surface_subsampled.vmd`` that can be loaded into VMD.

.. Note::

   Python (and MDAnalysis) stop indices are *exclusive* so the parameters
   ``start=1``, ``stop=9``, and ``step=2`` will analyze frames 1, 3, 5, 7.

.. _Loading-a-trajectory-into-VMD-with-subsampling:

Loading a trajectory into VMD with subsampling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Load your system into VMD. This can mean to load the topology file with
:menuselection:`File --> New Molecule` and adding the trajectory with
:menuselection:`File --> Load Data into Molecule` or just :menuselection:`File
--> New Molecule`.

When loading the trajectory, subsample the frames by setting parametes in in
the :guilabel:`Frames` section. Select *First: 1*, *Last: 7*, *Stride: 2*. Then
:guilabel:`Load` everything.

.. Note::

   VMD considers the stop/last frame to be *inclusive* so you need to typically
   choose one less than the ``stop`` value that you selected in MDAnalysis.

Then load the surface trajectory:

.. code-block:: tcl

   source hole_surface_subsampled.vmd

You should see a different surface for each frame in the trajectory. [#vmd_extra_frame]_


Creating a subsampled trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Instead of having VMD subsample the trajectory as described in
:ref:`Loading-a-trajectory-into-VMD-with-subsampling` we can write a subsampled
trajectory to a file. Although it requires more disk space, it can be
convenient if we want to visualize the system repeatedly.

The example trajectory comes as a multi-PDB file so we need a suitable topology
file. If you already have a topology file such as a PSF, TPR, or PRMTOP file
then skip this step. We write frame 0 as a PDB :file:`frame0.pdb` (which we
will use as the topology in VMD)::

    u.atoms.write("frame0.pdb")

Then write the actual trajectory in a convenient format such as TRR (or
DCD). Note that we apply the same slicing (``start=1``, ``stop=9``, ``step=2``)
to the trajectory itself and then use it as the value for the ``frames``
parameter of :meth:`AtomGroup.write<MDAnalysis.core.groups.AtomGroup.write>`
method::

    u.atoms.write("subsampled.trr", frames=u.trajectory[1:9:2])

This command creates the subsampled trajectory file :file:`subsampled.trr` in
TRR format.

In VMD we load the topology and the trajectory and then load the surface. In
our example we have a PDB file (:file:`frame0.pdb`) as topology so we need to
remove the first frame [#vmd_extra_frame]_ (skip the "trim" step below if you
are using a true topology file such as PSF, TPR, or PRMTOP). To keep this
example compact, we are using the tcl command line interface in VMD_
(:menuselection:`Extensions --> Tk Console`) for loading and trimming the
trajectory; you can use the menu commands if you prefer.

.. code-block:: tcl

   # load topology and subsampled trajectory
   mol load pdb frame0.pdb trr subsampled.trr

   # trim first frame (frame0) -- SKIP if using PSF, TPR, PRMTOP
   animate delete beg 0 end 0

   # load the HOLE surface trajectory
   source hole_surface_subsampled.vmd

You can now animate your molecule together with the surface and render it.


.. _HOLE: http://www.holeprogram.org
.. _VMD: https://www.ks.uiuc.edu/Research/vmd/

Functions and classes
---------------------

.. autofunction:: hole

.. autoclass:: HoleAnalysis
   :members:


.. rubric:: References

.. bibliography::
    :filter: False
    :style: MDA

    Smart1993
    Smart1996
    Stelzl2014

.. rubric:: Footnotes

.. Footnotes

.. [#create_vmd_surface_function] If you use the :class:`hole` class to run
              :program:`hole` on a single PDB file then you can use
              :func:`MDAnalysis.analysis.hole2.utils.create_vmd_surface`
              function to manually run :program:`sph_process` and
              :program:`sos_triangle` on the output files andcr eate a surface
              file.

.. [#vmd_extra_frame] If you loaded your system in VMD_ from separate topology
              and trajectory files and the topology file contained coordinates
              (such as a PDB or GRO) file then your trajectory will have an
              extra initial frame containing the coordinates from your topology
              file. Delete the initial frame with :menuselection:`Molecule -->
              Delete Frames` by setting *First* to 0 and *Last* to 0 and
              selecting :guilabel:`Delete`.

.. [#HOLEDCD] PDB files are not the only files that :program:`hole` can
              read. In principle, it is also able to read CHARMM DCD
              trajectories and generate a hole profile for each frame. However,
              native support for DCD in :program:`hole` is patchy and not every
              DCD is recognized. In particular, At the moment, DCDs generated
              with MDAnalysis are not accepted by HOLE. To overcome this
              PDB / DCD limitation, use :class:`HoleAnalysis` which creates
              temporary PDB files for each frame of a
              :class:`~MDAnalysis.core.universe.Universe` or
              :class:`~MDAnalysis.core.universe.AtomGroup` and runs
              :func:``hole`` on each of them.

"""
from ...due import due, Doi
from .hole import hole, HoleAnalysis
from .utils import create_vmd_surface

due.cite(Doi("10.1016/S0006-3495(93)81293-1"),
         description="HOLE program",
         path="MDAnalysis.analysis.hole2",
         cite_module=True)
due.cite(Doi("10.1016/S0263-7855(97)00009-X"),
         description="HOLE program",
         path="MDAnalysis.analysis.hole2",
         cite_module=True)
due.cite(Doi("10.1016/j.jmb.2013.10.024"),
         description="HOLE trajectory analysis with orderparameters",
         path="MDAnalysis.analysis.hole2",
         cite_module=True)
del Doi
