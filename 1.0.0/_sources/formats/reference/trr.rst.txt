.. -*- coding: utf-8 -*-
.. _TRR-format:

======================================
TRR (GROMACS lossless trajectory file)
======================================

.. include:: classes/TRR.txt

The GROMACS TRR trajectory is a lossless format. This file format can store coordinates, velocities, and forces. 

.. important::

    MDAnalysis currently treats trajectories with damaged frames by truncating them at the frame before. Check that you are loading a valid file with `gmx check`_.

.. _`gmx check` : http://manual.gromacs.org/documentation/current/onlinehelp/gmx-check.html


Reading in
==========

MDAnalysis uses XDR based readers for GROMACS formats, which store offsets on the disk. The offsets are used to enable access to random frames efficiently. These offsets will be generated automatically the first time the trajectory is opened, and offsets are generally stored in hidden ``*_offsets.npz`` files. [#error]_ 

Trajectories split across multiple files can be :ref:`read continuously into MDAnalysis <chainreader>` with ``continuous=True``, in the style of `gmx trjcat`_.

.. _`gmx trjcat`: http://manual.gromacs.org/documentation/2018/onlinehelp/gmx-trjcat.html

Writing out
===========

If the ``data`` dictionary of a :class:`~MDAnalysis.coordinates.base.Timestep` contains a ``lambda`` value, this will be used for the written TRR file. Otherwise, lambda is set to 0.

Developer notes
===============

It sometimes can happen that the stored offsets get out off sync with the trajectory they refer to. For this the offsets also store the number of atoms, size of the file and last modification time. If any of them change, the offsets are recalculated. Writing of the offset file can fail when the directory where the trajectory file resides is not writable or if the disk is full. In this case a warning message will be shown but the offsets will nevertheless be used during the lifetime of the trajectory Reader. However, the next time the trajectory is opened, the offsets will have to be rebuilt again.

.. [#error] Occasionally, MDAnalysis fails to read XDR offsets, resulting in an error. The workaround for this is to create the Universe with regenerated offsets by using the keyword argument ``refresh_offsets=True``, as documented in `Issue 1893`_.

.. _`Issue 1893`: https://github.com/MDAnalysis/mdanalysis/issues/1893