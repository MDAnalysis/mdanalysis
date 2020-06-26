.. -*- coding: utf-8 -*-
.. _XYZ-format:

==============
XYZ trajectory
==============

.. include:: classes/XYZ.txt

The :ref:`XYZ format <XYZ-format>` is a loosely defined, simple
coordinate trajectory format. The implemented format definition was
taken from the `VMD xyzplugin`_ and is therefore compatible with VMD.

Reading in
==========

As XYZ files only have atom name information, the atoms are all assigned to the same residue and segment.

The default timestep in MDAnalysis is 1 ps. A different timestep can be defined :ref:`by passing in the dt argument to Universe <universe-kwargs>`.



XYZ specification
=================

Definiton used by the :class:`XYZReader` and :class:`XYZWriter` (and
the `VMD xyzplugin`_ from whence the definition was taken)::

    [ comment line            ] !! NOT IMPLEMENTED !! DO NOT INCLUDE
    [ N                       ] # of atoms, required by this xyz reader plugin  line 1
    [ molecule name           ] name of molecule (can be blank)                 line 2
    atom1 x y z [optional data] atom name followed by xyz coords                line 3
    atom2 x y z [ ...         ] and (optionally) other data.
    ...
    atomN x y z [ ...         ]                                                 line N+2


Note
----
* comment lines not implemented (do not include them)
* molecule name: the line is required but the content is ignored
  at the moment
* optional data (after the coordinates) are presently ignored


.. Links
.. _`VMD xyzplugin`: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html