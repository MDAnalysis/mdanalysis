.. -*- coding: utf-8 -*-
.. _auxiliary:

===============
Auxiliary files
===============


Auxiliary readers allow you to read in timeseries data accompanying a trajectory, that is not stored in the regular trajectory file. 

Supported formats
-----------------

   +-------------------------+--------+-----------+---------------------------------------+
   | Reader                  | Format | Extension | Remarks                               |
   |                         |        | (if file) |                                       |
   +=========================+========+===========+=======================================+
   | :class:`.XVGReader`     | XVG    | xvg       | Produced by Gromacs during simulation | 
   |                         |        |           | or analysis.                          |
   |                         |        | (default) | Reads full file on initialisation.    |
   +-------------------------+--------+-----------+---------------------------------------+
   | :class:`.XVGFileReader` | XVG-F  | xvg       | Alternate xvg file reader, reading    |
   |                         |        |           | each step from the file in turn for a |
   |                         |        |           | lower memory footprint.               |
   |                         |        |           | :class:`.XVGReader` is the default    |
   |                         |        |           | reader for .xvg files.                |
   +-------------------------+--------+-----------+---------------------------------------+

Reading data directly
---------------------

.. ipython:: python

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import XVG_BZ2  # cobrotoxin protein forces
    aux = mda.auxiliary.core.auxreader(XVG_BZ2)
    aux

In stand-alone use, an auxiliary reader allows you to iterate over each step in a set of auxiliary data.

.. ipython:: python

    for step in aux:
        print(step.data)

Use slicing to skip steps.

.. ipython:: python

    for step in aux[1:2]:
        print(step.time)

The :func:`~MDAnalysis.auxiliary.core.auxreader` function uses the :func:`~MDAnalysis.auxiliary.core.get_auxreader_for` to return an appropriate class. This can guess the format either from a filename, '

.. ipython:: python

    mda.auxiliary.core.get_auxreader_for(XVG_BZ2)

or return the reader for a specified format.

.. ipython:: python

    mda.auxiliary.core.get_auxreader_for(format='XVG-F')


Loading data into a Universe
----------------------------

Auxiliary data may be added to a trajectory Reader through the 
:meth:`~MDAnalysis.coordinates.base.ProtoReader.add_auxiliary` method. Auxiliary data
may be passed in as a AuxReader instance, or directly as e.g. a filename, in 
which case :func:`~MDAnalysis.auxiliary.core.get_auxreader_for` is used to 
guess an appropriate reader.

.. ipython:: python

    from MDAnalysis.tests.datafiles import PDB_xvf, TRR_xvf
    u = mda.Universe(PDB_xvf, TRR_xvf)
    u.trajectory.add_auxiliary('protein_force', XVG_BZ2)
    for ts in u.trajectory:
        print(ts.aux.protein_force)

Passing arguments to auxiliary data
-----------------------------------

For alignment with trajectory data, auxiliary readers provide methods to 
assign each auxiliary step to the nearest trajectory timestep, read all steps
assigned to a trajectory timestep and calculate 'representative' value(s) of 
the auxiliary data for that timestep. 

To set a timestep or ??

'Assignment' of auxiliary steps to trajectory timesteps is determined from the time 
of the auxiliary step, ``dt`` of the trajectory and time at the first frame of the 
trajectory. If there are no auxiliary steps assigned to a given timestep (or none within 
``cutoff``, if set), the representative value(s) are set to ``np.nan``.

Iterating over auxiliary data
-----------------------------

Auxiliary data may not perfectly line up with the trajectory, or have missing data. 

.. ipython:: python
    :okwarning:

    from MDAnalysis.tests.datafiles import PDB, TRR
    u_long = mda.Universe(PDB, TRR)
    u_long.trajectory.add_auxiliary('protein_force', XVG_BZ2, dt=200)
    for ts in u_long.trajectory:
        print(ts.time, ts.aux.protein_force[:4])

The trajectory :class:`~MDAnalysis.coordinates.base.ProtoReader` methods 
:meth:`~MDAnalysis.coordinates.base.ProtoReader.next_as_aux` and
:meth:`~MDAnalysis.coordinates.base.ProtoReader.iter_as_aux` allow for movement
through only trajectory timesteps for which auxiliary data is available.

.. ipython:: python
    
    for ts in u_long.trajectory.iter_as_aux('protein_force'):
        print(ts.time, ts.aux.protein_force[:4])

This may be used to 
avoid representative values set to ``np.nan``, particularly when auxiliary data
is less frequent.

Sometimes the auxiliary data is longer than the trajectory. 

.. ipython:: python
    :okwarning:

    u_short = mda.Universe(PDB)
    u_short.trajectory.add_auxiliary('protein_force', XVG_BZ2)
    for ts in u_short.trajectory:
        print(ts.time, ts.aux.protein_force)


In order to acess auxiliary values at every individual step, including those 
outside the time range of the trajectory, 
:meth:`~MDAnalysis.coordinates.base.ProtoReader.iter_auxiliary` allows iteration
over the auxiliary independent of the trajectory.

.. ipython:: python
    :okwarning:

    for step in u_short.trajectory.iter_auxiliary('protein_force'):
        print(step.data)


To iterate over only a certain section of the auxiliary:

.. ipython:: python
    :okwarning:

    for step in u_short.trajectory.iter_auxiliary('protein_force', start=1, step=2):
        # every 2nd step from 1
        print(step.time)

The trajectory remains unchanged, and the auxiliary will be returned to the current 
timestep after iteration is complete.

Accessing auxiliary attributes
------------------------------

To check the values of attributes of an added auxiliary, use 
:meth:`~MDAnalysis.coordinates.base.ProtoReader.get_aux_attribute`.

.. ipython:: python

    u.trajectory.get_aux_attribute('protein_force', 'dt')


If attributes are settable, they can be changed using
:meth:`~MDAnalysis.coordinates.base.ProtoReader.set_aux_attribute`.

.. ipython:: python

    u.trajectory.set_aux_attribute('protein_force', 'data_selector', [1])


The auxiliary may be renamed using ``set_aux_attribute`` with 'auxname', or more
directly by using :meth:`~MDAnalysis.coordinates.base.ProtoReader.rename_aux`.

.. ipython:: python

    u.trajectory.ts.aux.protein_force
    u.trajectory.rename_aux('protein_force', 'f')
    u.trajectory.ts.aux.f

Recreating auxiliaries
----------------------

To recreate an auxiliary, the set of attributes necessary to replicate it can 
first be obtained with :meth:`~MDAnalysis.auxiliary.base.AuxReader.get_description`. 
The returned dictionary can then be passed to 
:func:`~MDAnalysis.auxiliary.core.auxreader` to load a new copy of the 
original auxiliary reader.

.. ipython:: python

    description = aux.get_description()
    list(description.keys())
    del aux
    mda.auxiliary.core.auxreader(**description)

The 'description' of any or all the auxiliaries added to a trajectory can be
obtained using :meth:`~MDAnalysis.coordinates.base.ProtoReader.get_aux_descriptions`.

.. ipython:: python

    descriptions = u.trajectory.get_aux_descriptions(['f'])

To reload, pass the dictionary into :meth:`~MDAnalysis.coordinates.base.ProtoReader.add_auxiliary`.

.. ipython:: python

    u2 = mda.Universe(PDB, TRR)
    for desc in descriptions:
        u2.trajectory.add_auxiliary(**desc)
