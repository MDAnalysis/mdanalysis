# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
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

"""
Auxiliary Readers --- :mod:`MDAnalysis.auxiliary`
=================================================

The auxiliary submodule contains code for reading 'auxiliary' data from a 
trajectory and allowing alignment with trajectory timesteps. Additional 
methods in :class:`MDAnalysis.coordinates.base.ProtoReader` allow auxiliary data
to be added and read alongside a trajectory.

Auxiliary data are timeseries accompanying a trajectory not stored in the 
regular trajectory file to be read by the trajectory Reader. They may be stored 
internally (e.g. in an array) or read from a file. In general, auxiliary data is 
assumed to be time ordered and contain no duplicates.

Supported formats
-----------------
Currently supported formats:

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
   | :class:`.EDRReader`     | EDR    | edr       | Produced by Gromacs during simulation.|
   |                         |        |           | Reads full file on initialisation     |
   +-------------------------+--------+-----------+---------------------------------------+

.. _Auxiliary API:

Auxiliary API
-------------
Auxiliary readers inherit from the base 
:class:`~MDAnalysis.auxiliary.base.AuxReader`. In stand-alone use they 
allow iteration over each step in a set of auxiliary data::

    aux = MDAnalysis.auxiliary.XVG.XVGReader('auxdata.xvg')
    for auxstep in aux:
        print(auxstep)

To iterate over only certain sections of the auxiliary::

    for auxstep in aux[100:200]:  
        # only steps 100-200
        do_something(auxstep)

Or to skip steps::

    for auxstep in aux[100::10]:  
        # every 10th step starting at 100
        do_something(auxstep)

A base :class:`~MDAnalysis.auxiliary.base.AuxFileReader` is also provided,
extending :class:`~MDAnalysis.auxiliary.base.AuxReader` with attributes/methods 
for when auxiliary data is to be read from a file by keeping the file open and
reading steps one-at-time as needed.

A :func:`~MDAnalysis.auxiliary.core.get_auxreader_for` function is available to 
return an appropriate :class:`~MDAnalysis.auxiliary.base.AuxReader` class for
a supplied format or set of auxiliary data, guessing the format from the 
datatype/file extension::

  auxreader = MDAnalysis.auxiliary.core.get_auxreader_for('auxdata.xvg')
  # will return the default XVGReader class
  auxreader = MDAnalysis.auxiliary.core.get_auxreader_for(format='XVG-F')
  # will return the XVGFileReader class

To directly load an instance of the guessed auxiliary reader class given the
supplied auxdata and any additional auxiliary reader options, the 
function :func:`~MDAnalysis.auxiliary.core.auxreader` can be used::

  aux = MDAnalysis.auxiliary.auxreader('auxdata.xvg', dt=2)
  for auxstep in aux:
     do_something(auxstep)


Auxiliaries and trajectories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For alignment with trajectory data, auxiliary readers provide methods to 
assign each auxiliary step to the nearest trajectory timestep, read all steps
assigned to a trajectory timestep and calculate 'representative' value(s) of 
the auxiliary data for that timestep. 

'Assignment' of auxiliary steps to trajectory timesteps is determined from the time 
of the auxiliary step, ``dt`` of the trajectory and time at the first frame of the 
trajectory (obtained through a :class:`~MDAnalysis.coordinates.timestep.Timestep` 
instance from the trajectory), as::

    frame = floor((time_at_step - time_at_frame_0 + dt/2)/dt)

If there are no auxiliary steps assigned to a given timestep (or none within 
``cutoff``, if set), the representative value(s) are set to ``np.nan``.

Adding an auxiliary to a trajectory
...................................
Auxiliary data may be added to a trajectory Reader through the 
:meth:`~MDAnalysis.coordinates.base.ProtoReader.add_auxiliary` method. Auxiliary data
may be passed in as a AuxReader instance, or directly as e.g. a filename, in 
which case :func:`~MDAnalysis.auxiliary.core.get_auxreader_for` is used to 
guess an appropriate reader. e.g.::

    u = MDAnalysis.Universe(PDB, XTC)
    u.trajectory.add_auxiliary('./pull_force.xvg', 'pullforce')

As the trajectory frame is updated, the auxiliary data will be read 
correspondingly, 
updated, and the representative auxiliary value(s) will be available as e.g.
``u.trajectory.ts.aux.pullforce``.


Iterating by an auxiliary
.........................
The trajectory :class:`~MDAnalysis.coordinates.base.ProtoReader` methods 
:meth:`~MDAnalysis.coordinates.base.ProtoReader.next_as_aux` and
:meth:`~MDAnalysis.coordinates.base.ProtoReader.iter_as_aux` allow for movement
through only trajectory timesteps to which one or more steps that fall within 
``cutoff`` from a given auxiliary have been assigned. This may be used to 
avoid representative values set to ``np.nan``, particularly when auxiliary data
is less frequent::

    u.trajectory.add_auxiliary('low_freq_aux_data.xvg', 'low_f')
    for ts in u.trajectory.iter_as_aux('low_f'):
        do_something(ts.aux.low_f) # do something with 'low_f' auxiliary data without
                                   # worrying about having to deal with np.nan

If the auxiliary data are more frequent and the cutoff (if set) is 
sufficiently high, :meth:`~MDAnalysis.coordinates.base.ProtoReader.next_as_aux` 
and :meth:`~MDAnalysis.coordinates.base.ProtoReader.iter_as_aux` behave the same 
as the :class:`~MDAnalysis.coordinates.base.ProtoReader` 
:meth:`~MDAnalysis.coordinates.base.ProtoReader.next` and 
:meth:`~MDAnalysis.coordinates.base.ProtoReader.__iter__`.

In order to access auxiliary values at every individual step, including those 
outside the time range of the trajectory, 
:meth:`~MDAnalysis.coordinates.base.ProtoReader.iter_auxiliary` allows iteration
over the auxiliary independent of the trajectory::

    for auxstep in u.trajectory.iter_auxiliary('pullforce'):
        do_something(auxstep)

To iterate over only a certain section of the auxiliary::

    for auxstep in u.trajectory.iter_auxiliary('pullforce', start=100, step=10):
        # every 10th step from 100
        do_something(auxstep)

The trajectory remains unchanged, and the auxiliary will be returned to the current 
timestep after iteration is complete.

Accessing auxiliary attributes
..............................
To check the values of attributes of an added auxiliary, use 
:meth:`~MDAnalysis.coordinates.base.ProtoReader.get_aux_attribute`, e.g.::

    pullf_dt = u.trajectory.get_aux_attribute('pullforce', 'dt')

If attributes are settable, they can be changed using
:meth:`~MDAnalysis.coordinates.base.ProtoReader.set_aux_attribute`, e.g.::

    u.trajectory.set_aux_attribute('pullforce', 'data_selector', [1])

The auxiliary may be renamed using ``set_aux_attribute`` with 'auxname', or more
directly by using :meth:`~MDAnalysis.coordinates.base.ProtoReader.rename_aux`::

    u.trajectory.rename_aux('pullforce', 'pullf')


Recreating auxiliaries
~~~~~~~~~~~~~~~~~~~~~~~
To recreate an auxiliary, the set of attributes necessary to replicate it can 
first be obtained with :meth:`~MDAnalysis.auxiliary.base.AuxReader.get_description`. 
The returned dictionary can then be passed to 
:func:`~MDAnalysis.auxiliary.core.auxreader` to load a new copy of the 
original auxiliary reader::

    description = aux.get_description()
    del(aux)
    reload_aux = MDAnalysis.auxiliary.auxreader(**description)

The 'description' of any or all the auxiliaries added to a trajectory can be
obtained using :meth:`~MDAnalysis.coordinates.base.ProtoReader.get_aux_descriptions`::

    descriptions = u.trajectory.get_aux_descriptions()

Get descriptions for selected auxiliaries only::

    descriptions = u.trajectory.get_aux_descriptions(['pullf', 'pullx'])

And to reload::

    for descr in descriptions:
        new_u.new_trajectory.add_auxiliary(**descr)


.. _AuxStep API:

AuxStep class
~~~~~~~~~~~~~
An AuxStep instance holds the auxiliary data for the current step. It is
updated whenever a new auxiliary step is read.

AuxStep classes are derived from the base class 
:class:`~MDAnalysis.auxiliary.base.AuxStep`. The appropriate AuxStep class for
a given auxiliary reader is identified in the reader by the `_Auxstep` attribute.


Attributes
..........
The following are inherited from :class:`~MDAnalysis.auxiliary.base.AuxStep`:

  ``step``
      Current auxiliary step (0-based).
  ``_data``
      All recorded data for the current step, as a numpy array.
  ``time``
      Time of current auxiliary step, as a float (in ps). Determined from the 
      ``_data`` if time selection is enabled and a valid ``time_selector`` 
      provided; otherwise calculated using ``dt`` and ``initial_time``.
  ``data``
      Auxiliary values of interest for the current step, as a numpy array. 
      Determined from ``_data`` id data selection is enabled at a valid 
      ``data_selector`` is provided; otherwise set to ``_data``.

The following are stored in AuxStep. The parent 
auxiliary reader has the equivalent attributes, though non-private (no _)
(see :ref:`AuxReader API` below). 

  ``_dt``
      Change in time between auxiliary steps (in ps). If not specified, will
      attempt to determine from auxiliary data; otherwise defaults to 1 ps.
  ``_initial_time``
      Time of first auxiliary step (in ps). If not specified, will attempt to
      determine from auxiliary data; otherwise defaults to 0 ps.
  ``_time_selector``
      Selection key to get time from full set of auxiliary data read with each
      step(``_data``) (ignored if time selection is not enabled by the reader). 
      Type depends on the auxiliary format - e.g. where data is stored in 
      columns, time_selector may be an index of 'time' column. 
      Default value is ``None``, in which case step time is calculated from
      ``dt``, ``initial_time`` and ``step``.
  ``_data_selector``
      Selection key(s) to get data of interest from full set of auxiliary data 
      read with each step (``_data``) (ignored if data selection is not enabled
      by the reader). As for ``time_selector``, type depends on the auxiliary 
      format. If ``None`` (default value), ``_data`` is returned.
  ``_constant_dt``
      Boolean of whether dt is constant throughout auxiliary data. Default is 
      ``True``.

Methods
.......
The following methods are inherited from 
:class:`~MDAnalysis.auxiliary.base.AuxStep`:

  ``__init__(**kwargs)``
    Setup appropriate attributes based on *kwargs*.

To enabled time/data selection in a particular AusStep, the following must be
provided:

  ``_select_time(key)``
    Return, as a float, the time value indicated by *key* from ``_data`` (the full 
    set of data read in from the current step). Raise ``ValueError`` if *key* is 
    not a valid time selector for the auxiliary format.

  ``_select_data(key)``
    Return, as a ndarray, the value(s) indicated by *key* (may be e.g. a list of 
    multiple individual 'keys') from ``_data``. Raise ``ValueError`` if *key* is 
    not a valid data selector for the auxiliary format.

Depending on the format of the auxiliary and hence the format of ``data``, it 
may be necessary to overload the following for a particular AuxStep:

  ``_empty_data()``
    Return a np.array in the same format as ``data`` with all values set to
    ``np.nan``. Used as the representative auxiliary value for a trajectory when 
    no auxiliary steps are assigned to the current frame.


.. _AuxReader API:

AuxReader class
~~~~~~~~~~~~~~~

Registry
........

In order to facilitate guessing of appropriate AuxReaders, all AuxReaders
should set as appropriate `format` attribute. For files, this will be the
expected file extension (in all caps).


Replicating auxiliaries
.......................

The names of the necessary attributes for replicating an auxiliary are stored in 
`required_attrs`, initially set in the base AuxReader. If a particular
AuxReader introduces additional attributes required to reload an auxiliary, 
these should be added.


Attributes
..........

The following attributes are inherited from 
:class:`~MDAnalysis.auxiliary.base.AuxReader`:

  ``auxname``
      Name under which auxiliary data will be stored in trajectory.
  ``represent_ts_as``
      Method to use in calculating representative auxiliary value for a 
      timestep. Default is 'closest'.
  ``cutoff``
      Cutoff (in ps) for ignoring auxiliary steps when calculating
      representative value(s) for a timestep.
  ``auxstep``
      The :class:`~MDAnalysis.auxiliary.base.AuxStep` object to store data for
      the current step. The particular AuxStep used will depend on the auxiliary
      format.
  ``n_steps``
      Total number of auxiliary steps.
  ``frame_data``
      `data` from each auxiliary step assigned to the last-read trajectory
      timestep.
  ``frame_rep``
      Representative value(s) of auxiliary data for last-read trajectory
      timestep.

The following are stored in ``auxstep`` as private attributes (with _)
but may be accessed from the auxiliary reader; see the :ref:`AuxStep API` above.

  ``step``

  ``time``

  ``dt``

  ``initial_time``

  ``time_selector``

  ``data_selector``

  ``constant_dt``

:class:`~MDAnalysis.auxiliary.base.AuxFileReader` provides:

  ``auxfile``
      File object for the auxiliary file.

Each auxiliary reader subclass will additionally set:

   ``_auxdata``
       Value of 'auxdata' used to load the auxiliary - e.g. path to the file 
       containing the auxiliary data. Will be used when getting the 'description'
       of the auxiliary for recreating.

Each auxiliary reader class must also identify an appropriate 
:class:`~MDAnalysis.auxiliary.base.AuxStep` with the `_Auxstep` class attribute
 

Methods
.......
The following methods are inherited from 
:class:`~MDAnalysis.auxiliary.base.AuxReader`:

  ``__init__(**kwargs)``
    Setup appropriate attributes based on *kwargs*.

  ``__len__()``
    Number of steps in auxiliary data.

  ``next()``
    Advance to next step.

  ``__iter__()``
    Allow iteration through each auxiliary step.

  ``__getitem__(step)``
    If *step* is a single index, move to that step and return the :class:`AuxStep`;
    if a list or slice, return an iterator over the specified auxiliary steps.
    See examples in :ref:`Auxiliary API` above.

  ``rewind()``
    Reposition to first step.

  ``update_ts(ts)``
    Read auxiliary steps assigned to trajectory step `ts`, calculate and update 
    the representative value in `ts`. Return `ts`.

  ``read_ts(ts)``
    Read auxiliary steps assigned to trajectory step `ts` and calculate
    representative value.

  ``move_to_ts(ts)``
    Move to the auxiliary step before the trajectory timestep *ts*, such that
    calling ``_read_next_step`` reads the first step assigned to *ts* or, when 
    auxiliary data are less frequent and no steps are assigned to *ts*, the 
    first step following *ts*.

  ``step_to_frame(step, ts, return_time_diff=False)``
    Return the frame number of the trajectory described by `ts` to which 
    the auxiliary step `step` is assigned. Optionally also return the difference
    it time between the step and the returned frame.

  ``step_to_time(step)``
    Return the time of the auxiliary step `step`. 

  ``next_nonempty_frame(ts)``
    Return the frame number of the next trajectory frame (after the current
    auxiliary time) for which a representative auxiliary value can be
    calculated (i.e., there is at least one assigned auxiliary step within
    ``cutoff``).

  ``calc_representative()``
    Return the representative value calculated from ``data`` following
    the method specified by ``represent_ts_as`` and ``cutoff``.

  ``__enter__()``
    Entry method for Context Manager (returns self).

  ``__exit__()``
    Exit method for Context Manager (calls ``close()``).

  ``__del__()``
    Calls ``close()``.

  ``get_description``
    Get the values of the attributes required for replicating an auxiliary (as
    listed in ``required_attrs``) and return as a dictionary.

  ``__eq__``
    Check for equality by checking each of the attributes required for
    replicating an auxiliary (as listed in ``required_attrs``) are equal.

Each AuxReader must subclass :class:`~MDAnalysis.auxiliary.base.AuxReader`
and addionally define:

  ``__init__(auxdata, **kwargs)``
    Additional processing of *kwargs* and any necessary initial processing of 
    *auxdata*. Must :func:`super` through
    :class:`~MDAnalysis.auxiliary.base.AuxReader`.

  ``_read_next_step()``
    Read the data from the next auxiliary step and update ``auxstep`` as
    appropriate. Raise `StopIteration` when attempting to read past last step.

  ``_go_to_step(i)``
    Move to and read step `i` (0-based) from the auxiliary data. Raise 
    ValueError when i is out of range. Update the ``auxstep`` and return.


Depending on the format of the auxiliary data, it may also be necessary to 
define/overwrite the following:

  ``read_all_times()``
    Return the list of times for each step (only required if ``constant_dt``
    may be false).

  ``count_n_steps()``
    Return the total number of steps (only required if `_n_steps` not otherwise
    set during `__init__`).

  ``_restart()``
    Reposition before the first step

  ``close()``

For convenience, when reading auxiliary data from an open file, one step at a
time, :class:`~MDAnalysis.auxiliary.base.AuxFileReader`
extends :class:`~MDAnalysis.auxiliary.base.AuxReader` by providing the
following (these may be overwritten by subclasses as appropriate):

  ``__init__(auxfile, **kwargs)``
    Open *auxfile* and any additional processing of *kwargs*.

  ``_restart()``,
    Seek to the start of ``auxfile``

  ``close()``
    Close ``auxfile``

  ``_go_to_step(i)``
    Iterate through all steps until step `i` is reached.

  ``_reopen()``
    Close ``auxfile`` and reopen.


"""

# registry of auxiliary readers
_AUXREADERS = {}

from . import base
from . import XVG
from . import EDR
from .core import auxreader
