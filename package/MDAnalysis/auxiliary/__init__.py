# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Auxiliary Readers --- :mod:`MDAnalysis.auxiliary`
=================================================

The auxiliary submodule contains code for reading 'auxiliary' data from a 
trajectory and allowing allignment with trajectory timesteps. Additional 
methods in :class:`MDAnalysis.coordinates.base.Reader` allow auxiliary data to 
be added and read alongside a trajectory.

Auxiliary data are timeseries accompanying a trajectory not stored/read with
the trajectory :class:`~MDAnalysis.coordinates.base.Reader`, and may be stored 
internally (e.g. in an array) or read for a file. In general, auxiliary data is 
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
   |                         |        |           | :class:`XVGReader` is the default     |
   |                         |        |           | reader for .xvg files.                |
   +-------------------------+--------+-----------+---------------------------------------+

.. _Auxiliary API:

Auxiliary API
-------------
Auxiliary readers inherit from the base 
:class:`~MDAnalysis.auxiliary.base.AuxReader`. In stand-alone use they 
allow iteration over each step in a set of auxiliary data::

    aux = MDAnalysis.auxiliary.xvg.XVGReader('auxdata.xvg')
    for step in aux:
        print step

A base :class:`~MDAnalysis.auxiliary.base.AuxFileReader` is also provided,
extending :class:`~MDAnalysis.auxiliary.base.AuxReader` with attributes/methods 
for when each step is to be read from a file in turn.

A :func:`MDAnalysis.auxiliary.core.get_auxreader_for` function is available to 
return an appropriate :class:`~MDAnalysis.auxiliary.base.AuxReader` instance for
a set of auxiliary data, guessed from the datatype/file extension.


Adding auxiliaries to trajectories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For alignment with trajectory data, auxiliary readers provide methods to 
assign each auxiliary step to the nearest trajectory timestep, read all steps
assigned to a trajectory timestep and calculate 'representative' value(s) of 
the auxiliary data for that timestep. 

'Assignment' of auxiliary steps to trajectory timesteps is determined from time 
of the auxiliary step and the ``dt`` and ``offset`` of the trajectory (provided 
through a :class:`~MDAnalysis.coordinates.base.Timestep` instance ''ts'' 
from the trajectory), as::

    frame = math.floor((time - offset+ts.dt/2)/ts.dt)

If there are no auxiliary steps assigned to a given timestep (or none within 
``cutoff``, if set), the representative 
value is set to ``np.nan``.

Auxiliary data may be added to a trajectory 
(:class:`MDAnalysis.coordinates.base.Reader` object) through the 
:meth:`~MDAnalysis.coordinates.base.ProtoReader.add_auxiliary` method. Auxiliary data
may be passed in as a AuxReader instance, or directly as e.g. a filename, in 
which case :func:`~MDAnalysis.auxiliary.core.get_auxreader_for` is used to 
guess an appropriate reader. e.g.::

    u = MDAnalysis.Universe(PDB, XTC)
    u.trajectory.add_auxiliary('pullforce', './pull_force.xvg')

The auxiliary data will be read appropriately as the trajectory frame is 
updated, and the representative auxiliary value(s) will be available as e.g.
``u.trajectory.ts.aux.pullforce``.

The trajectory :class:`~MDAnalysis.coordinates.base.Reader` methods 
:meth:`~MDAnalysis.coordinates.base.ProtoReader.next_as_aux` and
:meth:`~MDAnalysis.coordinates.base.ProtoReader.iter_as_aux` allow for movement
through only trajectory timesteps to which one or more steps from a given 
auxiliary have been assigned. This may be used when auxiliary data are less
frequent to avoid representative values set to ``np.nan``::

    u.trajectory.add_auxiliary('low_f', 'low_freq_aux_data.xvg')
    for ts in u.iter_as_aux():
        run_analysis(ts.aux.low_f) # do something with auxiliary data

If the auxiliary data are more frequent, 
:meth:`~MDAnalysis.coordinates.base.ProtoReader.next_as_aux` and
:meth:`~MDAnalysis.coordinates.base.ProtoReader.iter_as_aux` behave the same as :meth:`~MDAnalysis.coordinates.base.ProtoReader.next` and 
:meth:`~MDAnalysis.coordinates.base.ProtoReader.__iter__`.


.. _AuxStep API:

AuxStep class
~~~~~~~~~~~~~
An AuxStep instance holds the auxiliary data for the current step. It is
updated whenever the a new auxiliary step is read.

AuxStep classes are derived from the base class 
:class:`~MDAnalysis.auxiliary.base.AuxStep`. The appropriate AuxStep class for
a given auxiliary reader is identified using the `_Auxstep` attribute.


Attributes
..........
The following are inherited from :class:`~MDAnalysis.auxiliary.base.AuxStep`:

  ``step``
      Current auxiliary step (0-based).
  ``_data``
      All recorded data for the current step, as a numpy array.
  ``time``
      Time of current auxiliary step, determined from the ``_data`` using 
      ``_select_time()``/``time_selector``, or calculated using ``dt`` and 
      ``initial_time``.
  ``data``
      Auxiliary values of interest for the current step, determined from 
      ``_data`` using ``_select_data()``/``data_selector``.

The following are stored in AuxStep but can be acessed through the parent 
auxiliary reader (see :ref:`AuxReader API` below). 

  ``_dt``
      Change in time between auxiliary steps (in ps). If not specified, will
      attempt to determine from auxiliary data; otherwise defaults to 1 ps.
  ``_initial_time``
      Time of first auxiliary step (in ps). If not specified, will attempt to
      determine from auxiliary data; otherwise defaults to 0 ps.
  ``_time_selector``
      Selection key to get time from full set of auxilairy data read with each
      step(``_data``), if time selection is enabled by the reader (by defining 
      a ``_select_time` method). Type depends on the auxiliary format - e.g. 
      where data is stored in columns, time_selector may be an index of 'time' column.
      Default value is ``None``, in which case step time is calculated from ``dt``,
      ``initial_time`` and ``step``.
  ``_data_selector``
      Selection key(s) to get data of interest from full set of auxilairy data read with 
      each step (``_data``), if data selection is enabled by the reader (by 
      defining a ``_select_data`` method). As for ``time_selector``, type 
      depends on the auxiliary format.
      If ``None`` (default value), ``_data`` is returned.


Methods
.......
The following methods are inherited from :class:`~MDAnalaysis.auxiliary.base.AuxStep`:

  ``__init__(**kwargs)``
    Setup appropriate attributes based on *kwargs*.

The following must be defined by each AuxStep:

  ``_empty_data()``
    Return a np.array in the same format as ``data``, but with all values 
    ``np.nan``; used as the auxiliary value for a trajectory when no 
    auxiliary steps are assigned to the current frame.

To allow selection of time/data of interest from the full set of auxiliary
data, the following must be provided:

  ``_select_time(key)``
    Return the value indicated by *key* from ``_data`` (the full set of data read 
    in from the current step). Raise ``ValueError`` if *key* is not a valid 
    time selector for the auxiliary format.

  ``_select_data(key)``
    Return, as a ndarray, the value(s) indicated by *key* (may be e.g. a list of 
    multiple individual 'keys') from ``_data``. Raise ``ValueError`` if *key* is 
    not a valid data selector for the auxiliary format.


.. _AuxReader API:

AuxReader class
~~~~~~~~~~~~~~~

Registry
........

In order to facilitate guessing of appropriate AuxReaders, all AuxReaders
should set as appropriate `format` attribute. For files, this will be the
expected file extension.


Attributes
..........

The following attributes are inherited from 
:class:`~MDAnalaysis.auxiliary.base.AuxReader`:

  ``name``
      Name under which auxiliary data will be stored in trajectory/Timestep.
  ``represent_ts_as``
      Method to use in calculating representative auxiliary value for a 
      timestep. Default is 'closest'.
  ``cutoff``
      Cutoff (in ps) for ignoring auxiliary steps when calculating 
      representative value(s) for a timestep.
  ``auxstep``
      The :class:`~MDAnalysis.auxiliary.base.AuxStep` object to store data for
      the current step. Customised for each auxiliary format to allow selection
      of data.
  ``n_steps``
      Total number of auxiliary steps
  ``constant_dt``
      Boolean of whether dt is constant throughout auxiliary data. Default is 
      ``True``.
  ``frame_data``
      `data` from each auxiliary step assigned to the 
      last-read trajectory timestep.
  ``frame_rep``
      Represenatative value(s) of auxiliary data for last-read trajectory timestep.

The following are stored in ``auxstep`` but may be accessed from the auxiliary
reader; see the :ref:`AuxStep API` above.

  ``step``
  ``time``
  ``dt``
  ``initial_time``
  ``time_selector``
  ``data_selector``

:class:`~MDAnalaysis.auxiliary.base.AuxFileReader` additionally provides:
  ``auxfilename``
      Name of the auxiliary file.
  ``auxfile``
      File object for auxiliary file.

Each auxiliary reader class must also identify an appropriate 
:class:`~MDAnalysis.auxiliary.basae.AuxStep` with the `_Auxstep` attribute.
 

Methods
.......
The following methods are inherited from 
:class:`~MDAnalaysis.auxiliary.base.AuxReader`:

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

  ``rewind()``
    Reposition to first step.

  ``read_ts(ts)``
    Read auxiliary steps assigned to trajectory step `ts`; if ``name`` is set, 
    update the representative value stored in `ts`, and return.

  ``move_to_ts(ts)``
    Move to the auxiliary step before the trajectory timestep *ts*, such that
    calling 
    :meth:`~MDAnalysis.auxiliary.base.AuxFileReader.auxfile`_read_next_step` 
    reads the first step assigned to *ts* or, when auxiliary data are less 
    frequent and no steps are assigned to *ts*, the first step following *ts*.

  ``step_to_frame(step, ts)``
    Return the frame number of the trajectory described by `ts` to which 
    the auxiliary step `step` is assigned.

  ``step_to_time(step)``
    Return the time of the auxiliary step `step`. 

  ``calc_representative()``
    Return the representative value calculated from ``data`` following
    the method specified by ``represent_ts_as`` and ``cutoff``.

  ``__enter__()``
    Entry method for Context Manager (returns self).

  ``__exit__()``
    Exit method for Context Manager (calls ``close()``).

  ``__del__()``
    Calls ``close()``.


Each AuxReader must subclass :class:`~MDAnalysis.auxiliary.base.AuxReader`
and addionally define:

  ``__init__(auxdata, **kwargs)``
    Additional process *kwargs* and any necessary initial processing of 
    *auxdata*. Must :func:`super` through  
    :class:`~MDAnalysis.auxiliary.base.AuxReader`.

  ``_read_next_step()``
    Read the data from the next auxiliary step and set ``_data`` to the list of
    value(s) from each column. Raise `StopIteration` when attempting to read 
    past last step.

  ``_go_to_step(i)``
    Move to and read step `i` (0-based) from the auxiliary data. Raise 
    ValueError when i is out of range. 
    Sets ``_data`` in the AuxStep to an appropriate numpy array of the full set 
    of data read in, and return the AuxStep.


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

For convinience, when reading auxiliary data from a file, step at a time 
(rather than reading the full file at once), 
:class:`~MDAnalysis.auxiliary.base.AuxFileReader` may be used, which 
extends :class:`~MDAnalysis.auxiliary.base.AuxReader`, providing the
following (though these may be overwritten by subclasses as appropriate):
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
