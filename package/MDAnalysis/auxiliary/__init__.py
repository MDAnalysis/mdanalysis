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
trajectory and allowing allignment with existing trajectory data.

Auxiliary data are timeseries accompanying a trajectory not stored/read with
the trajectory readers, and may be stored internally (e.g. in an array) or read
for a file.

Supported formats
-----------------
Currently supported format:
 * XVG (file extension .xvg). Produced by Gromacs during simulation or analysis. 


AuxReaders
----------
Auxiliary readers inherit from the base 
:class:`~MDAnalysis.auxiliary.base.BaseAuxReader`. In stand-alone use they 
allow iteration over each step in a set of auxiliary data::

  aux = MDAnalysis.auxiliary.xvg.XVGReader('auxdata.xvg')
  for step_data in aux:
      print step_data

A :class:`~MDAnalysis.auxiliary.base.BaseAuxFileReader` is also provided,
extending the base AuxReader with attributes/methods when each step is to 
be read from a file in turn.

A :func:`MDAnalysis.auxiliary.core.get_auxreader_for` function is available to 
return an appropriate :class:`MDAnalysis.auxiliary.base.AuxReader` instance,
guessed based on datatype/file extension.

For reading alongside a trajectory, auxiliary readers provide methods to 
assign each auxiliary step to the nearest trajectory timestep, read all steps
assigned to a trajectory timestep and calculate 'representative' value(s) of 
the auxiliary data for that timestep. 

'Assignment' of auxiliary steps to trajectory timesteps is calculated from time 
of the auxiliary step and the *dt* and *offset* of the trajectory (provided 
through a :class:`~MDAnalysis.coordinates.base.Timestep` instance representing 
the trajectory), as:

  assigned_traj_frame = math.floor((aux_time - traj_offset+traj_dt/2)/traj_dt)

Auxiliary data may be added to a trajectory 
(:class:`MDAnalysis.coordinates.base.Reader` object) through the 
:meth:`MDAnalysis.coordinates.base.Reader.add_auxiliary` method. Auxiliary data
may be passed in as a AuxReader instance, or directly as e.g. a filename, in 
which case :func:`MDAnalysis.auxiliary.core.get_auxreader_for` is used to 
guess an appropriate reader. e.g.::

  u = MDAnalysis.Universe(PDB, XTC)
  u.trajectory.add_auxiliary('pullforce', './pull_force.xvg')

The auxiliary data will be read appropriately as the trajectory frame is 
updated, and the representative auxiliary value(s) will be available as e.g.
``u.trajectory.ts.aux.pullforce``.


Attributes
~~~~~~~~~~
In order to facilitate guessing of appropriate AuxReaders, all AuxReaders
should set as appropriate `format` attribute. For files, this will be the
expected file extension.

The following attributes are inherited from 
:class:`~MDAnalaysis.auxiliary.base.BaseAuxReader`:
  ``name``
      Name under which auxiliary data will be stored in the trajectory.
  ``represent_ts_as``
      Method to use in calculating representative auxiliary value for a 
      timestep. Default is 'closest'.
  ``cutoff``
      Cutoff (in ps) for ignoring auxiliary steps when calculating 
      representative value(s) for a timestep.
  ``dt``
      Change in time between auxiliary steps (in ps). If not specified, will
      attempt to determine from auxiliary data; otherwise defaults to 1ps.
  ``initial_time``
      Time of first auxilairy step (in ps). If not specified, will attempt to
      determine from auxiliary data; otherwise defaults to 0ps.
  ``time_col``
      Index of column in auxiliary data storing time (default value ``None``).
  ``data_cols``
      List of indices of columns containing data of interest. If not set,
      will default to all columns excluding that containing time (if present).
  ``constant_dt``
      Boolean of whether dt is constant throughout auxiliary data.
  ``step``
      Current auxiliary step, starting at 0.
  ``n_steps``
      Total number of auxiliary steps
  ``n_cols``
      Number of columns of data for each auxiliary step.
  ``_data``
      All recorded data for the current step (split into columns).
  ``time``
      Time of current auxiliary step, from the appropriate column of `_data` or
      calculated using `dt` and `initial_time`.
  ``times``
      List of the time for each auxiliary step.
  ``step_data``
      List of auxiliary values of interest for the current step, from the 
      appropriate column(s) of `_data`
  ``ts_data``
      Numpy array of `step_data` from each auxiliary step assigned to the 
      last-read trajectory timestep.
  ``ts_diff``
      List of difference in time between the last-read trajectory timestep
      and each auxiliary step assigned to it.
   ``ts_rep``
      Represenatative value of auxiliary data for last-read trajectory timestep.

If reading from file, AuxFileReader additionally provides:
  ``auxfilename``
      Name of the auxiliary file.
  ``auxfile``
      File object for auxiliary file.


Methods
~~~~~~~
The following methods are inherited from 
:class:`~MDAnalaysis.auxiliary.base.BaseAuxReader`:

  ``__init__(**kwargs)``
    Setup appropriate attributes based on *kwargs*.

  ``__len__()``
    Number of steps in auxiliary data.

  ``next()``
    Advance to next step.

  ``__iter__()``
    Allow iteration through each auxiliary step.

  ``go_to_first_step()``
    Reposition to first step.

  ``read_ts(ts)``
    Read auxiliary steps assigned to trajectory step `ts`, update the
    representative value stored in `ts` and return.

  ``step_to_frame(step, ts)``
    Return the frame number (of the trajectory described by `ts` to which 
    the auxiliary step `step` is assigned.

  ``calc_representative()``
    Return the representative value calculated from ``ts_data`` following
    the method specified by ``represent_ts_as`` and ``cutoff``.

  ``__enter__()``
    Entry method for Context Manager (returns self).

  ``__exit__()``
    Exit method for Context Manager (calls ``close()``).

  ``__del__()``
    Calls ``close()``.


Each AuxReader must subclass :class:`~MDAnalaysis.auxiliary.base.BaseAuxReader`
and addionally define:

  ``__init__(auxdata, **kwargs)
    Additional process *kwargs* and any necessary initial processing of 
    *auxdata*. Must :func:`super` through  
    :class:`~MDAnalaysis.auxiliary.base.BaseAuxReader`.

  ``_read_next_step()``
    Read and set ``_data`` to the data from the next auxiliary step. Used
    by ``next()`` and ``__iter__()``. Return ``_data`` or raise 
    :exc:`StopIteration` when attempting to read past last step.

  ``move_to_ts(ts)``
    Move to the auxiliary step before the trajectory timestep *ts*, such that
    calling next() reads the first step assigned to *ts* or, when auxiliary data
    are less frequent and no steps are assigned to *ts*, the first step 
    following *ts*.

  ``count_n_steps()``
    Return the total number of steps.

  ``read_all_times()``
    Set the ``times`` list.

  ``go_to_step(i)``
    Move to and read step i (0-based) from the auxiliary data. Return ``_data``
    for that step.

Depending on the format of the auxiliary data, it may also be necessary to 
overwrite the following:
  ``_restart()``
    Reposition before the first step

  ``close()``

For convinience, when reading auxiliary data from a file, step at a time 
(rather than reading the full file at once), 
:class:`~MDAnalaysis.auxiliary.base.BaseAuxFileReader` may be used, which 
extends :class:`~MDAnalaysis.auxiliary.base.BaseAuxReader`, providing the
following (though these may be overwritten by subclasses as appropriate):
  ``__init__(auxfile, **kwargs)``
    Open *auxfile* and any additional processing of **kwargs**.

  ``_restart()``, 
    Seek to the start of *auxfile*

  ``close()``
    Close *auxfile*

  ``move_to_ts(ts)``
    Iterate through all steps until last step before *ts* is reached.

  ``go_to_step(i)``
    Iterate through all steps until step i is reached.

  ``_reopen()``
    Close *auxfile* and reopen.

"""

# registry of auxiliary readers
_AUXREADERS = {}

from . import base
from . import XVG
