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
Auxiliary Readers --- :mod:`MDAnalysis.auxiliary.base`
======================================================

Base classes for deriving all auxiliary data readers. See the API in :mod:`~MDAnalysis.auxiliary.__init__`.

.. autoclass:: AuxReader
   :members:

.. autoclass:: AuxFileReader
   :members:

"""

import six

import numpy as np
import math
import warnings

from ..lib.util import asiterable

from . import _AUXREADERS


class _AuxReaderMeta(type):
    # auto register on class creation
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            fmt = asiterable(classdict['format'])
        except KeyError:
            pass
        else:
            for f in fmt:
                _AUXREADERS[f] = cls


class AuxReader(six.with_metaclass(_AuxReaderMeta)):
    """ Base class for auxiliary readers.

    Allows iteration over a set of data from a trajectory, additional 
    ('auxiliary') to the regular positions/velocities/etc. This auxiliary 
    data may be stored in e.g. an array or a separate file.

    Auxiliary data may be added to a trajectory by 
    :meth:`MDAnalysis.coordinates.base.Reader.add_auxiliary`, passing either an 
    AuxReader instance or the data/filename, in which case an appropriate reader 
    will be selected based on type/file extension.
   
    The AuxReader will handle alignment of the auxiliary and trajectory 
    timesteps (auxiliary data steps are 'assigned' to the closest trajectory 
    timestep), and for each trajectory timestep will provide a 'representative'
    auxiliary value (or values) based on the steps assigned to that timestep.


    Paramaters
    ----------
    name : str, optional
        Name for auxiliary data. When added to a trajectory, the representative 
        auxiliary value(s) for the timestep are stored as ``ts.aux.name``.
    represent_ts_as : {'closest', 'average'}
        How to calculated representative value of auxiliary data for a 
        trajectory timestep. Currently available:
          *'closest': value from step closest to the trajectory timestep
          *'average': average of values from auxiliary steps assigned to 
                      the trajectory timestep.
    cutoff : float, optional
        Auxiliary steps further from the trajectory timestep than *cutoff* 
        will be ignored when calculating representative values (the default
        value is -1, which indicates all auxiliary steps assigned to that 
        timestep will be used).
    dt : float, optional
        Change in time between auxiliary steps (in ps). If not specified, will
        attempt to be determined from auxiliary data; otherwise defaults to 1ps.
    initial_time : float, optional 
        Time of first auxilairy step (in ps). If not specified, will attempt to
        be determined from auxiliary data; otherwise defaults to 0ps.
    time_col : int, optional
        Index of column in auxiliary data storing time (default value ``None``).
    data_cols : list of str, optional
        Indices of columns containing data of interest to be stored in 
        `step_data` (defaults to all columns).
    constant_dt : bool, optional
        If true, will use dt/initial_time to calculate time even when time
        stored in auxiliary data (default value is ``True``).



    Attributes
    ----------
    step : int
        Number of the current auxiliary step, starting at 0.
    n_steps : int
        Total number of auxiliary steps
    n_cols : int
        Number of columns of data for each auxiliary step.
    time : float
        Time of current auxiliary step, as read from data (if present) or 
        calculated using `dt` and `initial_time`.
    times : list of float
        List of the times of each auxiliary step.
    step_data : ndarray
        Value(s) from the auxiliary data column(s) of interest (per `data_cols`) 
        for the current step.
      
    ts_data : ndarray
        List of 'step_data' from each auxiliary step assigned to the last-read
        trajectory timestep.
    ts_diff : list
        List of difference in time between the last-read trajectory timestep
        and each auxiliary step assigned to it.
    ts_rep : list of float
        Represenatative value of auxiliary data for current trajectory timestep.

    """

    # update when add new options
    represent_options = ['closest', 'average']
      
    def __init__(self, represent_ts_as='closest', name=None, cutoff=-1, 
                 dt=1, initial_time=0, time_col=None, data_cols=None, 
                 constant_dt=True):
        # allow name to be optional for when using reader separate from
        # trajectory.
        self.name = name

        self.represent_ts_as = represent_ts_as

        self.cutoff = cutoff

        # set initially to avoid error on first read
        self.n_cols = None
        self._data_cols = []
        self._time_col = None

        self.ts_data = None
        self.ts_rep = None

        self._initial_time = initial_time
        self._dt = dt
        self.constant_dt = constant_dt

        self.step = -1
        self._read_next_step()
        self.n_cols = len(self._data)

        self.time_col = time_col
        if data_cols:
            self.data_cols = data_cols
        else:
            self.data_cols = [i for i in range(self.n_cols) 
                              if i != self.time_col]

        # get dt and initial time from auxiliary data (if time included)
        if self.time_col is not None and self.constant_dt:
            self._initial_time = self.time
            self._read_next_step()
            self._dt = self.time - self._initial_time
            self.go_to_first_step()

    def __len__(self):
        """ Number of steps in auxiliary data. """
        return self.n_steps

    def next(self):
        """ Move to next step of auxiliary data. """
        self._read_next_step()
        return {'time': self.time, 'data': self.step_data}


    def __next__(self):
        """ Move to next step of auxiliary data. """
        return self.next()

    def __iter__(self):
        """ Iterate over all auxiliary steps. """
        self._restart()
        return self

    def _restart(self):
        """ Reset back to start; calling next() should read first step. """
        # Overwrite as appropriate
        self.step = -1
                
    def go_to_first_step(self):
        """ Return to and read first step. """
        # Overwrite as appropriate
        #  could also use go_to_step(0)
        self._restart()
        self._read_next_step()
                
    def _read_next_step(self):
        """ Move to next step and update data. """
        # Define in each auxiliary reader
        raise NotImplementedError(
            "BUG: Override _read_next_timestep() in auxiliary reader!")

    def read_ts(self, ts):
        """ Read auxiliary data corresponding to the trajectory timestep *ts*.

        Read the auxiliary steps 'assigned' to *ts* (the steps that are within
        *ts.dt*/2 of of the trajectory timestep/frame - ie. closer to *ts*
        than either the preceeding or following frame). Then calculate a 
        'representative value' for the timestep from the data in each of these 
        auxiliary steps, and add to *ts*.

        Paramaters
        ----------
        ts : :class:`~MDAnalysis.coordinates.base.Timestep` object
            The trajectory timestep for which corresponding auxiliary data is
            to be read.

        Returns
        -------
        The :class:`~MDAnalysis.coordinates.base.Timestep` object with the 
        representative value in *ts.aux* updated appropriately.

        Notes
        -----
        The auxiliary reader will end up positioned at the last step assigned
        to the trajectory frame or, if the frame includes no auxiliary steps,
        (as when auxiliary data are less frequent), the most recent auxiliary 
        step before the frame.

        """
        # Make sure our auxiliary step starts at the right point (just before
        # the frame being read): the current step should be assigned to a 
        # previous frame, and the next step to either the frame being read of a 
        # following frame. Move to right position if not.
        if not (self.step_to_frame(self.step, ts) < ts.frame
                and self.step_to_frame(self.step+1, ts) >= ts.frame):
            self.move_to_ts(ts)

        self._reset_ts() # clear previous ts data
        while self.step_to_frame(self.step+1, ts) == ts.frame:
            self._read_next_step()
            self._add_step_to_ts(ts.time)
        self.ts_rep = self.calc_representative()

        # check we have a name to set in ts
        if self.name:
            setattr(ts.aux, self.name, self.ts_rep)
        else:
            # assume we don't want to store in ts. Should only happen when
            # using AuxReader independantly of a trajectory.
            warnings.warn("Auxiliary name not set, assuming you don't want to "
                          "store representative value in *ts*.")
        return ts


    def step_to_frame(self, step, ts):
        """ Calculate closest trajectory frame for auxiliary step *step*.

        Calculated given dt and offset from *ts* as::
            frame = math.floor((self.times[step]-offset+ts.dt/2.)/ts.dt)

        Paramaters
        ----------
        step : int
            Step number to calculate closest trajectory frame for.
        ts : :class:`~MDAnalysis.coordinates.base.Timestep` object
            Timestep from the trajectory the calculated frame number is to 
            correspond to.

        Returns
        -------
        int or None
            Number of the trajectory frame closest (in time) to the given
            auxiliary step. If the step index is out of range for the auxiliary
            data, None is returned instead.

        Notes
        -----
        Note that the returned frame number may be out of range for the 
        trajectory.
        """
        if step not in range(self.n_steps):
            return None 
        offset = ts.data.get('time_offset', 0)
        return math.floor((self.times[step]-offset+ts.dt/2.)/ts.dt)

    def move_to_ts(self, ts):
        """ Position auxiliary reader just before trajectory timestep *ts*.

        Calling ``next()`` should read the first auxiliary step assigned to
        (closest to) the trajectory timestep *ts* or, if no auxiliary steps are 
        assigned to that timestep (as in the case of less frequent auxiliary
        data), the first auxiliary step after *ts*.

        Paramaters
        ----------
        ts : :class:`~MDAnalysis.coordinates.base.Timestep` object
            The trajectory timestep before which the auxiliary reader is to
            be positioned.
        """
        # Need to define in each auxiliary reader
        raise NotImplementedError(
            "BUG: Override move_to_ts() in auxiliary reader!")

    def _reset_ts(self):
        """ Clear existing timestep data. """
        self.ts_data = np.array([])
        self.ts_diffs = []

    def _add_step_to_ts(self, ts_time):
        """ Add data from the current step to *ts_data* and difference in time
        to *ts_time* (the trajectory timestep) of the current step to *ts_diffs*
        """
        if len(self.ts_data) == 0:
            self.ts_data = np.array([self.step_data])
        else:
            self.ts_data = np.append(self.ts_data, [self.step_data], axis=0)
        self.ts_diffs.append(abs(self.time - ts_time))

    def calc_representative(self):
        """ Calculate represenatative auxiliary value(s) from *ts_data*.
        
        Currently available options for calculating represenatative value are:
          *'closest': default; the value(s) from the step closest to in time to 
           the trajectory timestep
          *'average': average of the value(s) from steps 'assigned' to the 
           trajectory timestep.
        Additionally, if *cutoff* is specified, only steps within this time 
        of the trajectory timestep are considered in calculating the 
        represenatative.

        Returns
        -------
        List (of length *n_cols*) of auxiliary value(s) 'representing' the 
        timestep.
        """
        if self.cutoff != -1:
            cutoff_data = np.array([self.ts_data[i] 
                                    for i,d in enumerate(self.ts_diffs)
                                    if d < self.cutoff])
            cutoff_diffs = [d for d in self.ts_diffs if d < self.cutoff]
        else:
            cutoff_data = self.ts_data
            cutoff_diffs = self.ts_diffs
        if len(cutoff_data) == 0:
            value = np.array([np.nan]*len(self.data_cols))
        elif self.represent_ts_as == 'closest':
            value = cutoff_data[np.argmin(cutoff_diffs),:]
        elif self.represent_ts_as == 'average':
            value = np.mean(cutoff_data, axis=0)
        return value
    
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    def close(self):
        # Overwrite as appropriate
        pass

    @property
    def time(self):
        """ Time in ps of current auxiliary step.

        As read from the appropriate column of the auxiliary data, if present; 
        otherwise calcuated as ``step * dt + initial_time``
        """
        if self.time_col is not None:
            return self._data[self.time_col]
        else:
            return self.step * self.dt + self.initial_time

    @property
    def step_data(self):
        """ Auxiliary values of interest for the current step.

        As taken from the appropariate columns (identified in `data_cols`) of 
        the full auxiliary data read in for the current step.
        """
        return [self._data[i] for i in self.data_cols]

    @property
    def n_steps(self):
        """ Total number of steps in the auxiliary data. """
        try:
            return self._n_steps
        except AttributeError:
            self._n_steps = self.count_n_steps()
            return self._n_steps

    @property
    def times(self):
        """ List of times of each step in the auxiliary data. 

        Calculated using `dt` and `initial_time` if `constant_dt` is True; 
        otherwise as read from each auxiliary step in turn. 
        """
        try:
            return self._times
        except AttributeError:
            if self.constant_dt:
                self._times = [i*self.dt+self.initial_time 
                               for i in range(self.n_steps)]
            else:
                self._times = self.read_all_times()
            return self._times

    @property
    def dt(self):
        """ Change in time between steps. 

        Defaults to 1ps if not provided or read from auxiliary data. """
        return self._dt

    @property
    def initial_time(self):
        """ Time corresponding to first auxiliary step. 

        Defaults to 0ps if not provided or read from auxilairy data. """
        return self._initial_time

    @property
    def time_col(self):
        """ Get or set index of column in auxiliary data containing step time.
        """
        return self._time_col

    @time_col.setter
    def time_col(self, new): 
        try:
            self._data[new]
        except IndexError:
            raise ValueError("{0} is not a valid index for time column (num. "
                             "cols is {1})".format(new, self.n_cols))
        self._time_col = new

    @property
    def data_cols(self):
        """ Get or set list of column indices in auxiliary data containing 
        value(s) of interest. These are the values that will be stored in 
        *step_data*/*ts_data*/*ts_rep*.
        """ 
        return self._data_cols

    @data_cols.setter
    def data_cols(self, new):
        for col in new:
            try:
                self._data[col]
            except IndexError:
                raise ValueError("{0} not a valid data column index ("
                                 "num. cols is {1})".format(col, self.n_cols))
        self._data_cols = new

    @property
    def represent_ts_as(self):
        """ Get or set method by which 'representative' timestep values of 
        auxilairy data will be calculated.
        """
        return self._represent_ts_as

    @represent_ts_as.setter
    def represent_ts_as(self, new):
        if new not in self.represent_options:
            raise ValueError("{0} is not a valid option for calculating "
                             "representative value(s). Enabled options are: "
                             "{1}".format(new, self.represent_options))
        self._represent_ts_as = new


    def __del__(self):
        self.close()

    
class AuxFileReader(AuxReader):
    """ Base class for auxiliary readers that read from file.

    Extends AuxReader with attributes and methods particular to reading 
    auxiliary data from an open file, for use when auxiliary files may be too
    large to read in at once.

    Paramaters
    ----------
    filename : str
       Location of the file containing the auxiliary data.
    **kwargs
       Other AuxReader options.    

    Attributes
    ----------
    auxfile
       File object for the auxiliary file.

    See also
    --------
    :class:`AuxReader`
    """
    
    def __init__(self, filename, **kwargs):
        self.filename = filename
        self.auxfile = open(filename)
        
        super(AuxFileReader, self).__init__(**kwargs)

    def close(self):
        """ Close *auxfile*. """
        if self.auxfile == None:
            return
        self.auxfile.close()
        self.auxfile = None

    def _restart(self):
        """ Reposition to just before first step. """
        self.auxfile.seek(0)
        self.step = -1
        
    def _reopen(self):
        """ Close and then reopen *auxfile*. """
        self.auxfile.close()
        self.auxfile = open(self.filename)
        self.step = -1

    def move_to_ts(self, ts):
        """ Position auxiliary reader just before trajectory timestep *ts*.

        Calling ``next()`` should read the first auxiliary step assigned to
        (closest to) the trajectory timestep *ts* or, if no auxiliary steps are 
        assigned to that timestep (as in the case of less frequent auxiliary
        data), the first auxiliary step after *ts*.

        Paramaters
        ----------
        ts : :class:`~MDAnalysis.coordinates.base.Timestep` object
            The trajectory timestep before which the auxiliary reader is to
            be positioned.

        Notes
        -----
        Works by reading through all timesteps consecutively until correct 
        timestep is reached. Overwrite if this can be done more efficiently.
        """
        # only restart if we're currently beyond *ts*
        if self.step_to_frame(self.step, ts) >= ts.frame:
            self._restart()

        # read through each step till we reach the right place
        while self.step_to_frame(self.step+1, ts) < ts.frame:
            # avoid restarting when we _read_next past final step
            if self.step == self.n_steps-1:
                return
            self._read_next_step()

 

    def go_to_step(self, i):
        """ Move to and read i-th auxiliary step. 

        Paramaters
        ----------
        i : int
            Step number (0-indexed) to move to

        Raises
        ------
        ValueError
            If step index not in valid range.

        Notes
        -----
        Works by reading through all steps consecutively until correct step
        is reached. Overwrite if this can be done more efficiently.
        """
        ## could seek instead?
        if i not in range(self.n_steps):
            raise ValueError("Step index {0} is not valid for auxiliary "
                             "(num. steps {1}!".format(i, self.n_steps))
        value = self.go_to_first_step()
        while self.step != i:
            value = self.next()
        return value
