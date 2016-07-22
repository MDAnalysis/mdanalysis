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

Base classes for deriving all auxiliary data readers. See the API in :mod:`MDAnalysis.auxiliary.__init__`.

.. autoclass:: AuxReader
   :members:

.. autoclass:: AuxFileReader
   :members:

"""

import six
from six.moves import range

import os
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


class AuxStep(object):
    """ Base class for auxiliary timesteps.

    Stores the auxiliary data for the current auxiliary step.

    Parameters
    ----------
    dt : float, optional
        Change in time between auxiliary steps (in ps). If not specified, will
        attempt to be determined from auxiliary data; otherwise defaults to 1 ps.
    initial_time : float, optional 
        Time of first auxiliary step (in ps). If not specified, will attempt to
        be determined from auxiliary data; otherwise defaults to 0 ps.
    time_selector: optional
        Key to select 'time' value from the full set of data read for each step,
        if time selection is enabled; type will vary depending on the auxiliary
        data format (see individual AuxReader documentation). 
        If ``None`` (default value), time is instead calculated from ``dt``, 
        ``initial_time`` and ``step``.
    data_selector: optional
        Key(s) to select auxilairy data values of interest from the full set of
        data read for each step, if data selection is enabled by the reader; 
        type will vary depending on the auxiliary data formar (See individual
        AuxReader documentation). 
        If ``None`` (default value), the full set of data is returned.

    Attributes
    ----------
    step : int
        Number of the current auxiliary step (0-based).
    time : float
        Time of current auxiliary step, as read from data (if ``time_selector``
        specified) or calculated using `dt` and `initial_time`.
    data : ndarray
        Value(s) of auxiliary data of interest for the current step.
    """ 
    def __init__(self, dt=1, initial_time=0, time_selector=None, 
                 data_selector=None):
        self.step = -1
        self._initial_time = initial_time
        self._dt = dt
        self._time_selector_ = time_selector
        self._data_selector_ = data_selector

    @property
    def time(self):
        """ Time in ps of current auxiliary step.

        As read from the full set of data if present (specified by 
        ``time_selector``); otherwise calcuated as ``step * dt + initial_time``.

        Note
        ----
        Assumed to be in ps.
        """
        if self._time_selector is not None:
            return self._select_time(self._time_selector)
        else:
            return self.step * self._dt + self._initial_time    
 
    @property
    def data(self):
        """ Auxiliary values of interest for the current step.

        As selected from the full set of data read it for each step by  
        ``data_selector``.
        """
        if self._data_selector is not None:
            return self._select_data(self._data_selector)
        else:
            return self._data

    @property
    def _time_selector(self):
        """ 'Key' to select time from the full set of data read in each step.

        Will be passed to ``_select_time()``, defined separately for each 
        auxiliary format, when returning the time of the current step. 
        Format will hence depend on the auxiliary format. e.g. for the XVGReader,
        this is an index and `_select_time()` returns the value in that column 
        of the current step data.
        """
        return self._time_selector_

    @_time_selector.setter
    def _time_selector(self, new):
        # check *new* is valid before setting; _select_time should raise 
        # an error if not
        self._select_time(new)
        self._time_selector_ = new

    @property
    def _data_selector(self):
        """ 'Key' to select values of interest from full set of auxiliary data. 
        These are the values that will be stored in ``data``, 
        ``frame_data``, and ``frame_rep``.

        Will be passed to ``_select_data()``, defined separately for each 
        auxiliary format, when returning the data of interest for the current
        step (``data``). Format will depend on the auxiliary format; e.g.
        for the XVGReader, this is a list of indices and `_select_data()` returns
        the value(s) in those columns of the current step data.
        """ 
        return self._data_selector_

    @_data_selector.setter
    def _data_selector(self, new):
        # check *new* is valid before setting; _select_data should raise an 
        # error if not
        self._select_data(new)
        self._data_selector_ = new


class AuxReader(six.with_metaclass(_AuxReaderMeta)):
    """ Base class for auxiliary readers.

    Allows iteration over a set of data from a trajectory, additional 
    ('auxiliary') to the regular positions/velocities/etc. This auxiliary 
    data may be stored in e.g. an array or a separate file.

    Auxiliary data may be added to a trajectory by 
    :meth:`MDAnalysis.coordinates.base.ProtoReader.add_auxiliary`, passing either an 
    AuxReader instance or the data/filename, in which case an appropriate reader 
    will be selected based on type/file extension.
   
    The AuxReader will handle alignment of the auxiliary and trajectory 
    timesteps (auxiliary data steps are 'assigned' to the closest trajectory 
    timestep), and for each trajectory timestep will provide a 'representative'
    auxiliary value (or values) based on the steps assigned to that timestep.


    Parameters
    ----------
    auxname : str, optional
        Name for auxiliary data. When added to a trajectory, the representative 
        auxiliary value(s) for the timestep are accessed as ``ts.aux.auxname`` 
        or ``ts.aux['auxname']``.
    represent_ts_as : {'closest', 'average'}
        How to calculate representative value of auxiliary data for a 
        trajectory timestep. Currently available:

          * 'closest': value from step closest to the trajectory timestep

          * 'average': average of values from auxiliary steps assigned to the trajectory timestep.

    cutoff : float, optional
        Auxiliary steps further from the trajectory timestep than *cutoff* 
        will be ignored when calculating representative values (the default
        value is -1, which indicates all auxiliary steps assigned to that 
        timestep will be used).
    constant_dt : bool, optional
        If true, will use ``dt``/``initial_time`` to calculate step time even 
        when time stored in auxiliary data (default value is ``True``).
    **kwargs
        Options to be passed to :class:`~AuxStep`


    Attributes
    ----------
    auxstep : AuxStep object
       :class:`~AuxStep` object containing data for current step. 
    n_steps : int
        Total number of auxiliary steps.
    frame_data : 
        'data' from each auxiliary step assigned to the last-read
        trajectory timestep.
    frame_rep : list of float
        Represenatative value of auxiliary data for current trajectory timestep.


    Note
    ----
    Auxiliary data are assumed to be time ordered and contain no duplicates.
    It step time is read from the data, it is assumed to be in ps units.
    """

    _Auxstep = AuxStep
    # update when add new options
    represent_options = ['closest', 'average']

    # list of attributes required to recreate the auxiliary
    required_attrs = ['represent_ts_as', 'cutoff', 'dt', 'initial_time',
                      'time_selector', 'data_selector', 'constant_dt', 'auxname', 
                      'format', '_data_input']
      
    def __init__(self, represent_ts_as='closest', auxname=None, cutoff=-1, 
                 constant_dt=True, **kwargs):
        # allow auxname to be optional for when using reader separate from
        # trajectory.
        self.auxname = auxname
        self.represent_ts_as = represent_ts_as
        self.cutoff = cutoff
        self.constant_dt = constant_dt
        self.frame_data = None
        self.frame_rep = None

        self.auxstep = self._Auxstep(**kwargs)
        self._read_next_step()

        # get dt and initial time from auxiliary data (if time included)
        if self.time_selector is not None and self.constant_dt:
            self.auxstep._initial_time = self.time
            self._read_next_step()
            self.auxstep._dt = self.time - self.initial_time
            self.rewind()

    def __len__(self):
        """ Number of steps in auxiliary data. """
        return self.n_steps

    def next(self):
        """ Move to next step of auxiliary data. """
        return self._read_next_step()

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
        self.auxstep.step = -1
                
    def rewind(self):
        """ Return to and read first step. """
        # Overwrite as appropriate
        #  could also use _go_to_step(0)
        self._restart()
        return self._read_next_step()
                
    def _read_next_step(self):
        """ Move to next step and update data. 

        Should return an AuxStep instance, with values for the appropriate step"""
        # Define in each auxiliary reader
        raise NotImplementedError(
            "BUG: Override _read_next_timestep() in auxiliary reader!")

    def read_ts(self, ts):
        """ Read auxiliary data corresponding to the trajectory timestep *ts*.

        Read the auxiliary steps 'assigned' to *ts* (the steps that are within
        ``ts.dt/2`` of of the trajectory timestep/frame - ie. closer to *ts*
        than either the preceeding or following frame). Then calculate a 
        'representative value' for the timestep from the data in each of these 
        auxiliary steps.

        If ``auxname`` is set, the representative value will also be added to the
        timestep, accessible as ``ts.aux.auxname`` or ``ts.aux['auxname']``.

        Parameters
        ----------
        ts : :class:`~MDAnalysis.coordinates.base.Timestep` object
            The trajectory timestep for which corresponding auxiliary data is
            to be read.

        Returns
        -------
        :class:`~MDAnalysis.coordinates.base.Timestep`
            If ``auxname`` is set, the representative auxiliary value in ``ts.aux`` 
            will be updated appropriately.

        Note
        ----
        The auxiliary reader will end up positioned at the last step assigned
        to the trajectory frame or, if the frame includes no auxiliary steps,
        (as when auxiliary data are less frequent), the most recent auxiliary 
        step before the frame.
        """
        # Make sure our auxiliary step starts at the right point (just before
        # the frame being read): the current step should be assigned to a 
        # previous frame, and the next step to either the frame being read or a 
        # following frame. Move to right position if not.
        if not (self.step_to_frame(self.step, ts) < ts.frame
                and self.step_to_frame(self.step+1, ts) >= ts.frame):
            self.move_to_ts(ts)

        self._reset_frame_data() # clear previous frame data
        while self.step_to_frame(self.step+1, ts) == ts.frame:
            self._read_next_step()
            self._add_step_to_frame_data(ts.time)
        self.frame_rep = self.calc_representative()

        # check we have a name to set in ts
        if self.auxname:
            setattr(ts.aux, self.auxname, self.frame_rep)
        else:
            # assume we don't want to store in ts. Should only happen when
            # using AuxReader independantly of a trajectory.
            warnings.warn("Auxiliary auxname not set, assuming you don't want to "
                          "store representative value in *ts*.")
        return ts


    def step_to_frame(self, step, ts):
        """ Calculate closest trajectory frame for auxiliary step *step*.

        Calculated given dt and offset from *ts* as::
            frame = math.floor((self.step_to_time(step)-offset+ts.dt/2.)/ts.dt)

        Parameters
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
            data, ``None`` is returned instead.

        Notes
        -----
        Note that the returned frame number may be out of range for the 
        trajectory.
        """
        if step >= self.n_steps:
            return None 
        offset = ts.data.get('time_offset', 0)
        return int(math.floor((self.step_to_time(step)-offset+ts.dt/2.)/ts.dt))

    def move_to_ts(self, ts):
        """ Position auxiliary reader just before trajectory timestep *ts*.

        Calling ``next()`` should read the first auxiliary step assigned to
        (closest to) the trajectory timestep *ts* or, if no auxiliary steps are 
        assigned to that timestep (as in the case of less frequent auxiliary
        data), the first auxiliary step after *ts*.

        Parameters
        ----------
        ts : :class:`~MDAnalysis.coordinates.base.Timestep` object
            The trajectory timestep before which the auxiliary reader is to
            be positioned.
        """
        # figure out what step we want to end up at
        if self.constant_dt:
            # if dt constant, calculate from dt/offset/etc
            step = int(math.floor((ts.time-ts.dt/2-self.initial_time)/self.dt))
            # if we're out of range of the number of steps, reset back
            step = max(min(step, self.n_steps-1), -1)
        else:
            # otherwise, go through steps till we find the right one
            for i in range(self.n_steps+1):
                if self.step_to_frame(i) >= ts.frame:
                    break
            # we want the step before this
            step = i-1
        if step == -1:
            self._restart()
        else:
            self._go_to_step(step)

    def __getitem__(self, i):
        """ Return the AuxStep corresponding to the *i*-th auxiliary step(s)
        (0-based). Negative numbers are counted from the end.
        
        *i* may be an integer (in which case the corresponding AuxStep is 
        returned) or a list or integers or slice (in which case an iterator is 
        returned). So ::

             aux_reader[10]

        will load data from step 10 of the auxiliary and return the 
        :class:`AuxStep``. By using a slice/list, we can iterated over specified 
        parts of the trajectory ::

             for auxstep in aux_reader[100:200]:   # analyse only steps 100 to 200
                 run_analysis(auxstep)

             for auxstep in aux_reader[::10]       # analyse every 10th step        
                 run_analysis(auxstep)
        """
        if isinstance(i, int):
            i = self._check_index(i)
            return self._go_to_step(i)
        elif isinstance(i, (list, np.ndarray)):
            return self._list_iter([self._check_index(x) for x in i])
        elif isinstance(i, slice):
            start = i.start or 0
            stop = i.stop or self.n_steps-1
            start = self._check_index(start)
            stop = self._check_index(stop)
            step = i.step or 1
            if not isinstance(step, int) or step < 1:
                raise ValueError("Step must be positive integer") # allow -ve?
            if start > stop:
                raise IndexError("Stop frame is lower than start frame")
            return self._slice_iter(slice(start,stop,step))
        else:
            raise TypeError("Index must be integer, list of integers or slice")

    def _check_index(self, i):
        if not isinstance(i, (int, np.integer)):
                raise TypeError("Step indices must be integers")
        if i < 0:
            i = i + self.n_steps
        if i < 0 or i >= self.n_steps:
            raise IndexError("{} is out of range of auxiliary (num. steps "
                             "{})".format(i, self.n_steps))
        return i

    def _list_iter(self, i):
        for j in i:
            yield self._go_to_step(j)

    def _slice_iter(self, i):
        for j in range(i.start, i.stop, i.step):
            yield self._go_to_step(j)

    def _go_to_step(self, i):
        """ Move to and read i-th auxiliary step. """
        # Need to define in each auxiliary reader
        raise NotImplementedError(
            "BUG: Override _go_to_step() in auxiliary reader!")

    def _reset_frame_data(self):
        self.frame_data = {}

    def _add_step_to_frame_data(self, ts_time):
        """ Update ``frame_data`` with values for the current step.

        Parameters
        ----------
        ts_time : float
            the time of the timestep the current step is being 'added to'. Used
            to calculate difference in time between current step and timestep.
        """
        time_diff = self.time - ts_time
        self.frame_data[time_diff] = self.auxstep.data

    def calc_representative(self):
        """ Calculate represenatative auxiliary value(s) from *frame_data*.
        
        Currently available options for calculating represenatative value are:
          * 'closest': default; the value(s) from the step closest to in time to 
           the trajectory timestep
          * 'average': average of the value(s) from steps 'assigned' to the 
           trajectory timestep.
        Additionally, if ``cutoff`` is specified, only steps within this time 
        of the trajectory timestep are considered in calculating the 
        represenatative.

        If no auxiliary steps were assigned to the timestep, or none fall
        within the cutoff, representative values are set to ``np.nan``.

        Returns
        -------
        Numpy array of auxiliary value(s) 'representing' the timestep.
        """
        if self.cutoff == -1:
            cutoff_data = self.frame_data
        else:
            cutoff_data = {key: val for key, val in self.frame_data.items()
                           if abs(key) <= self.cutoff}
        if len(cutoff_data) == 0:
            # no steps are 'assigned' to this trajectory frame, so return
            # values of ``np.nan``
            value = self.auxstep._empty_data()
        elif self.represent_ts_as == 'closest':
            min_diff = min([abs(i) for i in cutoff_data])
            # temp fix to get right sign. Go with -ve over +ve.
            try:
                value = cutoff_data[-min_diff]
            except KeyError:
                value = cutoff_data[min_diff]
        elif self.represent_ts_as == 'average':
            value = np.mean(np.array([val for val in cutoff_data.values()]),axis=0)
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
    def n_steps(self):
        """ Total number of steps in the auxiliary data. """
        try:
            return self._n_steps
        except AttributeError:
            self._n_steps = self._count_n_steps()
            return self._n_steps

    def step_to_time(self, i):
        """ Return time of auxiliary step *i*.

        Calculated using ``dt`` and ``initial_time`` if ``constant_dt`` is True; 
        otherwise from the list of times as read from the auxiliary data for 
        each step. 

        Parameters
        ----------
        i : int
            Index (0-based) of step to return time for

        Returns
        -------
        float
            Time (in ps) of step *i*

        Raises
        ------
        ValueError
            When *i* not in valid range
        """
        if i >= self.n_steps:
            raise ValueError("{0} is not a valid step index (total number of "
                             "steps is {1})".format(i, self.n_steps))
        if self.constant_dt:
            return i*self.dt+self.initial_time
        else:
            try:
                return self._times[i]
            except AttributeError:
                self._times = self.read_all_times()
                return self._times[i]

    @property
    def represent_ts_as(self):
        """ Get or set method by which 'representative' timestep values of 
        auxiliary data will be calculated.
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

    def get_description(self):
        """Get the args/kwargs necessary for dupicating the AuxReader

        The resulting description may be passed to 'add_auxiliary' to 
        reload an auxiliary into a trajectory.

        Returns
        -------
        dict
            Of args/kwargs that can be used to recreate the AuxReader
        """
        description = {}
        for attr in self.required_attrs:
            # '_data_input' is passed to add_auxiliary as 'auxdata' so change
            # that here. Should probably fix so we don't have to do this...
            if attr == '_data_input':
                description['auxdata'] = getattr(self, attr)
            else:
                description[attr] = getattr(self, attr)
        return description

    def __eq__(self, other):
        for attr in self.required_attrs:
            if getattr(self, attr) != getattr(other, attr):
                return False
        return True

    @property
    def step(self):
        """Number of the current auxiliary step (as stored in ``auxstep``; 
        0-based)."""
        return self.auxstep.step

    @property
    def time(self):
        """Time of current auxiliary step (as stored in ``auxstep``)"""
        return self.auxstep.time

    @property
    def dt(self):
        """Change in time between auxiliary steps (as stored in ``auxstep``; 
        in ps)"""
        return self.auxstep._dt

    @property
    def initial_time(self):
        """Time of first auxiliary step (as stored in ``auxstep``; in ps)"""
        return self.auxstep._initial_time

    @property
    def time_selector(self):
        """Key to select 'time' value from the full set of data read for each step
        (as stored in ``austep``).

        Type differs between auxiliary formats - e.g. for .xvg files, 
        `time_selector` is a valid index. See each individual auxiliary reader.
        """ 
        return self.auxstep._time_selector

    @time_selector.setter
    def time_selector(self, new):
        self.auxstep._time_selector = new

    @property
    def data_selector(self):
        """Key(s) to select auxilairy data values of interest from the full set 
        of data read for each step (as stored in ``auxstep``).

        Type differs between auxiliary formats - e.g. for .xvg files, 
        `data_selector is a list of valid indices. See each individual 
        auxiliary reader."""
        return self.auxstep._data_selector

    @data_selector.setter
    def data_selector(self, new):
        self.auxstep._data_selector = new


class AuxFileReader(AuxReader):
    """ Base class for auxiliary readers that read from file.

    Extends AuxReader with attributes and methods particular to reading 
    auxiliary data from an open file, for use when auxiliary files may be too
    large to read in at once.

    Parameters
    ----------
    filename : str
       Location of the file containing the auxiliary data.
    **kwargs
       Other AuxReader options.    

    See also
    --------
    :class:`AuxReader`

    Attributes
    ----------
    auxfile
       File object for the auxiliary file.

    """
    
    def __init__(self, filename, **kwargs):
        self.filename = filename
        self.auxfile = open(filename)
        self._data_input = os.path.abspath(filename)
        super(AuxFileReader, self).__init__(**kwargs)

    def close(self):
        """ Close *auxfile*. """
        if self.auxfile is None:
            return
        self.auxfile.close()
        self.auxfile = None

    def _restart(self):
        """ Reposition to just before first step. """
        self.auxfile.seek(0)
        self.auxstep.step = -1
        
    def _reopen(self):
        """ Close and then reopen *auxfile*. """
        if self.auxfile != None:
            self.auxfile.close()
        self.auxfile = open(self.filename)
        self.auxstep.step = -1

    def _go_to_step(self, i):
        """ Move to and read i-th auxiliary step. 

        Parameters
        ----------
        i : int
            Step number (0-indexed) to move to

        Raises
        ------
        ValueError
            If step index not in valid range.

        Note
        ----
        Works by reading through all steps consecutively until correct step
        is reached. Overwrite if this can be done more efficiently.
        """
        ## could seek instead?
        if i >= self.n_steps:
            raise ValueError("Step index {0} is not valid for auxiliary "
                             "(num. steps {1}!".format(i, self.n_steps))
        value = self.rewind()
        while self.step != i:
            value = self.next()
        return value
