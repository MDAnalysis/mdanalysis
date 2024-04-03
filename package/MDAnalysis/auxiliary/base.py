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
Auxiliary Readers --- :mod:`MDAnalysis.auxiliary.base`
======================================================

Base classes for deriving all auxiliary data readers. See the API in :mod:`MDAnalysis.auxiliary.__init__`.

.. autoclass:: AuxStep
   :members:

.. autoclass:: AuxReader
   :members:

.. autoclass:: AuxFileReader
   :members:

"""

import os
import numbers
import math
import warnings
from typing import Union, Optional, Dict
from collections import defaultdict
import pickle

import numpy as np

from ..lib.util import asiterable, anyopen
from .. import units

from . import _AUXREADERS
from .core import auxreader


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
    """Base class for auxiliary timesteps.

    Stores the auxiliary data for the current auxiliary step. On creation,
    ``step`` is set to -1.

    .. versionchanged:: 2.4.0
        Added memory_limit parameter to control raising memory usage warnings.

    Parameters
    ----------
    dt : float, optional
        Change in time between auxiliary steps (in ps). If not specified, will
        attempt to determine from auxiliary data; otherwise defaults to 1 ps.
        Ignored if ``constant_dt`` is False.
    initial_time : float, optional
        Time of first auxiliary step (in ps). If not specified, will attempt to
        determine from auxiliary data; otherwise defaults to 0 ps. Ignored if
        ``constant_dt`` is False.
    time_selector: optional
        Key to select 'time' value from the full set of data read for each
        step, if time selection is enabled; type will vary depending on the
        auxiliary data format (see individual AuxReader documentation). If
        ``None`` (default value), time is instead calculated as: ``time = step
        * dt + initial_time``
    data_selector: optional
        Key(s) to select auxiliary data values of interest from the full set of
        data read for each step, if data selection is enabled by the reader;
        type will vary depending on the auxiliary data format (see individual
        AuxReader documentation).
        If ``None`` (default value), the full set of data is returned.
    constant_dt : bool, optional
        (Default: True) Set to False if ``dt`` is not constant
        throughout the auxiliary data set, in which case a valid
        ``time_selector`` must be provided.
    memory_limit : float, optional
        Sets the threshold of memory usage by auxiliary data  (in GB) at which
        to issue a warning. Default: 1 GB.

    Attributes
    ----------
    step : int
        Number of the current auxiliary step (0-based).
    """

    def __init__(self, dt=1, initial_time=0, time_selector=None,
                 data_selector=None, constant_dt=True, memory_limit=None):
        self.step = -1
        self._initial_time = initial_time
        self._dt = dt
        self._constant_dt = constant_dt
        # check for valid values when assigning _time/data_selector will fail
        # here as we don't have and _data yet, so set _time/data_selector_ directly;
        # if invalid, will catch later
        self._time_selector_ = time_selector
        self._data_selector_ = data_selector

    @property
    def time(self):
        """ Time in ps of current auxiliary step (as float).

        Read from the set of auxiliary data read each step if time selection
        is enabled and a valid ``time_selector`` is specified; otherwise
        calculated as ``step * dt + initial_time``.
        """
        if self._time_selector is not None:
            return self._select_time(self._time_selector)
        elif self._constant_dt:
            # default to calculting time...
            return self.step * self._dt + self._initial_time
        else:
            raise ValueError("If dt is not constant, must have a valid "
                             "time selector")


    @property
    def data(self):
        """ Auxiliary values of interest for the current step (as ndarray).

        Read from the full set of data read for each step if data selection is
        enabled and a valid ``data_selector`` is specified; otherwise
        defaults to the full set of data.
        """
        if self._data_selector is not None:
            return self._select_data(self._data_selector)
        # default to full set of data...
        return self._data

    @property
    def _time_selector(self):
        """ 'Key' to select time from the full set of data read in each step.

        Will be passed to ``_select_time()``, defined separately for each
        auxiliary format, when returning the time of the current step.
        Format will depend on the auxiliary format. e.g. for the XVGReader,
        this is an index and ``_select_time()`` returns the value in that column
        of the current step data.

        Defaults to 'None' if time selection is not enabled.
        """
        try:
            self._select_time
        except AttributeError:
            warnings.warn("{} does not support time selection. Reverting to "
                          "default".format(self.__class__.__name__))
            return None
        return self._time_selector_

    @_time_selector.setter
    def _time_selector(self, new):
        # check we have a select_time method
        try:
            select = self._select_time
        except AttributeError:
            warnings.warn("{} does not support time selection".format(
                                                       self.__class__.__name__))
        else:
            # check *new* is valid before setting; _select_time should raise
            # an error if not
            select(new)
            self._time_selector_ = new

    @property
    def _data_selector(self):
        """ 'Key' to select values of interest from full set of auxiliary data.
        These are the values that will be stored in ``data`` (and
        ``frame_data`` and ``frame_rep``).

        Will be passed to ``_select_data()``, defined separately for each
        auxiliary format, when returning the data of interest for the current
        step (``data``). Format will depend on the auxiliary format; e.g.
        for the XVGReader, this is a list of indices and `_select_data()` returns
        the value(s) in those columns of the current step data.

        Defaults to 'None' if data selection is not enabled.
        """
        try:
            self._select_data
        except AttributeError:
            warnings.warn("{} does not support data selection. Reverting to "
                          "default".format(self.__class__.__name__))
            return None
        return self._data_selector_

    @_data_selector.setter
    def _data_selector(self, new):
        # check we have a select_data method
        try:
            select = self._select_data
        except AttributeError:
            warnings.warn(
                "{} does not support data selection".format(self.__class__.__name__)
            )
        else:
            # check *new* is valid before setting; _select_data should raise an
            # error if not
            select(new)
            self._data_selector_ = new

    def _empty_data(self):
        """ Create an 'empty' ``data``-like placeholder.

        Returns an ndarray in the format of ``data`` with all values set to
        np.nan; to use at the 'representative value' when no auxiliary steps
        are assigned to a trajectory timestep/within the cutoff.

        Default behaviour here works when ``data`` is a ndarray of floats. May
        need to overwrite in individual format's AuxSteps.

        .. versionchanged:: 2.4.0
            dtype changed to np.float64 from np.float_
        """
        return np.full_like(self.data, np.nan, dtype=np.float64)


class AuxReader(metaclass=_AuxReaderMeta):
    """Base class for auxiliary readers.

    Allows iteration over a set of data from a trajectory, additional
    ('auxiliary') to the regular positions/velocities/etc. This auxiliary
    data may be stored in e.g. an array or a separate file.

    See the :ref:`Auxiliary API` for more on use.

    .. versionchanged:: 2.4.0
        Behaviour of ``cutoff`` changed, default parameter which specifies
        not cutoff is set is now None, not -1.

    .. versionchanged:: 2.4.0
        :class:`AuxReader` instances now have a :meth:`copy` method which
        creates a deep copy of the instance.

    Parameters
    ----------
    auxname : str, optional
        Name for auxiliary data. When added to a trajectory, the representative
        auxiliary value(s) for the timestep may be accessed as ``ts.aux.auxname``
        or ``ts.aux['auxname']``.
    represent_ts_as : str
        Method to use to calculate representative value of auxiliary data for a
        trajectory timestep. See :func:`calc_representative` for valid options.
    cutoff : float, optional
        Auxiliary steps further from the trajectory timestep than *cutoff*
        (in ps) will be ignored when calculating representative values. If None
        (default), all auxiliary steps assigned to that timestep will be used.
    **kwargs
        Options to be passed to :class:`~AuxStep`


    Attributes
    ----------
    auxstep :
       :class:`~AuxStep` object containing data for current step.
    frame_data : dict
        Dictionary containing ``data`` from each auxiliary step assigned to the
        current trajectory timestep, indexed by the difference in time between
        the step and trajectory timestep (i.e. ``auxstep.time - ts.time``; in ps)
    frame_rep : ndarray
        Representative value(s) of auxiliary data for current trajectory timestep.


    Note
    ----
    Auxiliary data are assumed to be time ordered and contain no duplicates.
    """

    _Auxstep = AuxStep
    # update when add new options
    represent_options = ['closest', 'average']

    # list of attributes required to recreate the auxiliary
    required_attrs = ['represent_ts_as', 'cutoff', 'dt', 'initial_time',
                      'time_selector', 'data_selector', 'constant_dt', 'auxname',
                      'format', '_auxdata']

    def __init__(self, represent_ts_as='closest', auxname=None, cutoff=None,
                 **kwargs):
        # allow auxname to be optional for when using reader separate from
        # trajectory.
        self.auxname = auxname
        self.represent_ts_as = represent_ts_as
        # no cutoff was previously passed as -1. This check is to maintain this
        # behaviour: negative number is appropriately turned to None
        self.cutoff = cutoff if cutoff is not None and cutoff >= 0 else None
        self.frame_data = None
        self.frame_rep = None

        self.auxstep = self._Auxstep(**kwargs)
        self._read_next_step()

        # if dt is constant and auxiliary data includes time, calculate
        # initial time and dt
        if self.time_selector is not None and self.constant_dt:
            self.auxstep._initial_time = self.time
            self._read_next_step()
            self.auxstep._dt = self.time - self.initial_time
            self.rewind()

    def copy(self):
        """Returns a deep copy of the AuxReader"""
        orig_dict = pickle.dumps(self)
        new_reader = pickle.loads(orig_dict)
        return new_reader

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
        """ Move to next step and update auxstep.

        Should return the AuxStep instance corresponding to the next step.
        """
        # Define in each auxiliary reader
        raise NotImplementedError(
            "BUG: Override _read_next_step() in auxiliary reader!")

    def update_ts(self, ts):
        """ Read auxiliary steps corresponding to and update the trajectory
        timestep *ts*.

        Calls :meth:`read_ts`, then updates *ts* with the representative value.
        ``auxname`` must be set; the representative value will be accessible in
        *ts* as ``ts.aux.auxname`` or ``ts.aux['auxname']``.

        Parameters
        ----------
        ts : :class:`~MDAnalysis.coordinates.timestep.Timestep` object
            The trajectory timestep for which corresponding auxiliary data is
            to be read and updated.

        Returns
        -------
        :class:`~MDAnalysis.coordinates.timestep.Timestep`
            *ts* with the representative auxiliary
            value in ``ts.aux`` be updated appropriately.

        Raises
        ------
        ValueError
            If ``auxname`` is not set.

        See Also
        --------
        :meth:`read_ts`
        """
        if not self.auxname:
            raise ValueError("Auxiliary name not set, cannot set representative "
                             "value in timestep. Name auxiliary or use read_ts "
                             "instead")
        self.read_ts(ts)
        setattr(ts.aux, self.auxname, self.frame_rep)
        return ts

    def read_ts(self, ts):
        """ Read auxiliary steps corresponding to the trajectory timestep *ts*.

        Read the auxiliary steps 'assigned' to *ts* (the steps that are within
        ``ts.dt/2`` of of the trajectory timestep/frame - ie. closer to *ts*
        than either the preceding or following frame). Then calculate a
        'representative value' for the timestep from the data in each of these
        auxiliary steps.

        To update *ts* with the representative value, use ``update_ts`` instead.

        Parameters
        ----------
        ts : :class:`~MDAnalysis.coordinates.timestep.Timestep` object
            The trajectory timestep for which corresponding auxiliary data is
            to be read.

        See Also
        --------
        :meth:`update_ts`

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
        frame_for_step = self.step_to_frame(self.step, ts)
        frame_for_next_step = self.step_to_frame(self.step+1, ts)
        if frame_for_step is not None:
            if frame_for_next_step is None:
                # self.step is the last auxiliary step in memory.
                if frame_for_step >= ts.frame:
                    self.move_to_ts(ts)
            elif not (frame_for_step < ts.frame <= frame_for_next_step):
                self.move_to_ts(ts)

        self._reset_frame_data() # clear previous frame data
        # read through all the steps 'assigned' to ts.frame + add to frame_data
        while self.step_to_frame(self.step+1, ts) == ts.frame:
            self._read_next_step()
            self._add_step_to_frame_data(ts.time)
        self.frame_rep = self.calc_representative()

    def attach_auxiliary(self,
                         coord_parent,
                         aux_spec: Optional[Union[str, Dict[str, str]]] = None,
                         format: Optional[str] = None,
                         **kwargs) -> None:
        """Attaches the data specified in `aux_spec` to the `coord_parent`

        This method is called from within
        :meth:`MDAnalysis.coordinates.base.ReaderBase.add_auxiliary()`.
        `add_auxiliary` should be agnostic to the type of AuxReader, so the
        method call leads here instead. First, some checks are done on
        the arguments to make sure the input is treated properly. Then,
        the AuxReader(s) with appropriate :attr:`data_selector` are
        associated with the `coord_parent` from which `add_auxiliary` was
        called.

        Parameters
        ----------
        coord_parent : MDAnalysis.coordinates.base.ReaderBase
            Reader object to which to attach the auxiliary data.

        aux_spec : str, Dict[str, str], None
            Specifies which data to add to `coord_parent`. String types are
            for :class:`MDAnalysis.auxiliary.XVG.XVGReader` only. Dictionaries
            are the standard way of providing `aux_spec` information (see also:
            :mod:`MDAnalysis.auxiliary.EDR`).
            Passing `None` causes all data to be added.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If trying to add data under an `aux_spec` key that is already
            assigned.

        """
        if "auxname" in kwargs:
            # This check is necessary for the tests in coordinates/base
            # `test_reload_auxiliaries_from_description`
            # because there, auxreaders are created from descriptions without
            # giving explicit `aux_spec`s
            aux_spec = {kwargs["auxname"]: None}

        elif aux_spec is None:
            # Add all terms found in the file if aux_spec is None
            aux_spec = {term: term for term in self.terms}

        elif isinstance(aux_spec, str):
            # This is to keep XVGReader functioning as-is, working with strings
            # for what used to be `auxname`. Often, no `auxterms` are specified
            # for XVG files. The below sets the data selector to None,
            # the default for XVGReader
            aux_spec = {aux_spec: None}

        for auxname in aux_spec:
            if auxname in coord_parent.aux_list:
                raise ValueError(f"Auxiliary data with name {auxname} already "
                                 "exists")
            if " " in auxname:
                warnings.warn(f"Auxiliary name '{auxname}' contains a space. "
                              "Only dictionary style accession, not attribute "
                              f"style accession of '{auxname}' will work.")
            description = self.get_description()
            # make sure kwargs are also used
            description_kwargs = {**description, **kwargs}
            # Make a copy of the auxreader for every attribute to add
            # This is necessary so all attributes can be iterated over
            aux = auxreader(**description_kwargs)
            aux.auxname = auxname
            if aux.data_selector is None:
                # When calling ReaderBase.copy(), aux_spec information is lost
                # but data_selector is retained
                aux.data_selector = aux_spec[auxname]
            coord_parent._auxs[auxname] = aux
            coord_parent.ts = aux.update_ts(coord_parent.ts)

        aux_memory_usage = 0
        # Check if testing, which needs lower memory limit
        memory_limit = kwargs.get("memory_limit", 1e+09)
        for reader in coord_parent._auxs.values():
            aux_memory_usage += reader._memory_usage()
        if aux_memory_usage > memory_limit:
            conv = 1e+09  # convert to GB
            warnings.warn("AuxReader: memory usage warning! "
                          f"Auxiliary data takes up {aux_memory_usage/conv} "
                          f"GB of memory (Warning limit: {memory_limit/conv} "
                          "GB).")

    def _memory_usage(self):
        raise NotImplementedError("BUG: Override _memory_usage() "
                                  "in auxiliary reader!")

    def step_to_frame(self, step, ts, return_time_diff=False):
        """ Calculate closest trajectory frame for auxiliary step *step*.

        Calculated given dt, time and frame from *ts*::

            time_frame_0 = ts.time - ts.frame*ts.dt   # time at frame 0

            frame = floor((step_to_time(step) - time_frame_0 + ts.dt/2)/ts.dt))

        The difference in time between the step and the calculated frame can
        also optionally be returned with *return_time_diff*.

        Parameters
        ----------
        step : int
            Number of the auxiliary step to calculate closest trajectory frame
            for.
        ts : :class:`~MDAnalysis.coordinates.timestep.Timestep` object
            (Any) timestep from the trajectory the calculated frame number is to
            correspond to.
        return_time_diff : bool, optional
            (Default: False) Additionally return the time difference between
            *step* and returned frame.

        Returns
        -------
        frame_index : int or None
            Number of the trajectory frame closest (in time) to the given
            auxiliary step. If the step index is out of range for the auxiliary
            data, ``None`` is returned instead.
        time_diff : float (optional)
            Difference in time between *step* and *frame_index*.

        Note
        ----
        Assumes trajectory dt is consant.
        The returned frame number may be out of range for the trajectory.
        """
        if step >= self.n_steps or step < 0:
            return None
        time_frame_0 = ts.time - ts.frame*ts.dt  # assumes ts.dt is constant
        time_step = self.step_to_time(step)
        frame_index = int(math.floor((time_step-time_frame_0+ts.dt/2.)/ts.dt))
        if not return_time_diff:
            return frame_index
        else:
            time_frame = time_frame_0 + frame_index*ts.dt
            time_diff = abs(time_frame - time_step)
            return frame_index, time_diff

    def move_to_ts(self, ts):
        """ Position auxiliary reader just before trajectory timestep *ts*.

        Calling ``next()`` should read the first auxiliary step 'assigned' to
        the trajectory timestep *ts* or, if no auxiliary steps are
        assigned to that timestep (as in the case of less frequent auxiliary
        data), the first auxiliary step after *ts*.

        Parameters
        ----------
        ts : :class:`~MDAnalysis.coordinates.timestep.Timestep` object
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

    def next_nonempty_frame(self, ts):
        """ Find the next trajectory frame for which a representative auxiliary
        value can be calculated.

        That is, the next trajectory frame to which one or more auxiliary steps
        are assigned and fall within the cutoff.

        Starts looking from the current step time. If the end of the auxiliary
        data is reached before a trajectory frame is found, None is returned.

        Parameters
        ----------
        ts : :class:`~MDAnalysis.coordinates.timestep.Timestep` object
            Any timestep from the trajectory for which the next 'non-empty'
            frame is to be found.

        Returns
        -------
        int
            Index of the next auxiliary-containing frame in the trajectory.

        Note
        ----
        The returned index may be out of range for the trajectory.
        """
        step = self.step
        while step < self.n_steps-1:
            next_frame, time_diff = self.step_to_frame(self.step+1, ts,
                                                       return_time_diff=True)
            if self.cutoff is not None and time_diff > self.cutoff:
                # 'representative values' will be NaN; check next step
                step = step + 1
            else:
                return next_frame
        # we ran out of auxiliary steps...
        return None

    def __getitem__(self, i):
        """ Return the AuxStep corresponding to the *i*-th auxiliary step(s)
        (0-based). Negative numbers are counted from the end.

        *i* may be an integer (in which case the corresponding AuxStep is
        returned) or a list of integers or slice (in which case an iterator is
        returned)::

             step_10 = aux_reader[10]

        will move to step 10 of the auxiliary and return the :class:`AuxStep`.
        By using a slice/list, we can iterated over specified steps in the
        auxiliary, e.g. when performing analysis ::

             for auxstep in aux_reader[100:200]:   # analyse only steps 100 to 200
                 run_analysis(auxstep)

             for auxstep in aux_reader[::10]       # analyse every 10th step
                 run_analysis(auxstep)
        """
        if isinstance(i, numbers.Integral):
            i = self._check_index(i)
            return self._go_to_step(i)

        elif isinstance(i, (list, np.ndarray)):
            return self._list_iter([self._check_index(x) for x in i])

        elif isinstance(i, slice):
            # default start to first frame (ie. 0)
            start = self._check_index(i.start) if i.start is not None else 0
            # default stop to after last frame (i.e. n_steps)
            # n_steps is a valid stop index but will fail _check_index;
            # deal with separately
            stop = (i.stop if i.stop == self.n_steps
                    else self._check_index(i.stop) if i.stop is not None
                    else self.n_steps)
            step = i.step or 1
            if not isinstance(step, numbers.Integral) or step < 1:
                raise ValueError("Step must be positive integer") # allow -ve?
            if start > stop:
                raise IndexError("Stop frame is lower than start frame")
            return self._slice_iter(slice(start,stop,step))
        else:
            raise TypeError("Index must be integer, list of integers or slice")

    def _check_index(self, i):
        if not isinstance(i, numbers.Integral):
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
        """ Calculate representative auxiliary value(s) from the data in
        *frame_data*.

        Currently implemented options for calculating representative value are:

          * `closest`: default; the value(s) from the step closest to in time
            to the trajectory timestep

          * `average`: average of the value(s) from steps 'assigned' to the
            trajectory timestep.

        Additionally, if ``cutoff`` is specified, only steps within this time
        of the trajectory timestep are considered in calculating the
        representative.

        If no auxiliary steps were assigned to the timestep, or none fall
        within the cutoff, representative values are set to ``np.nan``.

        Returns
        -------
        ndarray
            Array of auxiliary value(s) 'representative' for the timestep.
        """
        if self.cutoff is None:
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
            # we don't know the original sign, and might have two equally-spaced
            # steps; check the earlier time first
            try:
                value = cutoff_data[-min_diff]
            except KeyError:
                value = cutoff_data[min_diff]
        elif self.represent_ts_as == 'average':
            try:
                value = np.mean(np.array(
                                [val for val in cutoff_data.values()]
                                ), axis=0)
            except TypeError:
                # for readers like EDRReader, the above does not work
                # because each step contains a dictionary of numpy arrays
                # as data
                value = defaultdict(float)
                for dataset in cutoff_data:
                    for term in self.terms:
                        value[term] += cutoff_data[dataset][term]
                for term in value:
                    value[term] = value[term] / len(cutoff_data)
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
        time : float
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
        """ Method by which 'representative' timestep values of auxiliary data
        will be calculated.
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
        """ Get the values of the parameters necessary for replicating the
        AuxReader.

        An AuxReader can be duplicated using
        :func:`~MDAnalysis.auxiliary.core.auxreader`::

            description = original_aux.get_description()
            new_aux = MDAnalysis.auxiliary.auxreader(**description)

        The resulting dictionary may also be passed directly to
        :meth:`~MDAnalysis.coordinates.base.ProtoReader.add_auxiliary` to
        reload an auxiliary into a trajectory::

             trajectory.add_auxiliary(**description)

        Returns
        -------
        dict
            Key-word arguments and values that can be used to replicate the
            AuxReader.
        """
        description = {attr.strip('_'): getattr(self, attr)
                       for attr in self.required_attrs}
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
        """Time of current auxiliary step (as stored in ``auxstep``; in ps)"""
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
        """Key to select 'time' value from the full set of data read for each step.
        As stored in ``austep``.

        Type differs between auxiliary formats, depending how the data for each
        step is read in and stored; e.g. data from .xvg files is read in as a
        list and `time_selector` must be a valid index. If time selection is not
        enabled by the reader, ``time_selector`` will default to ``None``.

        See each individual auxiliary reader.
        """
        return self.auxstep._time_selector

    @time_selector.setter
    def time_selector(self, new):
        old = self.auxstep._time_selector
        self.auxstep._time_selector = new
        if old != new:
            # if constant_dt is False and so we're using a _times list, this will
            # now be made invalid
            try:
                del(self._times)
            except AttributeError:
                pass

    @property
    def data_selector(self):
        """Key(s) to select auxiliary data values of interest from the full set
        of data read for each step (as stored in ``auxstep``).

        Type differs between auxiliary formats, depending how the data for each
        step is read in and stored - e.g. data from .xvg files is read in as
        a list and `data_selector` must be a list of valid indicies. If data
        selection is not enabled by the reader, ``data_selector`` will default
        to ``None``.

        See each individual auxiliary reader.
        """
        return self.auxstep._data_selector

    @data_selector.setter
    def data_selector(self, new):
        self.auxstep._data_selector = new

    @property
    def constant_dt(self):
        """ True if ``dt`` is constant throughout the auxiliary (as stored in
        ``auxstep``) """
        return self.auxstep._constant_dt

    @constant_dt.setter
    def constant_dt(self, new):
        self.auxstep._constant_dt = new


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
        self.auxfile = anyopen(filename)
        self._auxdata = os.path.abspath(filename)
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
        if self.auxfile is not None:
            self.auxfile.close()
        self.auxfile = open(self._auxdata)
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
