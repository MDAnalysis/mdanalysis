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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
XVG auxiliary reader --- :mod:`MDAnalysis.auxiliary.XVG`
========================================================

xvg files are produced by Gromacs during simulation or analysis, formatted
for plotting data with Grace.

Data is column-formatted; time/data selection is enabled by providing column
indices.

Note
----
By default, the time of each step is assumed to be stored in the first column,
in units of ps.


.. autoclass:: XVGStep
   :members:


XVG Readers
-----------
The default :class:`XVGReader` reads and stores the full contents of the .xvg
file on initialisation, while a second reader (:class:`XVGFileReader`) that
reads steps one at a time as required is also provided for when a lower memory
footprint is desired.

Note
----
Data is assumed to be time-ordered.

Multiple datasets, separated in the .xvg file by '&', are currently not
supported (the readers will stop at the first line starting '&').


.. autoclass:: XVGReader
   :members:

.. autoclass:: XVGFileReader
   :members:


.. autofunction:: uncomment

"""
from __future__ import absolute_import

from six.moves import range

import numbers
import os
import numpy as np
from . import base
from ..lib.util import anyopen

def uncomment(lines):
    """ Remove comments from lines in an .xvg file

    Parameters
    ----------
    lines : list of str
        Lines as directly read from .xvg file

    Yields
    ------
    str
        The next non-comment line, with any trailing comments removed
    """
    for line in lines:
        stripped_line = line.strip()
        # ignore blank lines
        if not stripped_line:
            continue
        # '@' must be at the beginning of a line to be a grace instruction
        if stripped_line[0] == '@':
            continue
        # '#' can be anywhere in the line, everything after is a comment
        comment_position = stripped_line.find('#')
        if comment_position > 0 and stripped_line[:comment_position]:
            yield stripped_line[:comment_position]
        elif comment_position < 0 and stripped_line:
            yield stripped_line
        # if comment_position == 0, then the line is empty

class XVGStep(base.AuxStep):
    """ AuxStep class for .xvg file format.

    Extends the base AuxStep class to allow selection of time and
    data-of-interest fields (by column index) from the full set of data read
    each step.

    Parameters
    ----------
    time_selector : int | None, optional
        Index of column in .xvg file storing time, assumed to be in ps. Default
        value is 0 (i.e. first column).
    data_selector : list of int | None, optional
        List of indices of columns in .xvg file containing data of interest to
        be stored in ``data``. Default value is ``None``.
    **kwargs
        Other AuxStep options.

    See Also
    --------
    :class:`~MDAnalysis.auxiliary.base.AuxStep`
    """
    def __init__(self, time_selector=0, data_selector=None, **kwargs):
        super(XVGStep, self).__init__(time_selector=time_selector,
                                      data_selector=data_selector,
                                      **kwargs)

    def _select_time(self, key):
        if key is None:
            # here so that None is a valid value; just return
            return
        if isinstance(key, numbers.Integral):
            return self._select_data(key)
        else:
             raise ValueError('Time selector must be single index')

    def _select_data(self, key):
        if key is None:
            # here so that None is a valid value; just return
            return
        if isinstance(key, numbers.Integral):
            try:
                return self._data[key]
            except IndexError:
                raise ValueError('{} not a valid index for data with {} '
                                 'columns'.format(key, len(self._data)))
        else:
            return np.array([self._select_data(i) for i in key])


class XVGReader(base.AuxReader):
    """ Auxiliary reader to read data from an .xvg file.

    Detault reader for .xvg files. All data from the file will be read and stored
    on initialisation.

    Parameters
    ----------
    filename : str
        Location of the file containing the auxiliary data.
    **kwargs
       Other AuxReader options.

    See Also
    --------
    :class:`~MDAnalysis.auxiliary.base.AuxReader`

    Note
    ----
    The file is assumed to be of a size such that reading and storing the full
    contents is practical.
    """

    format = "XVG"
    _Auxstep = XVGStep

    def __init__(self, filename, **kwargs):
        self._auxdata = os.path.abspath(filename)
        with anyopen(filename) as xvg_file:
            lines = xvg_file.readlines()
        auxdata_values = []
        # remove comments before storing
        for i, line in enumerate(uncomment(lines)):
           if line.lstrip()[0] == '&':
               # multiple data sets not supported; stop at the end of the first
               break
           auxdata_values.append([float(l) for l in line.split()])
           # check the number of columns is consistent
           if len(auxdata_values[i]) != len(auxdata_values[0]):
                raise ValueError('Step {0} has {1} columns instead of '
                                 '{2}'.format(i, auxdata_values[i],
                                              auxdata_values[0]))
        self._auxdata_values = np.array(auxdata_values)
        self._n_steps = len(self._auxdata_values)
        super(XVGReader, self).__init__(**kwargs)

    def _read_next_step(self):
        """ Read next auxiliary step and update ``auxstep``.

        Returns
        -------
        AuxStep object
            Updated with the data for the new step.

        Raises
        ------
        StopIteration
            When end of auxiliary data set is reached.
        """
        auxstep = self.auxstep
        new_step = self.step + 1
        if new_step < self.n_steps:
            auxstep._data = self._auxdata_values[new_step]
            auxstep.step = new_step
            return auxstep
        else:
            self.rewind()
            raise StopIteration

    def _go_to_step(self, i):
        """ Move to and read i-th auxiliary step.

        Parameters
        ----------
        i: int
            Step number (0-indexed) to move to

        Returns
        -------
        :class:`XVGStep`

        Raises
        ------
        ValueError
            If step index not in valid range.
        """
        if i >= self.n_steps or i < 0:
            raise ValueError("Step index {0} is not valid for auxiliary "
                             "(num. steps {1})".format(i, self.n_steps))
        self.auxstep.step = i-1
        self.next()
        return self.auxstep

    def read_all_times(self):
        """ Get list of time at each step.

        Returns
        -------
        list of float
            Time at each step.
        """
        return self._auxdata_values[:,self.time_selector]


class XVGFileReader(base.AuxFileReader):
    """ Auxiliary reader to read (step at a time) from an .xvg file.

    An alternative XVG reader which reads each step from the .xvg file as
    needed (rather than reading and storing all from the start), for a lower
    memory footprint.

    Parameters
    ----------
    filename : str
       Location of the file containing the auxiliary data.
    **kwargs
       Other AuxReader options.

    See Also
    --------
    :class:`~MDAnalysis.auxiliary.base.AuxFileReader`


    Note
    ----
    The default reader for .xvg files is :class:`XVGReader`.
    """

    format = 'XVG-F'
    _Auxstep = XVGStep

    def __init__(self, filename, **kwargs):
        super(XVGFileReader, self).__init__(filename, **kwargs)

    def _read_next_step(self):
        """ Read next recorded step in xvg file and update ``austep``.

        Returns
        -------
        AuxStep object
            Updated with the data for the new step.

        Raises
        ------
        StopIteration
            When end of file or end of first data set is reached.
        """
        line = next(self.auxfile)
        while True:
            if not line or (line.strip() and line.strip()[0] == '&'):
                # at end of file or end of first set of data (multiple sets
                # currently not supported)
                self.rewind()
                raise StopIteration
            # uncomment the line
            for uncommented in uncomment([line]):
                # line has data in it; add to auxstep + return
                auxstep = self.auxstep
                auxstep.step = self.step + 1
                auxstep._data = [float(i) for i in uncommented.split()]
                # see if we've set n_cols yet...
                try:
                    auxstep._n_cols
                except AttributeError:
                    # haven't set n_cols yet; set now
                    auxstep._n_cols = len(auxstep._data)
                if len(auxstep._data) != auxstep._n_cols:
                    raise ValueError('Step {0} has {1} columns instead of '
                                     '{2}'.format(self.step, len(auxstep._data),
                                                  auxstep._n_cols))
                return auxstep
            # line is comment only - move to next
            line = next(self.auxfile)

    def _count_n_steps(self):
        """ Iterate through all steps to count total number.

        Returns
        -------
        int
            Total number of steps
        """
        if not self.constant_dt:
            # check if we've already iterated through to build _times list
            try:
                return len(self._times)
            except AttributeError:
                # might as well build _times now, since we'll need to iterate
                # through anyway
                self._times = self.read_all_times()
                return len(self.read_all_times())
        else:
            # don't need _times; iterate here instead
            self._restart()
            count = 0
            for step in self:
                count = count + 1
            return count

    def read_all_times(self):
        """ Iterate through all steps to build times list.

        Returns
        -------
        list of float
            Time of each step
        """
        self._restart()
        times = []
        for step in self:
            times.append(self.time)
        return np.array(times)
