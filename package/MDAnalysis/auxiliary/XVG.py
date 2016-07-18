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
XVG auxiliary reader --- :mod:`MDAnalysis.auxiliary.XVG`
========================================================

.. autofunction:: uncomment

.. autoclass:: XVGReader
   :members:

.. autoclass:: XVGFileReader
   :members:

"""

import os
import numpy as np
from . import base

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
        # '@' must be at the beginning of a line to be a grace instruction
        if line.lstrip()[0] == '@':
            continue
        # '#' can be anywhere in the line, everything after is a comment
        comment_position = line.find('#')
        if comment_position > 0 and line[:comment_position]:
            yield line[:comment_position]
        elif comment_position < 0 and line:
            yield line
        # if comment_position == 0, then the line is empty

class XVGStep(base.AuxStep):
    """ AuxStep class for .xvg file format 

    Parameters
    ----------
    time_selector : int
        Index of column in .xvg file storing time, assumed to be in ps. Default
        value is 0. 
    data_selector : list of int
        Indices of columns in .xvg file containing data of interest to be 
        stored in ``data``. Default value is ``None``, which will select 
        all columns but that containing time.
    """
    def _select_time(self, key):
        if key is None:
            # here so that None is a valid value; just return
            return
        if isinstance(key, int):
            return self._select_data(key)
        else:
             raise ValueError('Time selector must be single index')

    def _select_data(self, key):
        if key is None:
             key = [i for i in range(len(self._data)) if i != self._time_selector]
        if isinstance(key, int):
            try:
                return self._data[key]
            except IndexError:
                raise ValueError('{} not a valid index for data with {} '
                                 'columns'.format(key, len(self._data)))
        else:
            return np.array([self._select_data(i) for i in key])

    def _empty_data(self):
        return np.array([np.nan]*len(self.data))


class XVGReader(base.AuxReader):
    """ Auxiliary reader to read data from an .xvg file.

    All data from the file will be read and stored on initialisation. This is
    the default reader for .xvg files.
    
    xvg files are produced by Gromacs during simulation or analysis, formatted
    for plotting data with Grace.

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
    The time of each step is assumed to stored in the first column, in units of 
    ps, and data is assumed to be time-ordered.

    The file is assumed to be of a size such that reading and storing the full 
    contents is practical.
    """

    format = "XVG"
    _Auxstep = XVGStep

    def __init__(self, filename, time_selector=0, data_selector=None, **kwargs):
        self._data_input = os.path.abspath(filename)
        with open(filename) as xvg_file:
            lines = xvg_file.readlines()
        auxdata = []
        # remove comments before storing
        for i, line in enumerate(uncomment(lines)):
           auxdata.append([float(l) for l in line.split()])
           # check the number of columns is consistent
           if len(auxdata[i]) != len(auxdata[0]):
                raise ValueError('Step {0} has {1} columns instead of '
                                 '{2}'.format(i, auxdata[i], auxdata[0]))
        self._auxdata = np.array(auxdata)
        self._n_steps = len(self._auxdata)
        super(XVGReader, self).__init__(time_selector=time_selector, 
                                        data_selector=data_selector, **kwargs)

    def _read_next_step(self):
        """ Read next auxiliary step. """
        auxstep = self.auxstep
        new_step = self.step + 1
        if new_step < self.n_steps:
            auxstep._data = self._auxdata[new_step]
            auxstep.step = new_step
            return auxstep
        else:
            self.rewind()
            raise StopIteration

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
        """
        if i not in range(self.n_steps):
            raise ValueError("Step index {0} is not valid for auxiliary "
                             "(num. steps {1})".format(i, self.n_steps))
        self.auxstep.step = i-1
        self.next()
        return self.auxstep

    def read_all_times(self):
        """ Get list of time of each step.

        Returns
        -------
        list of float
            Time of each step
        """
        return self._auxdata[:,self.time_selector]


class XVGFileReader(base.AuxFileReader):
    """ Auxiliary reader to read (step at a time) from an .xvg file.

    An alternative XVG reader which reads each step from the .xvg file as 
    needed (rather than reading + storing all from the start), for a lower 
    memory footprint.     
    
    Parameters
    ----------
    filename : str
       Location of the file containing the auxiliary data.
    **kwargs
       Other AuxReader options.    

    See Also
    --------
    :class:`XVGReader`

    Note
    ----
    The default reader for .xvg files is :class:`XVGReader`.

    The time of each step is assumed to stored in the first column, in units of 
    ps, and data is assumed to be time-ordered.
    """

    format = 'XVG-F'
    _Auxstep = XVGStep

    def __init__(self, filename, time_selector=0, **kwargs):
        super(XVGFileReader, self).__init__(filename, time_selector=time_selector, 
                                            **kwargs)

    def _read_next_step(self):
        """ Read next recorded step in xvg file. """
        line = self.auxfile.readline()
        if line:
            auxstep = self.auxstep
            # xvg has both comments '#' and grace instructions '@'
            while line.lstrip()[0] in ['#', '@']:
                line = self.auxfile.readline()
            # remove end of line comments
            line_no_comment = line.split('#')[0]
            auxstep.step = self.step + 1
            auxstep._data = [float(i) for i in line_no_comment.split()]
            try: 
                n_cols = auxstep._n_cols
            except AttributeError:
                auxstep._n_cols = len(auxstep._data)
            if len(auxstep._data) != auxstep._n_cols:
                raise ValueError('Step {0} has {1} columns instead of '
                                 '{2}'.format(self.step, len(auxstep._data),
                                             auxstep._n_cols))
            return auxstep
        else:
            self.rewind()
            raise StopIteration

    def _count_n_steps(self):
        """ Iterate through all steps to count total number.

        Returns
        -------
        int
            Total number of steps
        """
        if self.constant_dt:
            # will have to iterate through all to built times list anyway
            return len(self.read_all_times())
        else:
            # iterate here instead
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
