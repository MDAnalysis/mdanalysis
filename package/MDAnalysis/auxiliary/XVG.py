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

.. autoclass:: XVGReader
   :members:

"""

import numpy as np
from . import base

class XVGReader(base.AuxFileReader):
    """ Auxiliary reader to read (step at a time) from an .xvg file.

    xvg files are produced by Gromacs during simulation or analysis, formatted
    for plotting data with Grace.
    
    Paramaters
    ----------
    filename : str
       Location of the file containing the auxiliary data.
    **kwargs
       Other AuxReader options.    

    Notes
    -----
    The time of each step is assumed to stored in the first column (*time_col*
    defaults to 0), and data is assumed to be time-ordered.
    """

    # TODO - swtich to reading file all at once

    format = "XVG"
 
    def __init__(self, filename, **kwargs):
        time_col = kwargs.pop('time_col', 0)
        super(XVGReader, self).__init__(filename, time_col=time_col, 
                                        **kwargs)

        
    def _read_next_step(self):
        """ Read next recorded step in xvg file. """
        line = self.auxfile.readline()
        if line:
            # xvg has both comments '#' and grace instructions '@'
            while line.lstrip()[0] in ['#', '@']:
                line = self.auxfile.readline()
            # remove end of line comments
            line_no_comment = line.split('#')[0]
            self.step = self.step + 1
            self._data = [float(i) for i in line_no_comment.split()]
            if self.n_cols and len(self._data) != self.n_cols:
                raise ValueError('Step {0} has {1} columns instead of '
                                 '{2}'.format(self.step, len(self._data),
                                             self.n_cols))
            return self._data
        else:
            self.go_to_first_step()
            raise StopIteration


    def count_n_steps(self):
        """ Iterate through all steps to count total number and make times list.

        Returns
        -------
        int
           Total number of steps
        """
        self._restart()
        times = []
        count = 0
        for step in self:
            count = count + 1
            times.append(self.time)
        self._times = times
        return count

    def read_all_times(self):
        self.count_n_steps()
