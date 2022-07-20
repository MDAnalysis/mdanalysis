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
EDR auxiliary reader --- :mod:`MDAnalysis.auxiliary.EDR`
========================================================

EDR files are binary files following the XDR protocol
(https://datatracker.ietf.org/doc/html/rfc1014). They are written by
GROMACS during simulations and contain the time-series energy data of the
system.

pyedr is a Python package ( https://github.com/mdanalysis/panedr ) that reads
EDR binary files and returns them human-readable form as a dictionary of NumPy
arrays. It is used by the EDR auxiliary reader to parse EDR files. As such, a
dictionary with string keys and numpy array values is loaded into the
EDRReader.

The EDR auxiliary reader takes the output from pyedr and loads the energy data
as auxiliary data into :class:`~MDAnalysis.core.universe.Universe`. Standalone
usage is also possible, where the energy terms are extracted without
associating them with the trajectory, for example, to allow easy plotting of
the energy terms.



"""
import os
from . import base
import pyedr


class EDRStep(base.AuxStep):
    """ AuxStep class for .edr file format.

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

    def __init__(self, time_selector="Time", data_selector=None, **kwargs):
        super(EDRStep, self).__init__(time_selector=time_selector,
                                      data_selector=data_selector,
                                      **kwargs)

    def _select_time(self, key):
        """'Time' is one of the entries in the dict returned by pyedr.
        The base AuxStep Class uses the time_selector 'Time' to return the
        time value of each step."""
        return self._select_data(key)

    def _select_data(self, key):
        try:
            return self._data[key]
        except KeyError:
            raise KeyError(f"'{key}' is not a key in the auxdata dictionary. "
                           "Check the EDRReader.terms attribute")


class EDRReader(base.AuxReader):
    """ Auxiliary reader to read data from a .edr file.

    Default reader for .edr files. All data from the file will be read and
    stored on initialisation.

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

    format = "EDR"
    _Auxstep = EDRStep

    def __init__(self, filename, **kwargs):
        self._auxdata = os.path.abspath(filename)
        self.auxdata = pyedr.edr_to_dict(filename)
        self._n_steps = len(self.auxdata["Time"])
        # attribute to communicate found energy terms to user
        self.terms = [key for key in self.auxdata.keys()]
        super(EDRReader, self).__init__(**kwargs)

    def _read_next_step(self):
        """Read next auxiliary step and update ``auxstep``.

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
            auxstep._data = {term: self.auxdata[term][self.step + 1]
                             for term in self.terms}
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
        :class:`EDRStep`

        Raises
        ------
        ValueError
            If step index not in valid range.
        """
        if i >= self.n_steps or i < 0:
            raise ValueError("Step index {0} is not valid for auxiliary "
                             "(num. steps {1})".format(i, self.n_steps))
        self.auxstep.step = i - 1
        self.next()
        return self.auxstep

    def read_all_times(self):
        """ Get list of time at each step.

        Returns
        -------
        NumPy array of float
            Time at each step.
        """
        return self.auxdata[self.time_selector]
