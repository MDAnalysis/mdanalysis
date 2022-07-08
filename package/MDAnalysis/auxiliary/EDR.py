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

EDR files are binary files following the XDR protocol. They are written by
GROMACS during simulations and contain  the time-series energy data of the
system.

panedr is a Python package ( https://github.com/mdanalysis/panedr ) that reads
these binary files and returns them human-readable form, either as a Pandas
DataFrame or as a dictionary of NumPy arrays. It is used by the EDR auxiliary
reader to parse EDR files.

The EDR auxiliary reader takes the output from panedr and loads the energy data
as auxiliary data into Universes. Standalone usage is also possible, where the
energy terms are extracted without associating them with the trajectory, for
example, to allow easy plotting of the energy terms.

"""
import numbers
import os
import numpy as np
from . import base
from ..lib.util import anyopen
import panedrlite as panedr


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

    def __init__(self, time_selector=0, data_selector=None, **kwargs):
        super(EDRStep, self).__init__(time_selector=time_selector,
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
                errmsg = (f'{key} not a valid index for data with '
                          f'{len(self._data)} columns')
                raise ValueError(errmsg) from None
        else:
            return np.array([self._select_data(i) for i in key])


class EDRReader(base.AuxReader):
    """ Auxiliary reader to read data from a .edr file.

    Detault reader for .edr files. All data from the file will be read and
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
        self._auxfile = os.path.abspath(filename)
        self.auxdata = panedr.edr_to_dict(filename)
        self._n_steps = len(self.auxdata["Time"])
        self.terms = [key for key in self.auxdata.keys()]
        super(EDRReader, self).__init__(**kwargs)

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
            auxstep._data = self.auxdata[new_step]
            auxstep.step = new_step
            return auxstep
        else:
            self.rewind()
            raise StopIteration

