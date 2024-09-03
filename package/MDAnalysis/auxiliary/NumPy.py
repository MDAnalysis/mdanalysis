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
NumPy auxiliary reader --- :mod:`MDAnalysis.auxiliary.NumPy`
========================================================
.. versionadded:: 2.4.0

Background
----------


Classes
-------

.. autoclass:: NumPyReader
   :members:


The actual data for each step is stored by instances of NumPyStep.

.. autoclass:: NumPyStep
   :members:


"""
from pathlib import Path
import warnings
from typing import Optional, Union, Dict, List

import numpy as np

from . import base
from .. import units


class NumPyStep(base.AuxStep):
    """:class:`AuxStep` class for NumPy arrays.

    Extends the base AuxStep class to allow selection of time and
    data-of-interest fields (by dictionary key) from the full set of data read
    each step.

    Parameters
    ----------
    time_selector : str, optional
        Name of the dictionary key that links to the time values (assumed to
        be in ps). Default value is "Time"
    data_selector : str | list of str | None, optional
        List of dictionary keys linking to data of interest to be stored in
        ``data``. Default value is ``None``.
    **kwargs
        Other AuxStep options.

    See Also
    --------
    :class:`MDAnalysis.auxiliary.base.AuxStep`
    """

    def __init__(self, time_selector: str = "Time",
                 data_selector: Optional[str] = None, **kwargs):
        super(NumPyStep, self).__init__(time_selector=time_selector,
                                      data_selector=data_selector,
                                      **kwargs)

    def _select_time(self, key: str) -> np.float64:


    def _select_data(self, key: Union[str, None]) -> np.float64:
        if key is None:
            return
        try:
            return self._data[key]
        except KeyError:
            raise KeyError(f"'{key}' is not a key in the data_dict dictionary."
                           " Check the NumPyReader.terms attribute")


class NumPyReader(base.AuxReader):
    """ Auxiliary reader to handle NumPy array data. 

    Parameters
    ----------
    filename : str
        Location of the file containing the auxiliary data.
    **kwargs
       Other AuxReader options.

    Attributes
    ----------
    _auxdata : pathlib.PosixPath
        path at which the auxiliary data file is located
    data_dict : dict
        dictionary that contains the auxiliary data
    _n_steps : int
        Number of steps for which auxdata is available
    terms : list
        Names of the auxiliary data entries available in `data_dict`.

    See Also
    --------
    :class:`MDAnalysis.auxiliary.base.AuxReader`
    :meth:`MDAnalysis.coordinates.base.ReaderBase.add_auxiliary`
    """

    format = "NumPy"
    _Auxstep = NumPyStep

    def __init__(self, filename: str, convert_units: bool = True, **kwargs):
        self._auxdata = Path(filename).resolve()
        self.data_dict = 
        self._n_steps = len(self.data_dict["Time"])
        # attribute to communicate found energy terms to user
        self.terms = list(self.data_dict.keys())
        super(NumPyReader, self).__init__(**kwargs)

    def _memory_usage(self):
        size = 0
        for array in self.data_dict.values():
            size += array.nbytes
        return size

    def _read_next_step(self) -> NumPyStep:
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
            auxstep._data = {term: self.data_dict[term][self.step + 1]
                             for term in self.terms}
            auxstep.step = new_step
            return auxstep
        else:
            self.rewind()
            raise StopIteration

    def _go_to_step(self, i: int) -> NumPyStep:
        """ Move to and read i-th auxiliary step.

        Parameters
        ----------
        i: int
            Step number (0-indexed) to move to

        Returns
        -------
        :class:`NumPyStep`

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

    def read_all_times(self) -> np.ndarray:
        """ Get list of time at each step.

        Returns
        -------
        NumPy array of float
            Time at each step.
        """
        return self.data_dict[self.time_selector]

    def get_data(self, data_selector: Union[str, List[str], None] = None
                 ) -> Dict[str, np.ndarray]:
        """ Returns the auxiliary data contained in the :class:`NumPyReader`.
        Returns either all data or data specified as `data_selector` in form
        of a str or a list of any of :attr:`NumPyReader.terms`. `Time` is
        always returned to allow easy plotting.

        Parameters
        ----------
        data_selector: str, List[str], None
            Keys to be extracted from the auxiliary reader's data dictionary.
            If ``None``, returns all data found in :attr:`.data_dict`.
        Returns
        -------
        data_dict : dict
            Dictionary mapping `data_selector` keys to NumPy arrays of the
            auxiliary data.

        Raises
        ------
        KeyError
            if an invalid data_selector key is passed.
        """
        if data_selector is None:
            return self.data_dict

        def _get_data_term(term, datadict):
            try:
                return datadict[term]
            except KeyError:
                raise KeyError(f"data selector {term} is invalid. Check the "
                               "NumPyReader's `terms` attribute.")

        data_dict = {"Time": self.data_dict["Time"]}

        if isinstance(data_selector, list):
            for term in data_selector:
                data_dict[term] = _get_data_term(term, self.data_dict)
        else:
            term = data_selector
            data_dict[term] = _get_data_term(term, self.data_dict)

        return data_dict
