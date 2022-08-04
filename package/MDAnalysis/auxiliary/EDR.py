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
.. versionadded:: 2.3.0

EDR_ files are binary files following the XDR_ protocol. They are written by
GROMACS during simulations and contain the time-series energy data of the
system.

pyedr_ is a Python package that reads EDR binary files and returns them
human-readable form as a dictionary of NumPy arrays. It is used by the EDR
auxiliary reader to parse EDR files. As such, a dictionary with string keys and
numpy array values is loaded into the :class:`EDRReader`. It is basically a
Python-based version of the C++ code in GROMACS_.

The EDR auxiliary reader takes the output from pyedr and loads the energy data
as auxiliary data into :class:`MDAnalysis.core.universe.Universe`. Standalone
usage is also possible, where the energy terms are extracted without
associating them with the trajectory, for example, to allow easy plotting of
the energy terms.

For adding EDR data, 2-3 arguments need to be provided to
:func:`MDAnalysis.coordinates.base.add_auxiliary`: `auxname`, `auxdata`, and
(optionally) `auxterm`. `auxterm`s, if provided need to exactly match the names
of EDR terms to be added. When an auxterm is provided, `auxname` can be freely
chosen. If an `auxterm` is not provided, `auxname` must instead exactly match
EDR terms. `auxdata` can be either the path to an EDR file, or an instance of
:class:`MDAnalysis.auxiliary.EDR.EDRReader`.

Examples::

    import MDAnalysis as mda
    from MDAnalysisTests.datafiles import AUX_EDR, AUX_EDR_TPR, AUX_EDR_XTC
    import matplotlib.pyplot as plt

The EDRReader can be used to add auxiliary data from EDR files to
trajectories::

    # Create a Universe
    In [1]: u = mda.Universe(AUX_EDR_TPR, AUX_EDR_XTC)
    # Open the EDR file with an EDRReader
    In [2]: aux = mda.auxiliary.EDR.EDRReader(AUX_EDR)
    # See which terms are present in the file
    In [3]: aux.terms
    Out[3]: ['Time', 'Bond', 'Angle', ...]
    # add the Potential data under the name epot
    In [4]: u.trajectory.add_auxiliary("epot", aux, "Potential")
    # add the Angle energy without stating an `auxterm`
    In [5]: u.trajectory.add_auxiliary("Angle", aux)
    # access the value of Potential and Angle for a timestep
    In [6]: u.trajectory.ts.aux.epot
    Out[6]: -525164.0625
    In [7]: u.trajectory.ts.aux.Angle
    Out[7]: 3764.52734375
    # dictionary style accession also works
    In [8]: u.trajectory.ts.aux["epot"]
    Out[8]: -525164.0625


Note: Some `auxterms` have spaces in their names. Unless an `auxname` without a
space is provided, these terms will not be accessible via the attribute syntax.
Only the dictionary syntax will work in that case.


Adding multiple or all terms at the same time is possible as such::

    # add list of values
    u.trajectory.add_auxiliary(["Bond", "Angle"], aux)
    # add list of values, provide an `auxterm`
    u.trajectory.add_auxiliary(["epot", "ekin"], aux, ["Potential",
                                                       "Kinetic En."]
    # add all values
    u.trajectory.add_auxiliary("*", aux)



The EDRReader can also provide the data independently of trajectories. This is
useful, for example, for plotting. The data for all terms, a list of terms, or
a single term can be returned in dictionary form. "Time" is always returned in
this dictionary to make plotting easier::

    # extract temperature data from EDR file
    temp = aux.return_data("Temperature")
    plt.plot(temp["Time"], temp["Temperature"])

    # extract list of terms:
    some_terms = aux.return_data(["Potential", "Kinetic En.", "Box-X"])

    # extract all terms:
    all_terms = aux.return_data()
    # alternatively
    all_terms = aux.return_data("*")

.. _EDR: https://manual.gromacs.org/current/reference-manual/file-formats.html#edr
.. _XDR: https://datatracker.ietf.org/doc/html/rfc1014
.. _pyedr: https://github.com/mdanalysis/panedr
.. _GROMACS: https://github.com/gromacs/gromacs/blob/main/src/gromacs/fileio/enxio.cpp


.. autoclass:: EDRReader
   :members:


"""
from pathlib import Path
from . import base
import pyedr
import numpy as np
from typing import Optional, Union, Dict, List


class EDRStep(base.AuxStep):
    """:class:`AuxStep` class for .edr file format.

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

    def __init__(self, time_selector: str = "Time",
                 data_selector: Optional[str] = None, **kwargs):
        super(EDRStep, self).__init__(time_selector=time_selector,
                                      data_selector=data_selector,
                                      **kwargs)

    def _select_time(self, key: str) -> np.float64:
        """'Time' is one of the entries in the dict returned by pyedr.
        The base AuxStep Class uses the time_selector 'Time' to return the
        time value of each step."""
        return self._select_data(key)

    def _select_data(self, key: Union[str, None]) -> np.float64:
        if key is None:
            return
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
    :func:`~MDAnalysis.coordinates.base.add_auxiliary`

    Note
    ----
    The file is assumed to be of a size such that reading and storing the full
    contents is practical.
    """

    format = "EDR"
    _Auxstep = EDRStep

    def __init__(self, filename: str, **kwargs):
        self._auxdata = Path(filename).resolve()
        self.auxdata = pyedr.edr_to_dict(filename)
        self._n_steps = len(self.auxdata["Time"])
        # attribute to communicate found energy terms to user
        self.terms = [key for key in self.auxdata.keys()]
        super(EDRReader, self).__init__(**kwargs)

    def _read_next_step(self) -> EDRStep:
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

    def _go_to_step(self, i: int) -> EDRStep:
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

    def read_all_times(self) -> np.ndarray:
        """ Get list of time at each step.

        Returns
        -------
        NumPy array of float
            Time at each step.
        """
        return self.auxdata[self.time_selector]

    def return_data(self, data_selector: Union[str, List[str], None] = None) \
            -> Dict[str, np.ndarray]:
        """ Returns the auxiliary data contained in the :class:`EDRReader`.
        Returns either all data or data specified as `data_selector` in form
        of a str or a list of any of :attribute:`EDRReader.terms`. `Time` is
        always returned to allow easy plotting. """
        if not data_selector or data_selector == "*":
            return self.auxdata
        elif isinstance(data_selector, list):
            data_dict = {"Time": self.auxdata["Time"]}
            for term in data_selector:
                if term not in self.terms:
                    raise KeyError(f"data_selector {term} invalid. Check the "
                                   "EDRReader's `terms` attribute.")
                data_dict[term] = self.auxdata[term]
            return data_dict
        elif isinstance(data_selector, str):
            if data_selector not in self.terms:
                raise KeyError(f"data_selector {data_selector} invalid. Check "
                               "the EDRReader's `terms` attribute.")
            data_dict = {"Time": self.auxdata["Time"],
                         data_selector: self.auxdata[data_selector]}
            return data_dict
        else:
            raise ValueError(f"data_selector of type {type(data_selector)} is "
                             "not supported. Use list or str to indicate valid"
                             " terms. Check the EDRReader's `terms` "
                             "attribute.")

    def calc_representative(self) -> Dict[str, np.ndarray]:
        """ Calculate representative auxiliary value(s) from the data in
        *frame_data*.
        Overloaded here to accommodate the different data type. Now, this works
        for energy data dictionaries.


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
            # we don't know the original sign, and might have two
            # equally-spaced steps; check the earlier time first
            try:
                value = cutoff_data[-min_diff]
            except KeyError:
                value = cutoff_data[min_diff]
        elif self.represent_ts_as == 'average':
            value = {}
            for dataset in cutoff_data:
                for term in self.terms:
                    if term not in value:
                        value[term] = cutoff_data[dataset][term]
                    else:
                        value[term] += cutoff_data[dataset][term]
            for term in value:
                value[term] = value[term] / len(cutoff_data)
        return value
