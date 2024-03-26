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
.. versionadded:: 2.4.0

Background
----------

`EDR files`_ are binary files following the XDR_ protocol. They are written by
GROMACS during simulations and contain time-series non-trajectory data on the
system, such as energies, temperature, or pressure.

pyedr_ is a Python package that reads EDR binary files and returns them
human-readable form as a dictionary of NumPy arrays. It is used by the EDR
auxiliary reader to parse EDR files. As such, a dictionary with string keys and
numpy array values is loaded into the :class:`EDRReader`. It is basically a
Python-based version of the C++ code in GROMACS_.


The classes in this module are based on the `pyedr`_ package. Pyedr is an
optional dependency and must be installed to use this Reader. Use of the reader
without pyedr installed will raise an `ImportError`. The variable `HAS_PYEDR`
indicates whether this module has pyedr availble.

The EDR auxiliary reader takes the output from pyedr and makes this data
available within MDAnalysis. The usual workflow starts with creating an
EDRReader and passing it the file to read as such::

    import MDAnalysis as mda
    aux = mda.auxiliary.EDR.EDRReader(some_edr_file)

The newly created `aux` object contains all the data found in the EDR file. It
is stored in the :attr:`.data_dict` dictionary, which maps the names GROMACS
gave each data entry to a NumPy array that holds the relevant data. These
GROMACS-given names are stored in and available through the :attr:`.terms`
attribute. In addition to the numeric data, the new `EDRReader` also stores the
units of each entry in the :attr:`.data_dict` dictionary in its
:attr:`.unit_dict` dictionary.

.. warning::
    Units are converted to `MDAnalysis base units`_ automatically unless
    otherwise specified. However, not all unit types have a defined base unit
    in MDAnalysis. (cf. :data:`MDAnalysis.units.MDANALYSIS_BASE_UNITS`).
    Pressure units, for example, are not currently defined, and
    will not be converted. This might cause inconsistencies between units!
    Conversion can be switched off by passing `convert_units=False` when the
    EDRReader is created::

        aux = mda.auxiliary.EDR.EDRReader(some_edr_file, convert_units=False)


Standalone Usage of the EDRReader
---------------------------------

The :class:`.EDRReader` can be used to access the data stored in EDR files on
its own, without association of data to trajectory frames.
This is useful, for example, for plotting. The data for a single term, a list
of terms, or for all terms can be returned in dictionary form. "Time" is always
returned in this dictionary to make plotting easier::


    temp = aux.get_data("Temperature")
    plt.plot(temp["Time"], temp["Temperature"])

    some_terms = aux.get_data(["Potential", "Kinetic En.", "Box-X"])
    plt.plot(some_terms["Time"], some_terms["Potential"])

    all_terms = aux.get_data()
    plt.plot(all_terms["Time"], all_terms["Pressure"])

Adding EDR data to trajectories
-------------------------------

Like other AuxReaders, the :class:`.EDRReader` can attach its data to a
trajectory by associating it to the appropriate time steps.
In general, to add EDR data to a trajectory, one needs to provide two
arguments.

.. note::
    The following will change with the release of MDAnalysis 3.0. From then on,
    the order of these two arguments will be reversed.

The first argument is the `aux_spec` dictionary. With it, users specify which
entries from the EDR file they want to add, and they give it a more convenient
name to be used in MDAnalysis (because GROMACS creates names like
"#Surf*SurfTen" or "'Constr. rmsd'" which may be inconvenient to use.)
This dictionary might look like this::

    aux_spec = {"epot": "Potential",
                "surf_tension": "#Surf*SurfTen"}

When provided as shown below, this would direct the :class:`.EDRReader` to take
the data it finds under the "Potential" key in its :attr:`.data_dict`
dictionary and attach it to the trajectory time steps under
`u.trajectory.ts.aux.epot` (and the same for the surface tension).


The second argument needed is the source of the EDR data itself. Here, either
the path to the EDR file or a previously created :class:`.EDRReader` object
can be provided.


Examples::

    import MDAnalysis as mda
    from MDAnalysisTests.datafiles import AUX_EDR, AUX_EDR_TPR, AUX_EDR_XTC
    import matplotlib.pyplot as plt

A :class:`Universe` and an :class:`.EDRReader` object are created and the data
available in the EDR file is printed::

    In [1]: u = mda.Universe(AUX_EDR_TPR, AUX_EDR_XTC)
    In [2]: aux = mda.auxiliary.EDR.EDRReader(AUX_EDR)
    In [3]: aux.terms
    Out[3]: ['Time', 'Bond', 'Angle', ...]

Data is associated with the trajectory, using an `aux_spec` dictionary to
specify which data to add under which name. Any number of terms can be added
using this method. The data is then accessible in the `ts.aux` namespace via
both attribute and dictionary syntax::

    In [4]: u.trajectory.add_auxiliary(aux,
                                        {"epot": "Potential",
                                        "angle": "Angle"}, 
                                        )
    In [5]: u.trajectory.ts.aux.epot
    Out[5]: -525164.0625
    In [6]: u.trajectory.ts.aux.Angle
    Out[6]: 3764.52734375
    In [7]: u.trajectory.ts.aux["epot"]
    Out[7]: -525164.0625


.. note::
    Some GROMACS-provided :attr:`terms` have spaces. Unless an attribute name
    without a space is provided, these terms will not be accessible via the
    attribute syntax. Only the dictionary syntax will work in that case.


Further, it is possible to add all data from the EDR file to the trajectory. To
do this, the `aux_spec` dictionary is omitted, and the data source (the second
argument as explained above) is provided explicitly as `auxdata`. When adding
data this way, the terms in :attr:`.terms` become the names used in `ts.aux`::

    In [7]: u.trajectory.add_auxiliary(auxdata=aux)
    In [8]: u.trajectory.ts.aux["#Surf*SurfTen"]
    Out[8]: -1857.519287109375


.. _EDR files: https://manual.gromacs.org/current/reference-manual/file-formats.html#edr
.. _XDR: https://datatracker.ietf.org/doc/html/rfc1014
.. _pyedr: https://github.com/mdanalysis/panedr
.. _GROMACS: https://github.com/gromacs/gromacs/blob/main/src/gromacs/fileio/enxio.cpp
.. _MDAnalysis base units: https://docs.mdanalysis.org/2.3.0/documentation_pages/units.html


Classes
-------

.. autoclass:: EDRReader
   :members:


The actual data for each step is stored by instances of EDRStep.

.. autoclass:: EDRStep
   :members:


"""
from pathlib import Path
import warnings
from typing import Optional, Union, Dict, List

import numpy as np

from . import base
from .. import units

try:
    import pyedr
except ImportError:
    # Indicates whether pyedr is found
    HAS_PYEDR = False
else:
    # Indicates whether pyedr is found
    HAS_PYEDR = True


class EDRStep(base.AuxStep):
    """:class:`AuxStep` class for the .edr file format.

    Extends the base AuxStep class to allow selection of time and
    data-of-interest fields (by dictionary key) from the full set of data read
    each step.

    Parameters
    ----------
    time_selector : str, optional
        Name of the dictionary key that links to the time values (assumed to
        be in ps). Default value is "Time"
    data_selector : str | list of str | None, optional
        List of dictionary keys linking to data of interest in the EDR file to
        be stored in ``data``. Default value is ``None``.
    **kwargs
        Other AuxStep options.

    See Also
    --------
    :class:`MDAnalysis.auxiliary.base.AuxStep`
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
            raise KeyError(f"'{key}' is not a key in the data_dict dictionary."
                           " Check the EDRReader.terms attribute")


class EDRReader(base.AuxReader):
    """ Auxiliary reader to read data from an .edr file.

    `EDR files`_
    are created by GROMACS during a simulation. They are binary files which
    contain time-series energy data and other data related to the simulation.

    Default reader for .edr files. All data from the file will be read and
    stored on initialisation.

    Parameters
    ----------
    filename : str
        Location of the file containing the auxiliary data.
    convert_units : bool, optional
        If True (default), units from the EDR file are automatically converted
        to MDAnalysis base units. If False, units are taken from the file
        as-is. Where no base unit is defined in MDAnalysis, no conversion takes
        place. Unit types in :data:`MDAnalysis.units.MDANALYSIS_BASE_UNITS`
        will be converted automatically by default.
    **kwargs
       Other AuxReader options.

    Attributes
    ----------
    _auxdata : pathlib.PosixPath
        path at which the auxiliary data file is located
    data_dict : dict
        dictionary that contains the auxiliary data, mapping the names GROMACS
        gave the entries in the EDR file to a NumPy array containing this data
    unit_dict : dict
        dictionary that contains the units of the auxiliary data, mapping the
        :attr:`data_selector` of the Reader (i.e. the name of the dataset in
        the EDR file) to its unit.
    _n_steps : int
        Number of steps for which auxdata is available
    terms : list
        Names of the auxiliary data entries available in `data_dict`. These are
        the names GROMACS set in the EDR file.

    See Also
    --------
    :class:`MDAnalysis.auxiliary.base.AuxReader`
    :meth:`MDAnalysis.coordinates.base.ReaderBase.add_auxiliary`

    Note
    ----
    The file is assumed to be of a size such that reading and storing the full
    contents is practical. A warning will be issued when memory usage exceeds
    1 GB. This warning limit can be changed via the ``memory_limit`` kwarg.
    """

    format = "EDR"
    _Auxstep = EDRStep

    def __init__(self, filename: str, convert_units: bool = True, **kwargs):
        if not HAS_PYEDR:
            raise ImportError("EDRReader: To read EDR files please install "
                              "pyedr.")
        self._auxdata = Path(filename).resolve()
        self.data_dict = pyedr.edr_to_dict(filename)
        self.unit_dict = pyedr.get_unit_dictionary(filename)
        self.convert_units = convert_units
        if self.convert_units:
            self._convert_units()
        self._n_steps = len(self.data_dict["Time"])
        # attribute to communicate found energy terms to user
        self.terms = list(self.data_dict.keys())
        super(EDRReader, self).__init__(**kwargs)

    def _convert_units(self):
        """Called during :func:`__init__` to convert the units found in the EDR
        file to MDAnalysis base units"""
        unknown_units = []
        for term, unit in self.unit_dict.items():
            try:
                unit_type = units.unit_types[unit]
            except KeyError:
                if unit not in unknown_units:
                    unknown_units.append(unit)
                continue  # skip conversion if unit not defined yet

            target_unit = units.MDANALYSIS_BASE_UNITS[unit_type]
            data = self.data_dict[term]
            self.data_dict[term] = units.convert(data, unit, target_unit)
            self.unit_dict[term] = units.MDANALYSIS_BASE_UNITS[unit_type]
        if unknown_units:
            warnings.warn("Could not find unit type for the following "
                          f"units: {unknown_units}")

    def _memory_usage(self):
        size = 0
        for array in self.data_dict.values():
            size += array.nbytes
        return size

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
            auxstep._data = {term: self.data_dict[term][self.step + 1]
                             for term in self.terms}
            auxstep.step = new_step
            return auxstep
        else:
            self.rewind()
            if self.n_steps > 1:
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
        return self.data_dict[self.time_selector]

    def get_data(self, data_selector: Union[str, List[str], None] = None
                 ) -> Dict[str, np.ndarray]:
        """ Returns the auxiliary data contained in the :class:`EDRReader`.
        Returns either all data or data specified as `data_selector` in form
        of a str or a list of any of :attr:`EDRReader.terms`. `Time` is
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
                               "EDRReader's `terms` attribute.")

        data_dict = {"Time": self.data_dict["Time"]}

        if isinstance(data_selector, list):
            for term in data_selector:
                data_dict[term] = _get_data_term(term, self.data_dict)
        else:
            term = data_selector
            data_dict[term] = _get_data_term(term, self.data_dict)

        return data_dict
