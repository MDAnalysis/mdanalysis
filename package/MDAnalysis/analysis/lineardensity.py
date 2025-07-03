# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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
Linear Density --- :mod:`MDAnalysis.analysis.lineardensity`
===========================================================

A tool to compute mass and charge density profiles along the three
cartesian axes [xyz] of the simulation cell. Works only for orthorombic,
fixed volume cells (thus for simulations in canonical NVT ensemble).
"""
import os.path as path

import numpy as np
import warnings

from MDAnalysis.analysis.base import AnalysisBase, Results
from MDAnalysis.units import constants
from MDAnalysis.lib.util import deprecate
from MDAnalysis.analysis.results import ResultsGroup


# TODO: Remove in version 3.0.0
class Results(Results):
    """From version 3.0.0 onwards, some entries in Results will be renamed. See
    the docstring for LinearDensity for details. The Results class is defined
    here to implement deprecation warnings for the user."""

    _deprecation_dict = {
        "pos": "mass_density",
        "pos_std": "mass_density_stddev",
        "char": "charge_density",
        "char_std": "charge_density_stddev",
    }

    def _deprecation_warning(self, key):
        warnings.warn(
            f"`{key}` is deprecated and will be removed in version 3.0.0. "
            f"Please use `{self._deprecation_dict[key]}` instead.",
            DeprecationWarning,
        )

    def __getitem__(self, key):
        if key in self._deprecation_dict.keys():
            self._deprecation_warning(key)
            return super(Results, self).__getitem__(
                self._deprecation_dict[key]
            )
        return super(Results, self).__getitem__(key)

    def __getattr__(self, attr):
        if attr in self._deprecation_dict.keys():
            self._deprecation_warning(attr)
            attr = self._deprecation_dict[attr]
        return super(Results, self).__getattr__(attr)


class LinearDensity(AnalysisBase):
    r"""Linear density profile

    Parameters
    ----------
    select : AtomGroup
          any atomgroup
    grouping : str {'atoms', 'residues', 'segments', 'fragments'}
          Density profiles will be computed either on the atom positions (in
          the case of 'atoms') or on the center of mass of the specified
          grouping unit ('residues', 'segments', or 'fragments').
    binsize : float
          Bin width in Angstrom used to build linear density
          histograms. Defines the resolution of the resulting density
          profile (smaller --> higher resolution)
    verbose : bool, optional
          Show detailed progress of the calculation if set to ``True``

    Attributes
    ----------
    results.x.dim : int
           index of the [xyz] axes
    results.x.mass_density : numpy.ndarray
           mass density in :math:`g \cdot cm^{-3}` in [xyz] direction
    results.x.mass_density_stddev : numpy.ndarray
           standard deviation of the mass density in [xyz] direction
    results.x.charge_density : numpy.ndarray
           charge density in :math:`\mathrm{e} \cdot mol \cdot cm^{-3}` in
           [xyz] direction
    results.x.charge_density_stddev : numpy.ndarray
           standard deviation of the charge density in [xyz] direction
    results.x.pos: numpy.ndarray
        Alias to the :attr:`results.x.mass_density` attribute.

        .. deprecated:: 2.2.0
           Will be removed in MDAnalysis 3.0.0. Please use
           :attr:`results.x.mass_density` instead.
    results.x.pos_std: numpy.ndarray
        Alias to the :attr:`results.x.mass_density_stddev` attribute.

        .. deprecated:: 2.2.0
           Will be removed in MDAnalysis 3.0.0. Please use
           :attr:`results.x.mass_density_stddev` instead.
    results.x.char: numpy.ndarray
        Alias to the :attr:`results.x.charge_density` attribute.

        .. deprecated:: 2.2.0
           Will be removed in MDAnalysis 3.0.0. Please use
           :attr:`results.x.charge_density` instead.
    results.x.char_std: numpy.ndarray
        Alias to the :attr:`results.x.charge_density_stddev` attribute.

        .. deprecated:: 2.2.0
           Will be removed in MDAnalysis 3.0.0. Please use
           :attr:`results.x.charge_density_stddev` instead.
    results.x.slice_volume : float
           volume of bin in [xyz] direction
    results.x.hist_bin_edges : numpy.ndarray
           edges of histogram bins for mass/charge densities, useful for, e.g.,
           plotting of histogram data.
    Note: These density units are likely to be changed in the future.

    Example
    -------
    First create a :class:`LinearDensity` object by supplying a selection,
    then use the :meth:`run` method. Finally access the results
    stored in results, i.e. the mass density in the x direction.

    .. code-block:: python

       ldens = LinearDensity(selection)
       ldens.run()
       print(ldens.results.x.mass_density)


    Alternatively, other types of grouping can be selected using the
    ``grouping`` keyword. For example to calculate the density based on
    a grouping of the :class:`~MDAnalysis.core.groups.ResidueGroup`
    of the input :class:`~MDAnalysis.core.groups.AtomGroup`.

    .. code-block:: python

       ldens = LinearDensity(selection, grouping='residues', binsize=1.0)
       ldens.run()



    .. versionadded:: 0.14.0

    .. versionchanged:: 1.0.0
       Support for the ``start``, ``stop``, and ``step`` keywords has been
       removed. These should instead be passed to :meth:`LinearDensity.run`.
       The ``save()`` method was also removed, you can use ``np.savetxt()`` or
       ``np.save()`` on the :attr:`LinearDensity.results` dictionary contents
       instead.

    .. versionchanged:: 1.0.0
       Changed `selection` keyword to `select`

    .. versionchanged:: 2.0.0
       Results are now instances of
       :class:`~MDAnalysis.core.analysis.Results` allowing access
       via key and attribute.

    .. versionchanged:: 2.2.0

       *  Fixed a bug that caused LinearDensity to fail if grouping="residues"
          or grouping="segments" were set.
       *  Residues, segments, and fragments will be analysed based on their
          centre of mass, not centre of geometry as previously stated.
       *  LinearDensity now works with updating atom groups.
       *  Added new result container :attr:`results.x.hist_bin_edges`.
          It contains the bin edges of the histrogram bins for calculated
          densities and can be used for easier plotting of histogram data.

    .. versionchanged:: 2.10.0
       *  Introduced :meth:`get_supported_backends` allowing for parallel execution
          on :mod:`multiprocessing` and :mod:`dask` backends.
       *  Removed undocumented and unused attribute :attr:`totalmass`.

    .. deprecated:: 2.2.0
       The `results` dictionary has been changed and the attributes
       :attr:`results.x.pos`, :attr:`results.x.pos_std`, :attr:`results.x.char`
       and :attr:`results.x.char_std` are now deprecated. They will be removed
       in 3.0.0. Please use :attr:`results.x.mass_density`,
       :attr:`results.x.mass_density_stddev`, :attr:`results.x.charge_density`,
       and :attr:`results.x.charge_density_stddev` instead.
    """

    _analysis_algorithm_is_parallelizable = True

    @classmethod
    def get_supported_backends(cls):
        return (
            "serial",
            "multiprocessing",
            "dask",
        )

    def __init__(self, select, grouping="atoms", binsize=0.25, **kwargs):
        super(LinearDensity, self).__init__(
            select.universe.trajectory, **kwargs
        )
        # allows use of run(parallel=True)
        self._ags = [select]
        self._universe = select.universe

        self.binsize = binsize

        # group of atoms on which to compute the COM (same as used in
        # AtomGroup.wrap())
        self.grouping = grouping

        # Initiate result instances
        self.results["x"] = Results(dim=0)
        self.results["y"] = Results(dim=1)
        self.results["z"] = Results(dim=2)

        # Box sides
        self.dimensions = self._universe.dimensions[:3]
        self.volume = np.prod(self.dimensions)
        # number of bins
        bins = (self.dimensions // self.binsize).astype(int)

        # Here we choose a number of bins of the largest cell side so that
        # x, y and z values can use the same "coord" column in the output file
        self.nbins = bins.max()
        slices_vol = self.volume / bins

        self.keys = [
            "mass_density",
            "mass_density_stddev",
            "charge_density",
            "charge_density_stddev",
        ]

        # Initialize results array with zeros
        for dim in self.results:
            idx = self.results[dim]["dim"]
            self.results[dim]["slice_volume"] = slices_vol[idx]
            for key in self.keys:
                self.results[dim][key] = np.zeros(self.nbins)

        # Get masses and charges for the selection (e.g. UpdatingAtomGroup)
        if self.grouping == "atoms":
            self.masses = self._ags[0].masses
            self.charges = self._ags[0].charges

        elif self.grouping in ["residues", "segments", "fragments"]:
            self.masses = self._ags[0].total_mass(compound=self.grouping)
            self.charges = self._ags[0].total_charge(compound=self.grouping)

        else:
            raise AttributeError(
                f"{self.grouping} is not a valid value for grouping."
            )

    @staticmethod
    def _custom_aggregator(results):
        # NB: the *stddev values here are not the standard deviation,
        # but the variance. The stddev is calculated in _conclude()
        mass_density = np.sum(
            [entry["mass_density"] for entry in results], axis=0
        )
        mass_density_stddev = np.sum(
            [entry["mass_density_stddev"] for entry in results], axis=0
        )
        charge_density = np.sum(
            [entry["charge_density"] for entry in results], axis=0
        )
        charge_density_stddev = np.sum(
            [entry["charge_density_stddev"] for entry in results], axis=0
        )
        return Results(
            dim=results[0]["dim"],
            slice_volume=results[0]["slice_volume"],
            hist_bin_edges=results[0]["hist_bin_edges"],
            mass_density=mass_density,
            mass_density_stddev=mass_density_stddev,
            charge_density=charge_density,
            charge_density_stddev=charge_density_stddev,
        )

    def _get_aggregator(self):
        return ResultsGroup(
            lookup={
                "x": self._custom_aggregator,
                "y": self._custom_aggregator,
                "z": self._custom_aggregator,
            }
        )

    def _single_frame(self):
        if self.grouping == "atoms":
            self.masses = self._ags[0].masses
            self.charges = self._ags[0].charges

        elif self.grouping in ["residues", "segments", "fragments"]:
            self.masses = self._ags[0].total_mass(compound=self.grouping)
            self.charges = self._ags[0].total_charge(compound=self.grouping)

        else:
            raise AttributeError(
                f"{self.grouping} is not a valid value for grouping."
            )

        self.group = getattr(self._ags[0], self.grouping)
        self._ags[0].wrap(compound=self.grouping)
        # Find position of atom/group of atoms
        if self.grouping == "atoms":
            positions = self._ags[0].positions  # faster for atoms
        else:
            # Centre of mass for residues, segments, fragments
            positions = self._ags[0].center_of_mass(compound=self.grouping)

        for dim in ["x", "y", "z"]:
            idx = self.results[dim]["dim"]

            key = "mass_density"
            key_std = "mass_density_stddev"
            # histogram for positions weighted on masses
            hist, _ = np.histogram(
                positions[:, idx],
                weights=self.masses,
                bins=self.nbins,
                range=(0.0, max(self.dimensions)),
            )

            self.results[dim][key] += hist
            self.results[dim][key_std] += np.square(hist)

            key = "charge_density"
            key_std = "charge_density_stddev"
            # histogram for positions weighted on charges
            hist, bin_edges = np.histogram(
                positions[:, idx],
                weights=self.charges,
                bins=self.nbins,
                range=(0.0, max(self.dimensions)),
            )

            self.results[dim][key] += hist
            self.results[dim][key_std] += np.square(hist)
            self.results[dim]["hist_bin_edges"] = bin_edges

    def _conclude(self):
        avogadro = constants["N_Avogadro"]  # unit: mol^{-1}
        volume_conversion = 1e-24  # unit: A^3/cm^3
        # divide result values by avodagro and convert from A3 to cm3
        k = avogadro * volume_conversion

        # Average results over the number of configurations
        for dim in ["x", "y", "z"]:
            for key in [
                "mass_density",
                "mass_density_stddev",
                "charge_density",
                "charge_density_stddev",
            ]:
                self.results[dim][key] /= self.n_frames
            # Compute standard deviation for the error
            # For certain tests in testsuite, floating point imprecision
            # can lead to negative radicands of tiny magnitude (yielding nan).
            # radicand_mass and radicand_charge are therefore calculated first
            # and negative values set to 0 before the square root
            # is calculated.
            radicand_mass = self.results[dim][
                "mass_density_stddev"
            ] - np.square(self.results[dim]["mass_density"])
            radicand_mass[radicand_mass < 0] = 0
            self.results[dim]["mass_density_stddev"] = np.sqrt(radicand_mass)

            radicand_charge = self.results[dim][
                "charge_density_stddev"
            ] - np.square(self.results[dim]["charge_density"])
            radicand_charge[radicand_charge < 0] = 0
            self.results[dim]["charge_density_stddev"] = np.sqrt(
                radicand_charge
            )

        for dim in ["x", "y", "z"]:
            # norming factor, units of mol^-1 cm^3
            norm = k * self.results[dim]["slice_volume"]
            for key in self.keys:
                self.results[dim][key] /= norm

    # TODO: Remove in 3.0.0
    @deprecate(
        release="2.2.0",
        remove="3.0.0",
        message="It will be replaced by a :meth:`_reduce` "
        "method in the future",
    )
    def _add_other_results(self, other):
        """For parallel analysis"""
        for dim in ["x", "y", "z"]:
            for key in self.keys:
                self.results[dim][key] += other.results[dim][key]
