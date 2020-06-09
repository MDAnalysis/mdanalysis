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

# MDAnalysis -- density analysis
# Copyright (c) 2007-2011 Oliver Beckstein <orbeckst@gmail.com>
# (based on code from Hop --- a framework to analyze solvation dynamics from MD simulations)

r"""Generating densities from trajectories --- :mod:`MDAnalysis.analysis.density`
=============================================================================

:Author: Oliver Beckstein
:Year: 2011
:Copyright: GNU Public License v3

The module provides classes and functions to generate and represent
volumetric data, in particular densities.

Generating a density from a MD trajectory
-----------------------------------------

A common use case is to analyze the solvent density around a protein of
interest. The density is calculated with :func:`density_from_Universe` in the
fixed coordinate system of the simulation unit cell. It is therefore necessary
to orient and fix the protein with respect to the box coordinate system. In
practice this means centering and superimposing the protein, frame by frame, on
a reference structure and translating and rotating all other components of the
simulation with the protein. In this way, the solvent will appear in the
reference frame of the protein.

An input trajectory must

1. have been centered on the protein of interest;
2. have all molecules made whole that have been broken across periodic
   boundaries [#pbc]_;
3. have the solvent molecules remapped so that they are closest to the
   solute (this is important when using triclinic unit cells such as
   a dodecahedron or a truncated octahedron) [#pbc]_.
4. have a fixed frame of reference; for instance, by superimposing a protein
   on a reference structure so that one can study the solvent density around
   it [#fit]_.

To generate the density of water molecules around a protein (assuming that the
trajectory is already appropriately treated for periodic boundary artifacts and
is suitably superimposed to provide a fixed reference frame) [#testraj]_ ::

  from MDAnalysis.analysis.density import DensityAnalysis
  u = Universe(TPR, XTC)
  ow = u.select_atoms("name OW")
  D = DensityAnalysis(ow, delta=1.0)
  D.run()
  D.density.convert_density('TIP4P')
  D.density.export("water.dx", type="double")

The positions of all water oxygens (the :class:`AtomGroup` `ow`) are
histogrammed on a grid with spacing *delta* = 1 Å. Initially the density is
measured in :math:`\text{Å}^{-3}`. With the :meth:`Density.convert_density`
method, the units of measurement are changed. In the example we are now
measuring the density relative to the literature value of the TIP4P water model
at ambient conditions (see the values in :data:`MDAnalysis.units.water` for
details). Finally, the density is written as an OpenDX_ compatible file that
can be read in VMD_, Chimera_, or PyMOL_.

The :class:`Density` object is accessible as the
:attr:`DensityAnalysis.density` attribute.  In particular, the data for the
density is stored as a NumPy array in :attr:`Density.grid`, which can be
processed in any manner.


Creating densities
------------------

The :class:`DensityAnalysis` class generates a :class:`Density` from an
atomgroup. Its function equivalent :func:`density_from_Universe` is deprecated.

:func:`density_from_PDB` generates a pseudo-density from a PDB file by
replacing each atom with a Gaussian density with a width that is computed from
the crystallographic temperature factors (B-factor) with :func:`Bfactor2RMSF`.

.. autoclass:: DensityAnalysis
   :members:
   :inherited-members: run

   .. attribute:: density

      After the analysis (see the :meth:`~DensityAnalysis.run` method), the resulting density is
      stored in the :attr:`density` attribute as a :class:`Density` instance.


Density object
--------------

The main output of the density creation functions is a :class:`Density`
instance, which is derived from a :class:`gridData.core.Grid`. A
:class:`Density` is essentially a 3D array with origin and lengths.

.. See Also:: :mod:`gridData`


.. autoclass:: Density
   :members:
   :inherited-members:
   :show-inheritance:


Deprecated functionality
------------------------

Use :class:`DensityAnalysis` instead of :func:`density_from_Universe`. The
:func:`density_from_PDB` function is no longer supported and will also be
removed in 2.0.0.

.. autofunction:: density_from_Universe

.. autofunction:: density_from_PDB

.. autofunction:: Bfactor2RMSF

.. autoclass:: BfactorDensityCreator
   :members:

.. autofunction:: notwithin_coordinates_factory



.. rubric:: Footnotes

.. [#pbc] Making molecules whole can be accomplished with the
          :meth:`MDAnalysis.core.groups.AtomGroup.wrap` of
          :attr:`Universe.atoms` (use ``compound="fragments"``).  or the
          PBC-wrapping transformations in
          :mod:`MDAnalysis.transformations.wrap`.

.. [#fit] Superposition can be performed with
          :class:`MDAnalysis.analysis.align.AlignTraj` or the fitting
          transformations in :mod:`MDAnalysis.transformations.fit`.

.. [#testraj] Note that the trajectory in the example (`XTC`) is *not*
              properly made whole and fitted to a reference structure;
              these steps were omitted to clearly show the steps necessary
              for the actual density calculation.

.. Links
.. -----

.. _OpenDX: http://www.opendx.org/
.. _VMD:   http://www.ks.uiuc.edu/Research/vmd/
.. _Chimera: https://www.cgl.ucsf.edu/chimera/
.. _PyMOL: http://www.pymol.org/
.. _Gromacs: http://www.gromacs.org
.. _`gmx trjconv`: http://manual.gromacs.org/programs/gmx-trjconv.html

"""

from __future__ import print_function, division, absolute_import
from six.moves import range, zip
from six import raise_from, string_types

import numpy as np
import sys
import os
import os.path
import errno
import warnings

from gridData import Grid

import MDAnalysis
from MDAnalysis.core import groups
from MDAnalysis.lib.util import fixedwidth_bins, iterable, asiterable
from MDAnalysis.lib import NeighborSearch as NS
from MDAnalysis import NoDataError, MissingDataWarning
from .. import units
from ..lib import distances
from MDAnalysis.lib.log import ProgressBar 
from MDAnalysis.lib.log import ProgressMeter   # remove in 2.0
from MDAnalysis.lib.util import deprecate

from .base import AnalysisBase

import logging

logger = logging.getLogger("MDAnalysis.analysis.density")

class DensityAnalysis(AnalysisBase):
    r"""Volumetric density analysis.

    The trajectory is read, frame by frame, and the atoms in `atomgroup` are
    histogrammed on a 3D grid with spacing `delta`.

    Parameters
    ----------
    atomgroup : AtomGroup or UpdatingAtomGroup
            Group of atoms (such as all the water oxygen atoms) being analyzed.
            This can be an :class:`~MDAnalysis.core.groups.UpdatingAtomGroup` for
            selections that change every time step.
    delta : float (optional)
            Bin size for the density grid in ångström (same in x,y,z).
    padding : float (optional)
            Increase histogram dimensions by padding (on top of initial box
            size) in ångström. Padding is ignored when setting a user defined
            grid.
    gridcenter : numpy ndarray, float32 (optional)
            3 element numpy array detailing the x, y and z coordinates of the
            center of a user defined grid box in ångström.
    xdim : float (optional)
            User defined x dimension box edge in ångström; ignored if
            gridcenter is "None".
    ydim : float (optional)
            User defined y dimension box edge in ångström; ignored if
            gridcenter is "None".
    zdim : float (optional)
            User defined z dimension box edge in ångström; ignored if
            gridcenter is "None".

    See Also
    --------
    pmda.density.DensityAnalysis for a parallel version

    Notes
    -----
    Normal :class:`AtomGroup` instances represent a static selection of
    atoms. If you want *dynamically changing selections* (such as "name OW and
    around 4.0 (protein and not name H*)", i.e., the water oxygen atoms that
    are within 4 Å of the protein heavy atoms) then create an
    :class:`~MDAnalysis.core.groups.UpdatingAtomGroup` (see Examples).


    Examples
    --------
    A common use case is to analyze the solvent density around a protein of
    interest. The density is calculated with :class:`DensityAnalysis` in the
    fixed coordinate system of the simulation unit cell. It is therefore
    necessary to orient and fix the protein with respect to the box coordinate
    system. In practice this means centering and superimposing the protein,
    frame by frame, on a reference structure and translating and rotating all
    other components of the simulation with the protein. In this way, the
    solvent will appear in the reference frame of the protein.

    An input trajectory must

    1. have been centered on the protein of interest;
    2. have all molecules made whole that have been broken across periodic
       boundaries [#pbc]_;
    3. have the solvent molecules remapped so that they are closest to the
       solute (this is important when using triclinic unit cells such as
       a dodecahedron or a truncated octahedron) [#pbc]_;
    4. have a fixed frame of reference; for instance, by superimposing a
       protein on a reference structure so that one can study the solvent
       density around it [#fit]_.

    .. rubric:: Generate the density

    To generate the density of water molecules around a protein (assuming that
    the trajectory is already appropriately treated for periodic boundary
    artifacts and is suitably superimposed to provide a fixed reference frame)
    [#testraj]_, first  create the :class:`DensityAnalysis` object by
    supplying an AtomGroup, then use  the :meth:`run` method::

        from MDAnalysis.analysis import density
        u = Universe(TPR, XTC)
        ow = u.select_atoms("name OW")
        D = density.DensityAnalysis(ow, delta=1.0)
        D.run()
        D.density.convert_density('TIP4P')

    The positions of all water oxygens are histogrammed on a grid with spacing
    *delta* = 1 Å and stored as a :class:`Density` object in the attribute
    :attr:`DensityAnalysis.density`.

    .. rubric:: Working with a density

    A :class:`Density` contains a large number of methods and attributes that
    are listed in the documentation. Here we use the
    :meth:`Density.convert_density` to convert the density from inverse cubic
    ångström to a density relative to the bulk density of TIP4P water at
    standard conditions. (MDAnalysis stores a number of literature values in
    :data:`MDAnalysis.units.water`.)

    One can directly access the density as a 3D NumPy array through
    :attr:`Density.grid`.

    By default, the :class:`Density` object returned contains a physical
    density in units of Å\ :sup:`-3`. If you are interested in recovering the
    underlying **probability density**, simply divide by the sum::

      probability_density = D.density.grid / D.density.grid.sum()

    Similarly, if you would like to recover a grid containing a **histogram of
    atom counts**, simply multiply by the volume `dV` of each bin (or voxel);
    in this case you need to ensure that the physical density is measured in
    Å\ :sup:`-3` by converting it::

      import numpy as np

      # ensure that the density is A^{-3}
      D.density.convert_density("A^{-3}")

      dV = np.prod(D.density.delta)
      atom_count_histogram = D.density.grid * dV


    .. rubric:: Writing the density to a file

    A density can be `exported to different formats
    <https://www.mdanalysis.org/GridDataFormats/gridData/formats.html>`_ with
    :meth:`Density.export` (thanks to the fact that :class:`Density` is a
    subclass :class:`gridData.core.Grid`, which provides the functionality).
    For example, to `write a DX file
    <https://www.mdanalysis.org/GridDataFormats/gridData/basic.html#writing-out-data>`_
    ``water.dx`` that can be read with VMD, PyMOL, or Chimera::

      D.density.export("water.dx", type="double")


    .. rubric:: Example: Water density in the whole simulation

    Basic use for creating a water density (just using the water oxygen
    atoms "OW")::

      D = DensityAnalysis(universe.select_atoms('name OW')).run()


    .. rubric:: Example: Water in a binding site (updating selection)

    If you are only interested in water within a certain region, e.g., within
    a vicinity around a binding site, you can use a selection that updates
    every step by using an :class:`~MDAnalysis.core.groups.UpdatingAtomGroup`::

      near_waters = universe.select_atoms('name OW and around 5 (resid 156 157 305)',
                    updating=True)
      D_site = DensityAnalysis(near_waters).run()


    .. rubric:: Example: Small region around a ligand (manual box selection)

    If you are interested in explicitly setting a grid box of a given edge size
    and origin, you can use the `gridcenter` and `xdim`/`ydim`/`zdim`
    arguments.  For example to plot the density of waters within 5 Å of a
    ligand (in this case the ligand has been assigned the residue name "LIG")
    in a cubic grid with 20 Å edges which is centered on the center of mass
    (COM) of the ligand::

      # Create a selection based on the ligand
      ligand_selection = universe.select_atoms("resname LIG")

      # Extract the COM of the ligand
      ligand_COM = ligand_selection.center_of_mass()

      # Create a density of waters on a cubic grid centered on the ligand COM
      # In this case, we update the atom selection as shown above.
      ligand_waters = universe.select_atoms('name OW and around 5 resname LIG',
                                            updating=True)
      D_water = DensityAnalysis(ligand_waters,
                                delta=1.0,
                                gridcenter=ligand_COM,
                                xdim=20, ydim=20, zdim=20)

    (It should be noted that the `padding` keyword is not used when a user
    defined grid is assigned).



    .. versionadded:: 1.0.0

    """
    def __init__(self, atomgroup, delta=1.0,
                 metadata=None, padding=2.0,
                 gridcenter=None,
                 xdim=None, ydim=None, zdim=None):
        u = atomgroup.universe
        super(DensityAnalysis, self).__init__(u.trajectory)
        self._atomgroup = atomgroup
        self._delta = delta
        self._padding = padding
        self._gridcenter = gridcenter
        self._xdim = xdim
        self._ydim = ydim
        self._zdim = zdim

    def _prepare(self):
        coord = self._atomgroup.positions
        if self._gridcenter is not None:
            # Issue 2372: padding is ignored, defaults to 2.0 therefore warn
            if self._padding > 0:
                msg = ("Box padding (currently set at {0}) "
                       "is not used in user defined grids.".format(self._padding))
                warnings.warn(msg)
                logger.warning(msg)
            # Generate a copy of smin/smax from coords to later check if the
            # defined box might be too small for the selection
            smin = np.min(coord, axis=0)
            smax = np.max(coord, axis=0)
            # Overwrite smin/smax with user defined values
            smin, smax = _set_user_grid(self._gridcenter, self._xdim,
                                        self._ydim, self._zdim, smin, smax)
        else:
            # Make the box bigger to avoid as much as possible 'outlier'. This
            # is important if the sites are defined at a high density: in this
            # case the bulk regions don't have to be close to 1 * n0 but can
            # be less. It's much more difficult to deal with outliers.  The
            # ideal solution would use images: implement 'looking across the
            # periodic boundaries' but that gets complicated when the box
            # rotates due to RMS fitting.
            smin = np.min(coord, axis=0) - self._padding
            smax = np.max(coord, axis=0) + self._padding
        BINS = fixedwidth_bins(self._delta, smin, smax)
        arange = np.transpose(np.vstack((BINS['min'], BINS['max'])))
        bins = BINS['Nbins']
        # create empty grid with the right dimensions (and get the edges)
        grid, edges = np.histogramdd(np.zeros((1, 3)), bins=bins,
                                     range=arange, normed=False)
        grid *= 0.0
        self._grid = grid
        self._edges = edges
        self._arange = arange
        self._bins = bins
        self.density = None

    def _single_frame(self):
        h, _ = np.histogramdd(self._atomgroup.positions,
                              bins=self._bins, range=self._arange,
                              normed=False)
        # reduce (proposed change #2542 to match the parallel version in pmda.density)
        # return self._reduce(self._grid, h)
        #
        # serial code can simply do
        self._grid += h

    def _conclude(self):
        # average:
        self._grid /= float(self.n_frames)
        density = Density(grid=self._grid, edges=self._edges,
                          units={'length': "Angstrom"},
                          parameters={'isDensity': False})
        density.make_density()
        self.density = density

    # _reduce is not strictly necessary for the serial version but is necessary for
    # pmda-style parallelism (see #2542)
    # @staticmethod
    # def _reduce(res, result_single_frame):
    #     """'accumulate' action for a time series
    #
    #     If `res` is a numpy array, the `result_single_frame` is added to it
    #     element-wise. If `res` and `result_single_frame` are lists then
    #     `result_single_frame` is appended to `res`.
    #     """
    #     if isinstance(res, list) and len(res) == 0:
    #         res = result_single_frame
    #     else:
    #         res += result_single_frame
    #     return res


class Density(Grid):
    r"""Class representing a density on a regular cartesian grid.

    Parameters
    ----------
    grid : array_like
        histogram or density, typically a :class:`numpy.ndarray`
    edges : list
        list of arrays, the lower and upper bin edges along the axes
    parameters : dict
        dictionary of class parameters; saved with
        :meth:`Density.save`. The following keys are meaningful to
        the class. Meaning of the values are listed:

         *isDensity*

            - ``False``: grid is a histogram with counts [default]
            - ``True``: a density

            Applying :meth:`Density.make_density`` sets it to ``True``.
    units : dict
        A dict with the keys

        - *length*:  physical unit of grid edges (Angstrom or nm) [Angstrom]
        - *density*: unit of the density if ``isDensity=True`` or ``None``
          otherwise; the default is "Angstrom^{-3}" for densities
          (meaning :math:`\text{Å}^{-3}`).
    metadata : dict
        a user defined dictionary of arbitrary values associated with the
        density; the class does not touch :attr:`Density.metadata` but
        stores it with :meth:`Density.save`

    Attributes
    ----------
    grid : array
        counts or density
    edges : list of 1d-arrays
        The boundaries of each cell in `grid` along all axes (equivalent
        to what :func:`numpy.histogramdd` returns).
    delta : array
        Cell size in each dimension.
    origin : array
        Coordinates of the *center* of the cell at index `grid[0, 0, 0, ...,
        0]`, which is considered to be the front lower left corner.
    units : dict
        The units for lengths and density; change units with the method
        :meth:`~Density.convert_length` or :meth:`~Density.convert_density`.


    Notes
    -----
    The data (:attr:`Density.grid`) can be manipulated as a standard numpy
    array. Changes can be saved to a file using the :meth:`Density.save` method. The
    grid can be restored using the :meth:`Density.load` method or by supplying the
    filename to the constructor.

    The attribute :attr:`Density.metadata` holds a user-defined dictionary that
    can be used to annotate the data. It is also saved with :meth:`Density.save`.

    The :meth:`Density.export` method always exports a 3D object (written in
    such a way to be readable in VMD_, Chimera_, and PyMOL_), the rest should
    work for an array of any dimension. Note that PyMOL_ only understands DX
    files with the DX data type "double" in the "array" object (see `known
    issues when writing OpenDX files`_ and issue
    `MDAnalysis/GridDataFormats#35`_ for details). Using the keyword
    ``type="double"`` for the method :meth:`Density.export`, the user can
    ensure that the DX file is written in a format suitable for PyMOL_.

    If the input histogram consists of counts per cell then the
    :meth:`Density.make_density` method converts the grid to a physical density. For
    a probability density, divide it by :meth:`Density.grid.sum` or use ``normed=True``
    right away in :func:`~numpy.histogramdd`.

    The user *should* set the *parameters* keyword (see docs for the
    constructor); in particular, if the data are already a density, one must
    set ``isDensity=True`` because there is no reliable way to detect if
    data represent counts or a density. As a special convenience, if data are
    read from a file and the user has not set ``isDensity`` then it is assumed
    that the data are in fact a density.

    .. _`MDAnalysis/GridDataFormats#35`:
       https://github.com/MDAnalysis/GridDataFormats/issues/35
    .. _`known issues when writing OpenDX files`:
       https://www.mdanalysis.org/GridDataFormats/gridData/formats/OpenDX.html#known-issues-for-writing-opendx-files

    See Also
    --------
    gridData.core.Grid is the base class of :class:`Density`.

    Examples
    --------
    Typical use:

    1. From a histogram (i.e. counts on a grid)::

        h,edges = numpy.histogramdd(...)
        D = Density(h, edges, parameters={'isDensity': False}, units={'length': 'A'})
        D.make_density()

    2. From a saved density file (e.g. in OpenDX format), where the lengths are
       in Angstrom and the density in 1/A**3::

         D = Density("density.dx")

    3. From a saved density file (e.g. in OpenDX format), where the lengths are
       in Angstrom and the density is measured relative to the density of water
       at ambient conditions::

         D = Density("density.dx", units={'density': 'water'})

    4. From a saved *histogram* (less common, but in order to demonstrate the
       *parameters* keyword) where the lengths are in nm::

         D = Density("counts.dx", parameters={'isDensity': False}, units={'length': 'nm'})
         D.make_density()
         D.convert_length('Angstrom^{-3}')
         D.convert_density('water')

       After the final step, ``D`` will contain a density on a grid measured in
       ångström, with the density values itself measured relative to the
       density of water.

    :class:`Density` objects can be algebraically manipulated (added,
    subtracted, multiplied, ...)  but there are *no sanity checks* in place to
    make sure that units, metadata, etc are compatible!


    .. Note::

       It is suggested to construct the Grid object from a histogram,
       to supply the appropriate length unit, and to use
       :meth:`Density.make_density` to obtain a density. This ensures
       that the length- and the density unit correspond to each other.

    """

    def __init__(self, *args, **kwargs):
        length_unit = 'Angstrom'

        parameters = kwargs.pop('parameters', {})
        if len(args) > 0 and isinstance(args[0], string_types) or isinstance(kwargs.get('grid', None), string_types):
            # try to be smart: when reading from a file then it is likely that this
            # is a density
            parameters.setdefault('isDensity', True)
        else:
            parameters.setdefault('isDensity', False)
        units = kwargs.pop('units', {})
        units.setdefault('length', length_unit)
        if parameters['isDensity']:
            units.setdefault('density', length_unit)
        else:
            units.setdefault('density', None)

        super(Density, self).__init__(*args, **kwargs)

        self.parameters = parameters  # isDensity: set by make_density()
        self.units = units

    def _check_set_unit(self, u):
        """Check and set units.

        First check that all units and their values in the dict `u` are valid
        and then set the object's units attribute.

        Parameters
        ----------
        u : dict
            ``{unit_type : value, ...}``

        Raises
        ------
        ValueError
            if unit types or unit values are not recognized or if required
            unit types are not in :attr:`units`
        """
        # all this unit crap should be a class...
        try:
            for unit_type, value in u.items():
                if value is None:  # check here, too iffy to use dictionary[None]=None
                    self.units[unit_type] = None
                    continue
                try:
                    units.conversion_factor[unit_type][value]
                    self.units[unit_type] = value
                except KeyError:
                    raise_from(
                        ValueError('Unit ' + str(value) + ' of type ' + str(unit_type) + ' is not recognized.'),
                        None,
                        )
        except AttributeError:
            errmsg = '"unit" must be a dictionary with keys "length" and "density.'
            logger.fatal(errmsg)
            raise_from(ValueError(errmsg), None)
        # need at least length and density (can be None)
        if 'length' not in self.units:
            raise ValueError('"unit" must contain a unit for "length".')
        if 'density' not in self.units:
            self.units['density'] = None

    def make_density(self):
        """Convert the grid (a histogram, counts in a cell) to a density (counts/volume).

        This method changes the grid irrevocably.

        For a probability density, manually divide by :meth:`grid.sum`.

        If this is already a density, then a warning is issued and nothing is
        done, so calling `make_density` multiple times does not do any harm.
        """
        # Make it a density by dividing by the volume of each grid cell
        # (from numpy.histogramdd, which is for general n-D grids)

        if self.parameters['isDensity']:
            msg = "Running make_density() makes no sense: Grid is already a density. Nothing done."
            logger.warning(msg)
            warnings.warn(msg)
            return

        dedges = [np.diff(edge) for edge in self.edges]
        D = len(self.edges)
        for i in range(D):
            shape = np.ones(D, int)
            shape[i] = len(dedges[i])
            self.grid /= dedges[i].reshape(shape)
        self.parameters['isDensity'] = True
        # see units.densityUnit_factor for units
        self.units['density'] = self.units['length'] + "^{-3}"

    def convert_length(self, unit='Angstrom'):
        """Convert Grid object to the new `unit`.

        Parameters
        ----------
        unit : str (optional)
              unit that the grid should be converted to: one of
              "Angstrom", "nm"

        Notes
        -----
        This changes the edges but will not change the density; it is the
        user's responsibility to supply the appropriate unit if the Grid object
        is constructed from a density. It is suggested to start from a
        histogram and a length unit and use :meth:`make_density`.

        """
        if unit == self.units['length']:
            return
        cvnfact = units.get_conversion_factor('length', self.units['length'], unit)
        self.edges = [x * cvnfact for x in self.edges]
        self.units['length'] = unit
        self._update()  # needed to recalculate midpoints and origin

    def convert_density(self, unit='Angstrom'):
        """Convert the density to the physical units given by `unit`.

        Parameters
        ----------
        unit : str (optional)
             The target unit that the density should be converted to.

             `unit` can be one of the following:

             =============  ===============================================================
             name           description of the unit
             =============  ===============================================================
             Angstrom^{-3}  particles/A**3
             nm^{-3}        particles/nm**3
             SPC            density of SPC water at standard conditions
             TIP3P          ... see :data:`MDAnalysis.units.water`
             TIP4P          ... see :data:`MDAnalysis.units.water`
             water          density of real water at standard conditions (0.997 g/cm**3)
             Molar          mol/l
             =============  ===============================================================

        Raises
        ------
        RuntimeError
             If the density does not have a unit associated with it to begin
             with (i.e., is not a density) then no conversion can take place.
        ValueError
             for unknown `unit`.

        Notes
        -----

        (1) This method only works if there is already a length unit associated with the
            density; otherwise raises :exc:`RuntimeError`
        (2) Conversions always go back to unity so there can be rounding
            and floating point artifacts for multiple conversions.

        """
        if not self.parameters['isDensity']:
            errmsg = 'The grid is not a density so converty_density() makes no sense.'
            logger.fatal(errmsg)
            raise RuntimeError(errmsg)
        if unit == self.units['density']:
            return
        try:
            self.grid *= units.get_conversion_factor('density',
                                                     self.units['density'], unit)
        except KeyError:
            raise_from(
                ValueError("The name of the unit ({0!r} supplied) must be one of:\n{1!r}".format(unit, units.conversion_factor['density'].keys())),
                None,
                )
        self.units['density'] = unit

    def __repr__(self):
        if self.parameters['isDensity']:
            grid_type = 'density'
        else:
            grid_type = 'histogram'
        return '<Density ' + grid_type + ' with ' + str(self.grid.shape) + ' bins>'


def _set_user_grid(gridcenter, xdim, ydim, zdim, smin, smax):
    """Helper function to set the grid dimensions to user defined values

    Parameters
    ----------
    gridcenter : numpy ndarray, float32
            3 element ndarray containing the x, y and z coordinates of the grid
            box center
    xdim : float
            Box edge length in the x dimension
    ydim : float
            Box edge length in the y dimension
    zdim : float
            Box edge length in the y dimension
    smin : numpy ndarray, float32
            Minimum x,y,z coordinates for the input selection
    smax : numpy ndarray, float32
            Maximum x,y,z coordinates for the input selection

    Returns
    -------
    umin : numpy ndarray, float32
            Minimum x,y,z coordinates of the user defined grid
    umax : numpy ndarray, float32
            Maximum x,y,z coordinates of the user defined grid
    """
    # Check user inputs
    try:
        gridcenter = np.asarray(gridcenter, dtype=np.float32)
    except ValueError:
        raise_from(ValueError("Non-number values assigned to gridcenter"), None)
    if gridcenter.shape != (3,):
        raise ValueError("gridcenter must be a 3D coordinate")
    try:
        xyzdim = np.array([xdim, ydim, zdim], dtype=np.float32)
    except ValueError:
        raise_from(ValueError("xdim, ydim, and zdim must be numbers"), None)

    # Set min/max by shifting by half the edge length of each dimension
    umin = gridcenter - xyzdim/2
    umax = gridcenter + xyzdim/2

    # Here we test if coords of selection fall outside of the defined grid
    # if this happens, we warn users they may want to resize their grids
    if any(smin < umin) or any(smax > umax):
        msg = ("Atom selection does not fit grid --- "
               "you may want to define a larger box")
        warnings.warn(msg)
        logger.warning(msg)
    return umin, umax


@deprecate(release="1.0.0", remove="2.0.0",
           message="Use ``DensityAnalysis(u, ...).run().density`` instead.")
def density_from_Universe(universe, delta=1.0, select='name OH2',
                          start=None, stop=None, step=None,
                          metadata=None, padding=2.0, cutoff=0, soluteselection=None,
                          use_kdtree=True, update_selection=False,
                          verbose=False, interval=1, quiet=None,
                          parameters=None,
                          gridcenter=None, xdim=None, ydim=None, zdim=None):
    """Create a density grid from a :class:`MDAnalysis.Universe` object.

    The trajectory is read, frame by frame, and the atoms selected with
    `select` are histogrammed on a grid with spacing `delta`.
    A physical density of units [Angstrom^{-3}] is returned (see
    :class:`Density` for more details).

    Parameters
    ----------
    universe : MDAnalysis.Universe
            :class:`MDAnalysis.Universe` object with a trajectory
    select : str (optional)
            selection string (MDAnalysis syntax) for the species to be analyzed
            ["name OH2"]
    delta : float (optional)
            bin size for the density grid in Angstrom (same in x,y,z) [1.0]
    start : int (optional)
    stop : int (optional)
    step : int (optional)
            Slice the trajectory as ``trajectory[start:stop:step]``; default
            is to read the whole trajectory.
    metadata : dict. optional
            `dict` of additional data to be saved with the object; the meta data
            are passed through as they are.
    padding : float (optional)
            increase histogram dimensions by padding (on top of initial box size)
            in Angstrom. Padding is ignored when setting a user defined grid. [2.0]
    soluteselection : str (optional)
            MDAnalysis selection for the solute, e.g. "protein" [``None``]
    cutoff : float (optional)
            With `cutoff`, select "<atomsel> NOT WITHIN <cutoff> OF <soluteselection>"
            (Special routines that are faster than the standard ``AROUND`` selection);
            any value that evaluates to ``False`` (such as the default 0) disables this
            special selection.
    update_selection : bool (optional)
            Should the selection of atoms be updated for every step? [``False``]

            - ``True``: atom selection is updated for each frame, can be slow
            - ``False``: atoms are only selected at the beginning
    verbose : bool (optional)
            Print status update to the screen for every *interval* frame? [``True``]

            - ``False``: no status updates when a new frame is processed
            - ``True``: status update every frame (including number of atoms
              processed, which is interesting with ``update_selection=True``)
    interval : int (optional)
           Show status update every `interval` frame [1]
    parameters : dict (optional)
            `dict` with some special parameters for :class:`Density` (see docs)
    gridcenter : numpy ndarray, float32 (optional)
            3 element numpy array detailing the x, y and z coordinates of the
            center of a user defined grid box in Angstrom [``None``]
    xdim : float (optional)
            User defined x dimension box edge in Angstrom; ignored if
            gridcenter is ``None``
    ydim : float (optional)
            User defined y dimension box edge in Angstrom; ignored if
            gridcenter is ``None``
    zdim : float (optional)
            User defined z dimension box edge in Angstrom; ignored if
            gridcenter is ``None``

    Returns
    -------
    :class:`Density`
            A :class:`Density` instance with the histogrammed data together
            with associated metadata.


    Notes
    -----

    By default, the `select` is static, i.e., atoms are only selected
    once at the beginning. If you want *dynamically changing selections* (such
    as "name OW and around 4.0 (protein and not name H*)", i.e., the water
    oxygen atoms that are within 4 Å of the protein heavy atoms) then set
    ``update_selection=True``. For the special case of calculating a density of
    the "bulk" solvent away from a solute use the optimized selections with
    keywords *cutoff* and *soluteselection* (see Examples below).

    Examples
    --------
    Basic use for creating a water density (just using the water oxygen atoms "OW")::

      density = density_from_Universe(universe, delta=1.0, select='name OW')

    If you are only interested in water within a certain region, e.g., within a
    vicinity around a binding site, you can use a selection that updates every
    step by setting the `update_selection` keyword argument::

      site_density = density_from_Universe(universe, delta=1.0,
                                           select='name OW and around 5 (resid 156 157 305)',
                                           update_selection=True)

    A special case for an updating selection is to create the "bulk density",
    i.e., the water outside the immediate solvation shell of a protein: Select
    all water oxygen atoms that are *farther away* than a given cut-off (say, 4
    Å) from the solute (here, heavy atoms of the protein)::

      bulk = density_from_Universe(universe, delta=1.0, select='name OW',
                                   solute="protein and not name H*",
                                   cutoff=4)

    (Using the special case for the bulk with `soluteselection` and `cutoff`
    improves performance over the simple `update_selection` approach.)

    If you are interested in explicitly setting a grid box of a given edge size
    and origin, you can use the gridcenter and x/y/zdim arguments. For example
    to plot the density of waters within 5 Å of a ligand (in this case the
    ligand has been assigned the residue name "LIG") in a cubic grid with 20 Å
    edges which is centered on the centre of mass (COM) of the ligand::

      # Create a selection based on the ligand
      ligand_selection = universe.select_atoms("resname LIG")

      # Extract the COM of the ligand
      ligand_COM = ligand_selection.center_of_mass()

      # Generate a density of waters on a cubic grid centered on the ligand COM
      # In this case, we update the atom selection as shown above.
      water_density = density_from_Universe(universe, delta=1.0,
                                            select='name OW around 5 resname LIG',
                                            update_selection=True,
                                            gridcenter=ligand_COM,
                                            xdim=20.0, ydim=20.0, zdim=20.0)

      (It should be noted that the `padding` keyword is not used when a user
      defined grid is assigned).

    As detailed above, the :class:`Density` object returned contains a
    physical density in units of Angstrom^{-3}. If you are interested in
    recovering the underlying probability density, simply divide by the sum::

      physical_density = density_from_Universe(universe, delta=1.0,
                                               select='name OW')

      probability_density = physical_density / physical_density.grid.sum()

    Similarly, if you would like to recover a grid containing a histogram of
    atom counts, simply multiply by the volume::

      # Here we assume that numpy is imported as np
      volume = np.prod(physical_density.delta)

      atom_count_histogram = physical_density * volume


    .. versionchanged:: 0.13.0
       *update_selection* and *quiet* keywords added
    .. deprecated:: 0.16
       The keyword argument *quiet* is deprecated in favor of *verbose*.
    .. versionchanged:: 0.19.0
       *gridcenter*, *xdim*, *ydim* and *zdim* keywords added to allow for user
       defined boxes
    .. versionchanged:: 0.20.0
       ProgressMeter now iterates over the number of frames analysed.
    .. versionchanged:: 1.0.0
       time_unit and length_unit default to ps and Angstrom now flags have
       been removed (same as previous flag defaults);
       warns users that `padding` value is not used in user defined grids
    .. deprecated:: 1.0.0
       `density_from_Universe` will removed in 2.0.0; use `DensityAnalysis` instead
    """
    u = universe

    if cutoff > 0 and soluteselection is not None:
        # special fast selection for '<atomsel> not within <cutoff> of <solutesel>'
        notwithin_coordinates = notwithin_coordinates_factory(
            u, select, soluteselection, cutoff,
            use_kdtree=use_kdtree, updating_selection=update_selection)
        def current_coordinates():
            return notwithin_coordinates()
    else:
        group = u.select_atoms(select, updating=update_selection)

        def current_coordinates():
            return group.positions

    coord = current_coordinates()
    logger.info(
        "Selected {0:d} atoms out of {1:d} atoms ({2!s}) from {3:d} total."
        "".format(coord.shape[0], len(u.select_atoms(select)),
                  select, len(u.atoms))
    )

    # mild warning; typically this is run on RMS-fitted trajectories and
    # so the box information is rather meaningless
    box, angles = u.trajectory.ts.dimensions[:3], u.trajectory.ts.dimensions[3:]
    if tuple(angles) != (90., 90., 90.):
        msg = ("Non-orthorhombic unit-cell --- "
               "make sure that it has been remapped properly!")
        warnings.warn(msg)
        logger.warning(msg)

    if gridcenter is not None:
        # Issue 2372: padding is ignored, defaults to 2.0 therefore warn
        if padding > 0:
            msg = ("Box padding (currently set at {0}) "
                   "is not used in user defined grids.".format(padding))
            warnings.warn(msg)
            logger.warning(msg)
        # Generate a copy of smin/smax from coords to later check if the
        # defined box might be too small for the selection
        smin = np.min(coord, axis=0)
        smax = np.max(coord, axis=0)
        # Overwrite smin/smax with user defined values
        smin, smax = _set_user_grid(gridcenter, xdim, ydim, zdim, smin, smax)
    else:
        # Make the box bigger to avoid as much as possible 'outlier'. This
        # is important if the sites are defined at a high density: in this
        # case the bulk regions don't have to be close to 1 * n0 but can
        # be less. It's much more difficult to deal with outliers.  The
        # ideal solution would use images: implement 'looking across the
        # periodic boundaries' but that gets complicate when the box
        # rotates due to RMS fitting.
        smin = np.min(coord, axis=0) - padding
        smax = np.max(coord, axis=0) + padding

    BINS = fixedwidth_bins(delta, smin, smax)
    arange = np.vstack((BINS['min'], BINS['max']))
    arange = np.transpose(arange)
    bins = BINS['Nbins']

    # create empty grid with the right dimensions (and get the edges)
    grid, edges = np.histogramdd(np.zeros((1, 3)), bins=bins, range=arange, normed=False)
    grid *= 0.0
    h = grid.copy()

    start, stop, step = u.trajectory.check_slice_indices(start, stop, step)
    n_frames = len(range(start, stop, step))

    for ts in ProgressBar(u.trajectory[start:stop:step],
                          verbose=verbose, desc="Histogramming"):
        coord = current_coordinates()

        if len(coord) == 0:
            continue

        h[:], edges[:] = np.histogramdd(coord, bins=bins, range=arange, normed=False)
        grid += h  # accumulate average histogram
    grid /= float(n_frames)

    metadata = metadata if metadata is not None else {}
    metadata['psf'] = u.filename
    metadata['dcd'] = u.trajectory.filename
    metadata['select'] = select
    metadata['n_frames'] = n_frames
    metadata['totaltime'] = round(u.trajectory.n_frames * u.trajectory.dt, 3)
    metadata['dt'] = u.trajectory.dt
    metadata['time_unit'] = 'ps'
    try:
        metadata['trajectory_skip'] = u.trajectory.skip_timestep  # frames
    except AttributeError:
        metadata['trajectory_skip'] = 1  # seems to not be used..
    try:
        metadata['trajectory_delta'] = u.trajectory.delta  # in native units
    except AttributeError:
        metadata['trajectory_delta'] = 1
    if cutoff > 0 and soluteselection is not None:
        metadata['soluteselection'] = soluteselection
        metadata['cutoff'] = cutoff  # in Angstrom

    parameters = parameters if parameters is not None else {}
    parameters['isDensity'] = False  # must override


    g = Density(grid=grid, edges=edges, units={'length': 'Angstrom'},
                parameters=parameters, metadata=metadata)
    g.make_density()
    logger.info("Density completed (initial density in Angstrom**-3)")

    return g


@deprecate(release="1.0.0", remove="2.0.0")
def notwithin_coordinates_factory(universe, sel1, sel2, cutoff,
                                  not_within=True, use_kdtree=True, updating_selection=False):
    """Generate optimized selection for '*sel1* not within *cutoff* of *sel2*'

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe object on which to operate
    sel1 : str
        Selection string for the *solvent* selection (should be
        the group with the *larger number of atoms* to make the
        KD-Tree search more efficient)
    sel2 : str
        Selection string for the *solute* selection
    cutoff : float
        Distance cutoff
    not_within : bool

        - ``True``: selection behaves as "not within" (As described above)
        - ``False``: selection is a "<sel1> WITHIN <cutoff> OF <sel2>"
    use_kdtree : bool

        - ``True``: use fast KD-Tree based selections
        - ``False``: use distance matrix approach
    updating_selection : bool
        If ``True``, re-evaluate the selection string each frame.

    Notes
    -----
    * Periodic boundary conditions are *not* taken into account: the naive
      minimum image convention employed in the distance check is currently
      not being applied to remap the coordinates themselves, and hence it
      would lead to counts in the wrong region.
    * With ``updating_selection=True``, the selection is evaluated every turn;
      do not use distance based selections (such as "AROUND") in your selection
      string because it will likely completely negate any gains from using
      this function factory in the first place.

    Examples
    --------
    :func:`notwithin_coordinates_factory` creates an optimized function that, when called,
    returns the coordinates of the "solvent" selection that are *not within* a given cut-off
    distance of the "solute". Because it is KD-tree based, it is cheap to query the KD-tree with
    a different cut-off::

      notwithin_coordinates = notwithin_coordinates_factory(universe, 'name OH2', 'protein and not name H*', 3.5)
      ...
      coord = notwithin_coordinates()        # get coordinates outside cutoff 3.5 A
      coord = notwithin_coordinates(cutoff2) # can use different cut off

    For programmatic convenience, the function can also function as a factory for a simple
    *within cutoff* query if the keyword ``not_within=False`` is set::

      within_coordinates = notwithin_coordinates_factory(universe, 'name OH2','protein and not name H*', 3.5,
                                                         not_within=False)
      ...
      coord = within_coordinates()        # get coordinates within cutoff 3.5 A
      coord = within_coordinates(cutoff2) # can use different cut off

    (Readability is enhanced by properly naming the generated function
    ``within_coordinates()``.)


    .. deprecated:: 1.0.0
       :func:`notwithin_coordinates_factory` is no longer supported and will be
       removed in 2.0.0.
    """
    # Benchmark of FABP system (solvent 3400 OH2, protein 2100 atoms) on G4 powerbook, 500 frames
    #                    cpu/s    relative   speedup       use_kdtree
    # distance matrix    633        1          1           False
    # AROUND + kdtree    420        0.66       1.5         n/a ('name OH2 around 4 protein')
    # manual + kdtree    182        0.29       3.5         True
    solvent = universe.select_atoms(sel1, updating=updating_selection)
    protein = universe.select_atoms(sel2, updating=updating_selection)

    if use_kdtree:
        # using faster hand-coded 'not within' selection with kd-tree
        if not_within is True:  # default
            def notwithin_coordinates(cutoff=cutoff):
                # must update every time step
                ns_w = NS.AtomNeighborSearch(solvent)  # build kd-tree on solvent (N_w > N_protein)
                solvation_shell = ns_w.search(protein, cutoff)  # solvent within CUTOFF of protein
                # Find indices in solvent NOT in solvation shell
                uniq_idx = np.setdiff1d(solvent.ix, solvation_shell.ix)
                # Then reselect these from Universe.atoms (as these indices are global)
                group = universe.atoms[uniq_idx]
                return group.positions
        else:
            def notwithin_coordinates(cutoff=cutoff):
                # acts as '<solvent> WITHIN <cutoff> OF <protein>'
                # must update every time step
                ns_w = NS.AtomNeighborSearch(solvent)  # build kd-tree on solvent (N_w > N_protein)
                group = ns_w.search(protein, cutoff)  # solvent within CUTOFF of protein
                return group.positions
    else:
        # slower distance matrix based (calculate all with all distances first)
        dist = np.zeros((len(solvent), len(protein)), dtype=np.float64)
        box = None  # as long as s_coor is not minimum-image remapped
        if not_within is True:  # default
            compare = np.greater
            aggregatefunc = np.all
        else:
            compare = np.less_equal
            aggregatefunc = np.any

        def notwithin_coordinates(cutoff=cutoff):
            s_coor = solvent.positions
            p_coor = protein.positions
            # Does water i satisfy d[i,j] > r for ALL j?
            d = distances.distance_array(s_coor, p_coor, box=box, result=dist)
            return s_coor[aggregatefunc(compare(d, cutoff), axis=1)]
    return notwithin_coordinates


# This is useful in principle but it is not needed here anymore. Maybe move
# to the rms module??
@deprecate(release="1.0.0", remove="2.0.0")
def Bfactor2RMSF(B):
    r"""Atomic root mean square fluctuation (in Angstrom) from the crystallographic B-factor

    RMSF and B-factor are related by [Willis1975]_

    .. math::

        B = \frac{8\pi^2}{3} \rm{RMSF}^2

    and this function returns

    .. math::

        \rm{RMSF} = \sqrt{\frac{3 B}{8\pi^2}}

    References
    ----------

    .. [Willis1975]  BTM Willis and AW Pryor. *Thermal vibrations in crystallography*. Cambridge Univ. Press, 1975


    .. deprecated:: 1.0.0
       :func:`Bfactor2RMSF` is no longer supported and will be removed in 2.0.0.
       as part of the removal of the :func:`density_from_PDB` function.
    """
    return np.sqrt(3. * B / 8.) / np.pi


@deprecate(release="1.0.0", remove="2.0.0")
def density_from_PDB(pdb, **kwargs):
    """Create a density from a single frame PDB.

    Typical use is to make a density from the crystal water
    molecules. The density is created from isotropic gaussians
    centered at each selected atoms. If B-factors are present in the
    file then they are used to calculate the width of the gaussian.

    Using the *sigma* keyword, one can override this choice and
    prescribe a gaussian width for all atoms (in Angstrom), which is
    calculated as described for :func:`Bfactor2RMSF`.


    Parameters
    ----------
    pdb : str
          PDB filename (should have the temperatureFactor set); ANISO
          records are currently *not* processed
    select : str
          selection string (MDAnalysis syntax) for the species to be analyzed
          ['resname HOH and name O']
    delta : float
          bin size for the density grid in Angstrom (same in x,y,z) [1.0]
    metadata : dict
          dictionary of additional data to be saved with the object [``None``]
    padding : float
          increase histogram dimensions by padding (on top of initial box size) [1.0]
    sigma : float
          width (in Angstrom) of the gaussians that are used to build up the
          density; if ``None`` then uses B-factors from *pdb* [``None``]

    Returns
    -------
    :class:`Density`
          object with a density measured relative to the water density at
          standard conditions

    Notes
    -----
    The current implementation is *painfully* slow.

    See Also
    --------
    :func:`Bfactor2RMSF` and :class:`BfactorDensityCreator`


    .. deprecated: 1.0.0
       This function is not well tested or optimized and will be removed
       in 2.0.0.
    """
    return BfactorDensityCreator(pdb, **kwargs).Density()


# REMOVE in 2.0.0
class BfactorDensityCreator(object):
    """Create a density grid from a pdb file using MDAnalysis.

    The main purpose of this function is to convert crystal waters in
    an X-ray structure into a density so that one can compare the
    experimental density with the one from molecular dynamics
    trajectories. Because a pdb is a single snapshot, the density is
    estimated by placing Gaussians of width sigma at the position of
    all selected atoms.

    Sigma can be fixed or taken from the B-factor field, in which case
    sigma is taken as sqrt(3.*B/8.)/pi (see :func:`BFactor2RMSF`).

    .. TODO
    .. * Make Gaussian convolution more efficient (at least for same
    ..   sigma) because right now it is *very* slow (which may be
    ..   acceptable if one only runs this once)
    .. * Using a temporary Creator class with the
    ..   :meth:`BfactorDensityCreator.Density` helper method is clumsy.


    .. deprecated: 1.0.0
       This function is not well tested or optimized and will be removed
       in 2.0.0.
    """

    @deprecate(release="1.0.0", remove="2.0.0")
    def __init__(self, pdb, delta=1.0, select='resname HOH and name O',
                 metadata=None, padding=1.0, sigma=None):
        """Construct the density from psf and pdb and the select.

        Parameters
        ----------
        pdb : str
            PDB file or :class:`MDAnalysis.Universe`;
        select : str
            selection string (MDAnalysis syntax) for the species to be analyzed
        delta : float
            bin size for the density grid in Angstrom (same in x,y,z) [1.0]
        metadata : dict
            dictionary of additional data to be saved with the object
        padding : float
            increase histogram dimensions by padding (on top of initial box size)
        sigma : float
            width (in Angstrom) of the gaussians that are used to build up the
            density; if ``None`` (the default) then uses B-factors from pdb


        Notes
        -----
        For assigning X-ray waters to MD densities one might have to use a sigma
        of about 0.5 A to obtain a well-defined and resolved x-ray water density
        that can be easily matched to a broader density distribution.

        .. versionchanged:: 1.0.0
           Changed `selection` keyword to `select`

        Examples
        --------
        The following creates the density with the B-factors from the pdb file::

          DC = BfactorDensityCreator(pdb, delta=1.0, select="name HOH",
                                     padding=2, sigma=None)
          density = DC.Density()

        See Also
        --------
        :func:`density_from_PDB` for a convenience function

        """
        u = MDAnalysis.as_Universe(pdb)
        group = u.select_atoms(select)
        coord = group.positions
        logger.info("Selected {0:d} atoms ({1!s}) out of {2:d} total.".format(coord.shape[0], select, len(u.atoms)))
        smin = np.min(coord, axis=0) - padding
        smax = np.max(coord, axis=0) + padding

        BINS = fixedwidth_bins(delta, smin, smax)
        arange = list(zip(BINS['min'], BINS['max']))
        bins = BINS['Nbins']

        # get edges by doing a fake run
        grid, self.edges = np.histogramdd(np.zeros((1, 3)),
                                             bins=bins, range=arange, normed=False)
        self.delta = np.diag([(e[-1] - e[0]) / (len(e) - 1) for e in self.edges])
        self.midpoints = [0.5 * (e[:-1] + e[1:]) for e in self.edges]
        self.origin = [m[0] for m in self.midpoints]
        n_frames = 1

        if sigma is None:
            # histogram individually, and smear out at the same time
            # with the appropriate B-factor
            if np.any(group.tempfactors == 0.0):
                wmsg = "Some B-factors are Zero (will be skipped)."
                logger.warning(wmsg)
                warnings.warn(wmsg, category=MissingDataWarning)
            rmsf = Bfactor2RMSF(group.tempfactors)
            grid *= 0.0  # reset grid
            self.g = self._smear_rmsf(coord, grid, self.edges, rmsf)
        else:
            # histogram 'delta functions'
            grid, self.edges = np.histogramdd(coord, bins=bins, range=arange, normed=False)
            logger.info("Histogrammed {0:6d} atoms from pdb.".format(len(group.atoms)))
            # just a convolution of the density with a Gaussian
            self.g = self._smear_sigma(grid, sigma)

        try:
            metadata['pdb'] = pdb
        except TypeError:
            metadata = {'pdb': pdb}
        metadata['select'] = select
        metadata['n_frames'] = n_frames
        metadata['sigma'] = sigma
        self.metadata = metadata

        logger.info("Histogram completed (initial density in Angstrom**-3)")

        # Density automatically converts histogram to density for isDensity=False -- ??[OB]

    def Density(self, threshold=None):
        """Returns a :class:`Density` object."""
        d = Density(grid=self.g, edges=self.edges, units=dict(length='Angstrom'),
                    parameters=dict(isDensity=False), metadata=self.metadata)
        d.make_density()
        d.convert_density('water')
        return d

    def _smear_sigma(self, grid, sigma):
        # smear out points
        # (not optimized -- just to test the principle; faster approach could use
        # convolution of the whole density with a single Gaussian via FFTs:
        # rho_smeared = F^-1[ F[g]*F[rho] ]
        g = np.zeros(grid.shape)  # holds the smeared out density
        pos = np.where(grid != 0)  # position in histogram (as bin numbers)
        for iwat in range(len(pos[0])):  # super-ugly loop
            p = tuple([wp[iwat] for wp in pos])
            g += grid[p] * np.fromfunction(self._gaussian, grid.shape, dtype=np.int, p=p, sigma=sigma)
            # print("Smearing out atom position {0:4d}/{1:5d} with RMSF {2:4.2f} A\r".format(iwat + 1, len(pos[0]), sigma),)
        return g

    def _smear_rmsf(self, coordinates, grid, edges, rmsf):
        # smear out each water with its individual Gaussian
        # (slower than smear_sigma)
        g = np.zeros(grid.shape)  # holds the smeared out density
        N, D = coordinates.shape
        for iwat, coord in enumerate(coordinates):
            if rmsf[iwat] == 0:
                continue
            g += np.fromfunction(self._gaussian_cartesian, grid.shape, dtype=np.int,
                                    c=coord, sigma=rmsf[iwat])
            # print("Smearing out atom position {0:4d}/{1:5d} with RMSF {2:4.2f} A\r".format(iwat + 1, N, rmsf[iwat]),)
        return g

    def _gaussian(self, i, j, k, p, sigma):
        # i,j,k can be numpy arrays
        # p is center of gaussian as grid index, sigma its width (in A)
        x = self.delta[0, 0] * (i - p[0])  # in Angstrom
        y = self.delta[1, 1] * (j - p[1])
        z = self.delta[2, 2] * (k - p[2])
        return (2 * np.pi * sigma) ** (-1.5) * np.exp(-(x * x + y * y + z * z) / (2 * sigma * sigma))

    def _gaussian_cartesian(self, i, j, k, c, sigma):
        # i,j,k can be numpy arrays
        # c is center of gaussian in cartesian coord (A), sigma its width (in A)
        x = self.origin[0] + self.delta[0, 0] * i - c[0]  # in Angstrom
        y = self.origin[1] + self.delta[1, 1] * j - c[1]
        z = self.origin[2] + self.delta[2, 2] * k - c[2]
        return (2 * np.pi * sigma) ** (-1.5) * np.exp(-(x * x + y * y + z * z) / (2 * sigma * sigma))
