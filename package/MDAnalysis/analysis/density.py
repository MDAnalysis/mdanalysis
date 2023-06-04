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


.. versionchanged:: 2.0.0
   Deprecated :func:`density_from_Universe`, :func:`density_from_PDB`, and
   :func:`Bfactor2RMSF` have now been removed.


Generating a density from a MD trajectory
-----------------------------------------

A common use case is to analyze the solvent density around a protein of
interest. The density is calculated with :class:`DensityAnalysis` in the
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
  D.results.density.convert_density('TIP4P')
  D.results.density.export("water.dx", type="double")

The positions of all water oxygens (the :class:`AtomGroup` `ow`) are
histogrammed on a grid with spacing *delta* = 1 Å. Initially the density is
measured in :math:`\text{Å}^{-3}`. With the :meth:`Density.convert_density`
method, the units of measurement are changed. In the example we are now
measuring the density relative to the literature value of the TIP4P water model
at ambient conditions (see the values in :data:`MDAnalysis.units.water` for
details). Finally, the density is written as an OpenDX_ compatible file that
can be read in VMD_, Chimera_, or PyMOL_.

The :class:`Density` object is accessible as the
:attr:`DensityAnalysis.results.density` attribute.  In particular, the data
for the density is stored as a NumPy array in :attr:`Density.grid`, which can
be processed in any manner.


Creating densities
------------------

The :class:`DensityAnalysis` class generates a :class:`Density` from an
atomgroup.

.. autoclass:: DensityAnalysis
   :members:
   :inherited-members: run

   .. automethod:: _set_user_grid


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
import numpy as np
import sys
import os
import os.path
import errno
import warnings

from gridData import Grid

import MDAnalysis
from MDAnalysis.core import groups
from MDAnalysis.lib.util import (fixedwidth_bins, iterable, asiterable,
                                 deprecate,)
from MDAnalysis.lib import NeighborSearch as NS
from MDAnalysis import NoDataError, MissingDataWarning
from .. import units
from ..lib import distances
from MDAnalysis.lib.log import ProgressBar

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
            User defined x dimension box edge in ångström.
    ydim : float (optional)
            User defined y dimension box edge in ångström.
    zdim : float (optional)
            User defined z dimension box edge in ångström.

    Attributes
    ----------
    results.density : :class:`Density`
            A :class:`Density` instance containing a physical density of units
            :math:`Angstrom^{-3}`.

            After the analysis (see the :meth:`~DensityAnalysis.run` method),
            the resulting density is stored in the :attr:`results.density`
            attribute as a :class:`Density` instance. Note: this replaces the
            now deprecated :attr:`density` attribute.

    density : :class:`Density`
            Alias to the :attr:`results.density`.

            .. deprecated:: 2.0.0
               Will be removed in MDAnalysis 3.0.0. Please use
               :attr:`results.density` instead.

    Raises
    ------
    ValueError
        if AtomGroup is empty and no user defined grid is provided, or
        if the user defined grid is not or incorrectly provided
    UserWarning
        if AtomGroup is empty and a user defined grid is provided


    See Also
    --------
    pmda.density.DensityAnalysis
        A parallel version of :class:`DensityAnalysis`

    Notes
    -----
    If the `gridcenter` and `x/y/zdim` arguments are not provided,
    :class:`DensityAnalysis` will attempt to automatically generate
    a gridbox from the atoms in 'atomgroup' (See Examples).

    Normal :class:`AtomGroup` instances represent a static selection of
    atoms. If you want *dynamically changing selections* (such as "name OW and
    around 4.0 (protein and not name H*)", i.e., the water oxygen atoms that
    are within 4 Å of the protein heavy atoms) then create an
    :class:`~MDAnalysis.core.groups.UpdatingAtomGroup` (see Examples).

    :class:`DensityAnalysis` will fail when the :class:`AtomGroup` instance
    does not contain any selection of atoms, even when `updating` is set to
    ``True``. In such a situation, user defined box limits can be provided to
    generate a `Density`. Although, it remains the user's responsibility
    to ensure that the provided grid limits encompass atoms to be selected
    on all trajectory frames.

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
        D.results.density.convert_density('TIP4P')

    The positions of all water oxygens are histogrammed on a grid with spacing
    *delta* = 1 Å and stored as a :class:`Density` object in the attribute
    :attr:`DensityAnalysis.results.density`.

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

      probability_density = D.results.density.grid / D.results.density.grid.sum()

    Similarly, if you would like to recover a grid containing a **histogram of
    atom counts**, simply multiply by the volume `dV` of each bin (or voxel);
    in this case you need to ensure that the physical density is measured in
    Å\ :sup:`-3` by converting it::

      import numpy as np

      # ensure that the density is A^{-3}
      D.results.density.convert_density("A^{-3}")

      dV = np.prod(D.results.density.delta)
      atom_count_histogram = D.results.density.grid * dV


    .. rubric:: Writing the density to a file

    A density can be `exported to different formats
    <https://www.mdanalysis.org/GridDataFormats/gridData/formats.html>`_ with
    :meth:`Density.export` (thanks to the fact that :class:`Density` is a
    subclass :class:`gridData.core.Grid`, which provides the functionality).
    For example, to `write a DX file
    <https://www.mdanalysis.org/GridDataFormats/gridData/basic.html#writing-out-data>`_
    ``water.dx`` that can be read with VMD, PyMOL, or Chimera::

      D.results.density.export("water.dx", type="double")


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
    .. versionchanged:: 2.0.0
       :func:`_set_user_grid` is now a method of :class:`DensityAnalysis`.
       :class:`Density` results are now stored in a
       :class:`MDAnalysis.analysis.base.Results` instance.

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
        if (self._gridcenter is not None or
                any([self._xdim, self._ydim, self._zdim])):
            # Issue 2372: padding is ignored, defaults to 2.0 therefore warn
            if self._padding > 0:
                msg = (f"Box padding (currently set at {self._padding}) "
                       f"is not used in user defined grids.")
                warnings.warn(msg)
                logger.warning(msg)
            # Generate a copy of smin/smax from coords to later check if the
            # defined box might be too small for the selection
            try:
                smin = np.min(coord, axis=0)
                smax = np.max(coord, axis=0)
            except ValueError as err:
                msg = ("No atoms in AtomGroup at input time frame. "
                       "This may be intended; please ensure that "
                       "your grid selection covers the atomic "
                       "positions you wish to capture.")
                warnings.warn(msg)
                logger.warning(msg)
                smin = self._gridcenter     #assigns limits to be later -
                smax = self._gridcenter     #overwritten by _set_user_grid
            # Overwrite smin/smax with user defined values
            smin, smax = self._set_user_grid(self._gridcenter, self._xdim,
                                             self._ydim, self._zdim, smin,
                                             smax)
        else:
            # Make the box bigger to avoid as much as possible 'outlier'. This
            # is important if the sites are defined at a high density: in this
            # case the bulk regions don't have to be close to 1 * n0 but can
            # be less. It's much more difficult to deal with outliers.  The
            # ideal solution would use images: implement 'looking across the
            # periodic boundaries' but that gets complicated when the box
            # rotates due to RMS fitting.
            try:
                smin = np.min(coord, axis=0) - self._padding
                smax = np.max(coord, axis=0) + self._padding
            except ValueError as err:
                errmsg = ("No atoms in AtomGroup at input time frame. "
                          "Grid for density could not be automatically"
                          " generated. If this is expected, a user"
                          " defined grid will need to be "
                          "provided instead.")
                raise ValueError(errmsg) from err
        BINS = fixedwidth_bins(self._delta, smin, smax)
        arange = np.transpose(np.vstack((BINS['min'], BINS['max'])))
        bins = BINS['Nbins']
        # create empty grid with the right dimensions (and get the edges)
        grid, edges = np.histogramdd(np.zeros((1, 3)), bins=bins,
                                     range=arange, density=False)
        grid *= 0.0
        self._grid = grid
        self._edges = edges
        self._arange = arange
        self._bins = bins

    def _single_frame(self):
        h, _ = np.histogramdd(self._atomgroup.positions,
                              bins=self._bins, range=self._arange,
                              density=False)
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
        self.results.density = density

    @property
    def density(self):
        wmsg = ("The `density` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.density` instead")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.density

    @staticmethod
    def _set_user_grid(gridcenter, xdim, ydim, zdim, smin, smax):
        """Helper function to set the grid dimensions to user defined values

        Parameters
        ----------
        gridcenter : numpy ndarray, float32
                3 element ndarray containing the x, y and z coordinates of the
                grid box center
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


        .. versionchanged:: 2.0.0
           Now a staticmethod of :class:`DensityAnalysis`.
        """
        # Check user inputs
        if any(x is None for x in [gridcenter, xdim, ydim, zdim]):
            errmsg = ("Gridcenter or grid dimensions are not provided")
            raise ValueError(errmsg)
        try:
            gridcenter = np.asarray(gridcenter, dtype=np.float32).reshape(3,)
        except ValueError as err:
            raise ValueError("Gridcenter must be a 3D coordinate") from err
        try:
            xyzdim = np.array([xdim, ydim, zdim], dtype=np.float32)
        except ValueError as err:
            raise ValueError("xdim, ydim, and zdim must be numbers") from err
        if any(np.isnan(gridcenter)) or any(np.isnan(xyzdim)):
            raise ValueError("Gridcenter or grid dimensions have NaN element")


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
    a probability density, divide it by :meth:`Density.grid.sum` or use ``density=True``
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
        if (len(args) > 0 and isinstance(args[0], str) or
            isinstance(kwargs.get('grid', None), str)):
            # try to be smart: when reading from a file then it is likely that
            # this is a density
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
                    errmsg = (f"Unit {value} of type {unit_type} is not "
                              f"recognized.")
                    raise ValueError(errmsg) from None
        except AttributeError:
            errmsg = '"unit" must be a dictionary with keys "length" and "density.'
            logger.fatal(errmsg)
            raise ValueError(errmsg) from None
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
            errmsg = 'The grid is not a density so convert_density() makes no sense.'
            logger.fatal(errmsg)
            raise RuntimeError(errmsg)
        if unit == self.units['density']:
            return
        try:
            self.grid *= units.get_conversion_factor('density',
                                                     self.units['density'], unit)
        except KeyError:
            errmsg = (f"The name of the unit ({unit} supplied) must be one "
                      f"of:\n{units.conversion_factor['density'].keys()}")
            raise ValueError(errmsg) from None
        self.units['density'] = unit

    def __repr__(self):
        if self.parameters['isDensity']:
            grid_type = 'density'
        else:
            grid_type = 'histogram'
        return '<Density ' + grid_type + ' with ' + str(self.grid.shape) + ' bins>'
