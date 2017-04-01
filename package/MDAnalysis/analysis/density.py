# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

# MDAnalysis -- density analysis
# Copyright (c) 2007-2011 Oliver Beckstein <orbeckst@gmail.com>
# (based on code from Hop --- a framework to analyze solvation dynamics from MD simulations)

r"""
Generating densities from trajectories --- :mod:`MDAnalysis.analysis.density`
=============================================================================

:Author: Oliver Beckstein
:Year: 2011
:Copyright: GNU Public License v3

The module provides classes and functions to generate and represent
volumetric data, in particular densities.

Generating a density from a MD trajectory
-----------------------------------------

An input trajectory must

1. have been centered on the protein of interest;
2. have all molecules made whole that have been broken across periodic
   boundaries;
3. have the solvent molecules remapped so that they are closest to the
   solute (this is important when using funky unit cells such as
   a dodecahedron or a truncated octahedron).

To generate the density of water molecules around a protein::

  from MDAnalysis.analysis.density import density_from_Universe
  u = Universe(TPR, XTC)
  D = density_from_Universe(u, delta=1.0, atomselection="name OW")
  D.convert_density('TIP4P')
  D.export("water.dx")

The positions of all water oxygens are histogrammed on a grid with spacing
*delta* = 1 Å. Initially the density is measured in :math:`\text{Å}^{-3}`. With the
:meth:`Density.convert_density` method, the units of measurement are
changed. In the example we are now measuring the density relative to the
literature value of the TIP4P water model at ambient conditions (see the values
in :data:`MDAnalysis.units.water` for details). Finally, the density is
writte as an OpenDX_ compatible file that can be read in VMD_ or PyMOL_.

See :class:`Density` for details. In particular, the density is stored
as a NumPy array in :attr:`Density.grid`, which can be processed in
any manner.


Creating densities
------------------

The following functions take trajectory or coordinate data and generate a
:class:`Density` object.

.. autofunction:: density_from_Universe
.. autofunction:: density_from_PDB
.. autofunction:: Bfactor2RMSF


Supporting classes and functions
--------------------------------

The main output of the density creation functions is a
:class:`Density` instance, which is derived from a
:class:`gridData.core.Grid`. A :class:`Density` is essentially, a 3D
array with origin and lengths together with associated metadata (which
can be used in downstream processing).

.. autoclass:: Density
   :members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: BfactorDensityCreator
   :members:

.. autofunction:: notwithin_coordinates_factory


.. Links
.. -----

.. _OpenDX: http://www.opendx.org/
.. _VMD:   http://www.ks.uiuc.edu/Research/vmd/
.. _PyMOL: http://www.pymol.org/
.. _GridDataFormats: https://github.com/MDAnalysis/GridDataFormats

"""

from __future__ import print_function
from six.moves import range

import numpy as np
import sys
import os
import os.path
import errno
import warnings

try:
    from gridData import Grid
except ImportError:
    raise ImportError(
        """ImportError: The GridDataFormats package can not be found!

        The 'gridData' module from GridDataFormats could not be
        imported. Please install it first.  You can try installing
        directly from the internet:

          pip install GridDataFormats

        or

          conda config --add channels conda-forge
          conda install griddataformats

        Alternatively, download the package from

          http://pypi.python.org/pypi/GridDataFormats/

        and install in the usual manner.
        """
    )

import MDAnalysis
from MDAnalysis.core import groups
from MDAnalysis.lib.util import fixedwidth_bins, iterable, asiterable
from MDAnalysis.lib import NeighborSearch as NS
from MDAnalysis import NoDataError, MissingDataWarning
from .. import units
from ..lib import distances
from MDAnalysis.lib.log import ProgressMeter

import logging

logger = logging.getLogger("MDAnalysis.analysis.density")


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

        (Actually, the default unit is the value of
        ``MDAnalysis.core.flags['length_unit']``; in most
        cases this is "Angstrom".)
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

    The :meth:`Density.export` method always exports a 3D object
    (written in such a way to be readable in VMD_ and PyMOL_), the
    rest should work for an array of any dimension.

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

    See Also
    --------
    gridData.core.Grid : the base class of :class:`Density`.


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
       Ångstrom, with the density values itself measured relative to the
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
        length_unit = MDAnalysis.core.flags['length_unit']

        parameters = kwargs.pop('parameters', {})
        if (len(args) > 0 and type(args[0]) is str) or (type(kwargs.get('grid', None) is str)):
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
                    raise ValueError('Unit ' + str(value) + ' of type ' + str(unit_type) + ' is not recognized.')
        except AttributeError:
            errmsg = '"unit" must be a dictionary with keys "length" and "density.'
            logger.fatal(errmsg)
            raise ValueError(errmsg)
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
            logger.warn(msg)
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
            raise ValueError("The name of the unit ({0!r} supplied) must be one of:\n{1!r}".format(unit, units.conversion_factor['density'].keys()))
        self.units['density'] = unit

    def __repr__(self):
        if self.parameters['isDensity']:
            grid_type = 'density'
        else:
            grid_type = 'histogram'
        return '<Density ' + grid_type + ' with ' + str(self.grid.shape) + ' bins>'

def density_from_Universe(universe, delta=1.0, atomselection='name OH2',
                          start=None, stop=None, step=None,
                          metadata=None, padding=2.0, cutoff=0, soluteselection=None,
                          use_kdtree=True, update_selection=False,
                          verbose=None, interval=1, quiet=None,
                          parameters=None):
    """Create a density grid from a :class:`MDAnalysis.Universe` object.

    The trajectory is read, frame by frame, and the atoms selected with `atomselection` are
    histogrammed on a grid with spacing `delta`.

    Parameters
    ----------
    universe : MDAnalysis.Universe
            :class:`MDAnalysis.Universe` object with a trajectory
    atomselection : str (optional)
            selection string (MDAnalysis syntax) for the species to be analyzed
            ["name OH2"]
    delta : float (optional)
            bin size for the density grid in Angstroem (same in x,y,z) [1.0]
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
            in Angstroem [2.0]
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

    Returns
    -------
    :class:`Density`
            A :class:`Density` instance with the histogrammed data together
            with associated metadata.


    Notes
    -----

    By default, the `atomselection` is static, i.e., atoms are only selected
    once at the beginning. If you want *dynamically changing selections* (such
    as "name OW and around 4.0 (protein and not name H*)", i.e., the water
    oxygen atoms that are within 4 Å of the protein heavy atoms) then set
    ``update_selection=True``. For the special case of calculating a density of
    the "bulk" solvent away from a solute use the optimized selections with
    keywords *cutoff* and *soluteselection* (see Examples below).

    Examples
    --------
    Basic use for creating a water density (just using the water oxygen atoms "OW")::

      density = density_from_Universe(universe, delta=1.0, atomselection='name OW')

    If you are only interested in water within a certain region, e.g., within a
    vicinity around a binding site, you can use a selection that updates every
    step by setting the `update_selection` keyword argument::

      site_density = density_from_Universe(universe, delta=1.0,
                                           atomselection='name OW and around 5 (resid 156 157 305)',
                                           update_selection=True)

    A special case for an updating selection is to create the "bulk density",
    i.e., the water outside the immediate solvation shell of a protein: Select
    all water oxygen atoms that are *farther away* than a given cut-off (say, 4
    Å) from the solute (here, heavy atoms of the protein)::

      bulk = density_from_Universe(universe, delta=1.0, atomselection='name OW',
                                   solute="protein and not name H*",
                                   cutoff=4)

    (Using the special case for the bulk with `soluteselection` and `cutoff`
    improves performance over the simple `update_selection` approach.)

    .. versionchanged:: 0.13.0
       *update_selection* and *quiet* keywords added

    .. deprecated:: 0.16
       The keyword argument *quiet* is deprecated in favor of *verbose*.

    """
    u = universe

    if cutoff > 0 and soluteselection is not None:
        # special fast selection for '<atomsel> not within <cutoff> of <solutesel>'
        notwithin_coordinates = notwithin_coordinates_factory(
            u, atomselection, soluteselection, cutoff,
            use_kdtree=use_kdtree, updating_selection=update_selection)
        def current_coordinates():
            return notwithin_coordinates()
    else:
        group = u.select_atoms(atomselection, updating=update_selection)

        def current_coordinates():
            return group.positions

    coord = current_coordinates()
    logger.info(
        "Selected {0:d} atoms out of {1:d} atoms ({2!s}) from {3:d} total."
        "".format(coord.shape[0], len(u.select_atoms(atomselection)),
                  atomselection, len(u.atoms))
    )

    # mild warning; typically this is run on RMS-fitted trajectories and
    # so the box information is rather meaningless
    box, angles = u.trajectory.ts.dimensions[:3], u.trajectory.ts.dimensions[3:]
    if tuple(angles) != (90., 90., 90.):
        msg = "Non-orthorhombic unit-cell --- make sure that it has been remapped properly!"
        warnings.warn(msg)
        logger.warning(msg)

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

    pm = ProgressMeter(u.trajectory.n_frames, interval=interval,
                       verbose=verbose, quiet=quiet,
                       format="Histogramming %(n_atoms)6d atoms in frame "
                       "%(step)5d/%(numsteps)d  [%(percentage)5.1f%%]\r")
    start, stop, step = u.trajectory.check_slice_indices(start, stop, step)
    for ts in u.trajectory[start:stop:step]:
        coord = current_coordinates()

        pm.echo(ts.frame, n_atoms=len(coord))
        if len(coord) == 0:
            continue

        h[:], edges[:] = np.histogramdd(coord, bins=bins, range=arange, normed=False)
        grid += h  # accumulate average histogram

    n_frames = len(range(start, stop, step))
    grid /= float(n_frames)

    metadata = metadata if metadata is not None else {}
    metadata['psf'] = u.filename
    metadata['dcd'] = u.trajectory.filename
    metadata['atomselection'] = atomselection
    metadata['n_frames'] = n_frames
    metadata['totaltime'] = round(u.trajectory.n_frames * u.trajectory.dt, 3)
    metadata['dt'] = u.trajectory.dt
    metadata['time_unit'] = MDAnalysis.core.flags['time_unit']
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


    g = Density(grid=grid, edges=edges, units={'length': MDAnalysis.core.flags['length_unit']},
                parameters=parameters, metadata=metadata)
    g.make_density()
    logger.info("Density completed (initial density in Angstrom**-3)")

    return g


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


def Bfactor2RMSF(B):
    r"""Atomic root mean square fluctuation (in Angstrom) from the crystallographic B-factor

    RMSF and B-factor are related by [Willis1975]_

    .. math::

        B = \frac{8\pi^2}{3} \rm{RMSF}^2

    and this function returns

    .. math::

        \rm{RMSF} = \sqrt{\frac{3 B}{8\pi^2}}

    .. rubric:: References

    .. [Willis1975]  BTM Willis and AW Pryor. *Thermal vibrations in crystallography*. Cambridge Univ. Press, 1975
    """
    return np.sqrt(3. * B / 8.) / np.pi


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
    atomselection : str
          selection string (MDAnalysis syntax) for the species to be analyzed
          ['resname HOH and name O']
    delta : float
          bin size for the density grid in Angstroem (same in x,y,z) [1.0]
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

    """
    return BfactorDensityCreator(pdb, **kwargs).Density()


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

    """

    def __init__(self, pdb, delta=1.0, atomselection='resname HOH and name O',
                 metadata=None, padding=1.0, sigma=None):
        """Construct the density from psf and pdb and the atomselection.

        Parameters
        ----------
        pdb : str
            PDB file or :class:`MDAnalysis.Universe`;
        atomselection : str
            selection string (MDAnalysis syntax) for the species to be analyzed
        delta : float
            bin size for the density grid in Angstroem (same in x,y,z) [1.0]
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

        Examples
        --------
        The following creates the density with the B-factors from the pdb file::

          DC = BfactorDensityCreator(pdb, delta=1.0, atomselection="name HOH",
                                     padding=2, sigma=None)
          density = DC.Density()

        See Also
        --------
        :func:`density_from_PDB` for a convenience function

        """
        u = MDAnalysis.as_Universe(pdb)
        group = u.select_atoms(atomselection)
        coord = group.positions
        logger.info("Selected {0:d} atoms ({1!s}) out of {2:d} total.".format(coord.shape[0], atomselection, len(u.atoms)))
        smin = np.min(coord, axis=0) - padding
        smax = np.max(coord, axis=0) + padding

        BINS = fixedwidth_bins(delta, smin, smax)
        arange = zip(BINS['min'], BINS['max'])
        bins = BINS['Nbins']

        # get edges by doing a fake run
        grid, self.edges = np.histogramdd(np.zeros((1, 3)),
                                             bins=bins, range=arange, normed=False)
        self.delta = np.diag(map(lambda e: (e[-1] - e[0]) / (len(e) - 1), self.edges))
        self.midpoints = map(lambda e: 0.5 * (e[:-1] + e[1:]), self.edges)
        self.origin = map(lambda m: m[0], self.midpoints)
        n_frames = 1

        if sigma is None:
            # histogram individually, and smear out at the same time
            # with the appropriate B-factor
            if np.any(group.bfactors == 0.0):
                wmsg = "Some B-factors are Zero (will be skipped)."
                logger.warn(wmsg)
                warnings.warn(wmsg, category=MissingDataWarning)
            rmsf = Bfactor2RMSF(group.bfactors)
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
        metadata['atomselection'] = atomselection
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
            print("Smearing out atom position {0:4d}/{1:5d} with RMSF {2:4.2f} A\r".format(iwat + 1, len(pos[0]), sigma),)
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
            print("Smearing out atom position {0:4d}/{1:5d} with RMSF {2:4.2f} A\r".format(iwat + 1, N, rmsf[iwat]),)
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
