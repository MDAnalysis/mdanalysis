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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

# Water Bridge Analysis
r"""Water Bridge analysis --- :mod:`MDAnalysis.analysis.hbonds.WaterBridgeAnalysis`
===============================================================================

:Author: Zhiyi Wu
:Year: 2017
:Copyright: GNU Public License v3
:Maintainer: Zhiyi Wu <zhiyi.wu@gtc.ox.ac.uk>,  `@xiki-tempula`_ on GitHub


.. _`@xiki-tempula`: https://github.com/xiki-tempula


Given a :class:`~MDAnalysis.core.universe.Universe` (simulation
trajectory with 1 or more frames) measure all water bridges for each
frame between selections 1 and 2.
Water bridge is defined as a bridging water which simultaneously forms
two hydrogen bonds with atoms from both selection 1 and selection 2.

A water bridge can form between two hydrogen bond acceptors.

e.g. -CO\ :sub:`2`\ :sup:`-`:···H−O−H···:\ :sup:`-`\ O\ :sub:`2`\ C-

A water bridge can also form between two hydrogen bond donors.

e.g. -NH···:O:···HN- (where O is the oxygen of a bridging water)

A hydrogen bond acceptor and another hydrogen bond donor can be bridged by a
water.

e.g. -CO\ :sub:`2`\ :sup:`-`:···H−O:···HN- (where H−O is part of **H−O**\ −H)

The :class:`WaterBridgeAnalysis` class is modeled after the  \
:class:`~MDAnalysis.analysis.hbonds.hbond_analysis.HydrogenBondAnalysis`.

The following keyword arguments are important to control the behavior of the
water bridge analysis:

 - *water_selection* (``resname SOL``): the selection string for the bridging
   water
 - donor-acceptor *distance* (Å): 3.0
 - Angle *cutoff* (degrees): 120.0
 - *forcefield* to switch between default values for different force fields
 - *donors* and *acceptors* atom types (to add additional atom names)

.. _wb_Analysis_Output:

Output
------

The results are a list of hydrogen bonds between the selection 1 or selection 2
and the bridging water.

Each list is formated similar to the \ :attr:`HydrogenBondAnalysis.timeseries
<MDAnalysis.analysis.hbonds.hbond_analysis.HydrogenBondAnalysis.timeseries>`
and contains

  - the **identities** of donor and acceptor heavy-atoms,
  - the **distance** between the heavy atom acceptor atom and the hydrogen atom
  - the **angle** donor-hydrogen-acceptor angle (180º is linear).

Water bridge data are returned per frame, which is stored in \
:attr:`WaterBridgeAnalysis.timeseries` (In the following description, ``#``
indicates comments that are not part of the output.)::

    results = [
        [ # frame 1
           # hbonds linking the selection 1 and selection 2 to the bridging
           # water 1
           [ # hbond 1 from selection 1 to the bridging water 1
              <donor index (0-based)>,
              <acceptor index (0-based)>, <donor string>, <acceptor string>,
              <distance>, <angle>
           ],
           [ # hbond 2 from selection 1 to the bridging water 1
              <donor index (0-based)>,
              <acceptor index (0-based)>, <donor string>, <acceptor string>,
              <distance>, <angle>
           ],
           [ # hbond 1 from selection 2 to the bridging water 1
              <donor index (0-based)>,
              <acceptor index (0-based)>, <donor string>, <acceptor string>,
              <distance>, <angle>
           ],
           [ # hbond 2 from selection 2 to the bridging water 1
              <donor index (0-based)>,
              <acceptor index (0-based)>, <donor string>, <acceptor string>,
              <distance>, <angle>
           ],

           # hbonds linking the selection 1 and selection 2 to the bridging
           # water 2
           [ # hbond 1 from selection 1 to the bridging water 2
              <donor index (0-based)>,
              <acceptor index (0-based)>, <donor string>, <acceptor string>,
              <distance>, <angle>
           ],
           [ # hbond 1 from selection 2 to the bridging water 2
              <donor index (0-based)>,
              <acceptor index (0-based)>, <donor string>, <acceptor string>,
              <distance>, <angle>
           ],
           ....
        ],
        [ # frame 2
          [ ... ], [ ... ], ...
        ],
        ...
    ]

Using the :meth:`WaterBridgeAnalysis.generate_table` method one can reformat
the results as a flat "normalised" table that is easier to import into a
database or dataframe for further processing.
:meth:`WaterBridgeAnalysis.save_table` saves the table to a pickled file. The
table itself is a :class:`numpy.recarray`.

Detection of water bridges
---------------------------
Water bridges are recorded if a bridging water simultaneously forms two
hydrogen bonds with selection 1 and selection 2.

Hydrogen bonds are detected as is described in \
:class:`~MDAnalysis.analysis.hbonds.hbond_analysis.HydrogenBondAnalysis`, see \
:ref:`Detection-of-hydrogen-bonds`.

The lists of donor and acceptor names can be extended by providing lists of
atom names in the `donors` and `acceptors` keywords to
:class:`WaterBridgeAnalysis`. If the lists are entirely inappropriate
(e.g. when analysing simulations done with a force field that uses very
different atom names) then one should either use the value "other" for
`forcefield` to set no default values, or derive a new class and set the
default list oneself::

 class WaterBridgeAnalysis_OtherFF(WaterBridgeAnalysis):
       DEFAULT_DONORS = {"OtherFF": tuple(set([...]))}
       DEFAULT_ACCEPTORS = {"OtherFF": tuple(set([...]))}

Then simply use the new class instead of the parent class and call it with
`forcefield` = "OtherFF". Please also consider contributing the list of heavy
atom names to MDAnalysis.

How to perform WaterBridgeAnalysis
-----------------------------------

All water bridges between arginine and aspartic acid can be analysed with ::

  import MDAnalysis
  import MDAnalysis.analysis.hbonds

  u = MDAnalysis.Universe('topology', 'trajectory')
  w = MDAnalysis.analysis.hbonds.WaterBridgeAnalysis(u, 'resname ARG', 'resname ASP')
  w.run()

The results are stored as the attribute
:attr:`WaterBridgeAnalysis.timeseries`; see :ref:`wb_Analysis_Output` for the
format.

An example of using the :attr:`~WaterBridgeAnalysis.timeseries` would be
detecting the percentage of time a certain water bridge exits.

Trajectory :code:`u` has two frames, where the first frame contains a water
bridge from the oxygen of the first arginine to the oxygen of the third
aspartate. No water bridge is detected in the second frame. ::

  print(w.timeseries)

prints out (the comments are not part of the data structure but are added here
for clarity): ::

  [ # frame 1
    # A water bridge SOL2 links O from ARG1 and ASP3
   [[0,1,'ARG1:O',  'SOL2:HW1',3.0,180],
    [2,3,'SOL2:HW2','ASP3:O',  3.0,180],
   ],
    # frame 2
    # No water bridge detected
   []
  ]

To calculate the percentage, we can iterate through :code:`w.timeseries`. ::

  water_bridge_presence = []
  for frame in w.timeseries:
      if frame:
          water_bridge_presence.append(True)
      else:
          water_bridge_presence.append(False)
  p_bridge = float(sum(water_bridge_presence))/len(water_bridge_presence)
  print("Fraction of time with water bridge present: {}".format(p_bridge))

In the example above, :code:`p_bridge` would become 0.5, i.e., for 50% of the
trajectory a water bridge was detected between the selected residues.

Alternatively, :meth:`~WaterBridgeAnalysis.count_by_type` can also be used to
generate the frequence of all water bridges in the simulation. ::

  w.count_by_type()

Returns ::

  [(0, 3, 'ARG', 1, 'O', 'ASP', 3, 'O', 0.5)]

For further data analysis, it is convenient to process the
:attr:`~WaterBridgeAnalysis.timeseries` data into a normalized table with the
:meth:`~WaterBridgeAnalysis.generate_table` method, which creates a new data
structure :attr:`WaterBridgeAnalysis.table` that contains one row for each
observation of a hydrogen bond::

  w.generate_table()

This table can then be easily turned into, e.g., a `pandas.DataFrame`_, and
further analyzed::

  import pandas as pd
  df = pd.DataFrame.from_records(w.table)


.. _pandas.DataFrame: http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html

Classes
-------

.. autoclass:: WaterBridgeAnalysis
   :members:

   .. attribute:: timesteps

      List of the times of each timestep. This can be used together with
      :attr:`~WaterBridgeAnalysis.timeseries` to find the specific time point
      of a water bridge existence, or see :attr:`~WaterBridgeAnalysis.table`.

   .. attribute:: table

      A normalised table of the data in
      :attr:`WaterBridgeAnalysis.timeseries`, generated by
      :meth:`WaterBridgeAnalysis.generate_table`. It is a
      :class:`numpy.recarray` with the following columns:

          0. "time"
          1. "donor_index"
          2. "acceptor_index"
          3. "donor_resnm"
          4. "donor_resid"
          5. "donor_atom"
          6. "acceptor_resnm"
          7. "acceptor_resid"
          8. "acceptor_atom"
          9. "distance"
          10. "angle"

      It takes up more space than :attr:`~WaterBridgeAnalysis.timeseries` but
      it is easier to analyze and to import into databases or dataframes.


      .. rubric:: Example

      For example, to create a `pandas.DataFrame`_ from ``h.table``::

         import pandas as pd
         df = pd.DataFrame.from_records(w.table)
"""
from __future__ import absolute_import, division
import six

from collections import defaultdict
import logging
import warnings

import numpy as np

from .hbond_analysis import HydrogenBondAnalysis
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
from MDAnalysis.lib.log import ProgressMeter, _set_verbose
from MDAnalysis.lib import distances
from MDAnalysis import SelectionWarning

logger = logging.getLogger('MDAnalysis.analysis.wbridges')

class WaterBridgeAnalysis(HydrogenBondAnalysis):
    """Perform a water bridge analysis

    The analysis of the trajectory is performed with the
    :meth:`WaterBridgeAnalysis.run` method. The result is stored in
    :attr:`WaterBridgeAnalysis.timeseries`. See
    :meth:`~WaterBridgeAnalysis.run` for the format.

    :class:`WaterBridgeAnalysis` uses the same default atom names as
    :class:`~MDAnalysis.analysis.hbonds.hbond_analysis.HydrogenBondAnalysis`,
    see :ref:`Default atom names for hydrogen bonding analysis`


    .. versionadded:: 0.17.0
    """
    def __init__(self, universe, selection1='protein',
                 selection2='not resname SOL', water_selection='resname SOL',
                 selection1_type='both', update_selection1=False,
                 update_selection2=False, update_water_selection=True,
                 filter_first=True, distance_type='hydrogen', distance=3.0,
                 angle=120.0, forcefield='CHARMM27', donors=None,
                 acceptors=None, start=None, stop=None, step=None, debug=None,
                 verbose=None):
        """Set up the calculation of water bridges between two selections in a
        universe.

        The timeseries is accessible as the attribute
        :attr:`WaterBridgeAnalysis.timeseries`.

        If no hydrogen bonds are detected or if the initial check fails, look
        at the log output (enable with :func:`MDAnalysis.start_logging` and set
        `verbose` ``=True``). It is likely that the default names for donors
        and acceptors are not suitable (especially for non-standard
        ligands). In this case, either change the `forcefield` or use
        customized `donors` and/or `acceptors`.

        Parameters
        ----------
        universe : Universe
            Universe object
        selection1 : str (optional)
            Selection string for first selection ['protein']
        selection2 : str (optional)
            Selection string for second selection ['not resname SOL']
            This string selects everything except water where water is assumed
            to have a residue name as SOL.
        water_selection : str (optional)
            Selection string for bridging water selection ['resname SOL']
            The default selection assumes that the water molecules have residue
            name "SOL". Change it to the appropriate selection for your
            specific force field.

            However, in theory this selection can be anything which forms
            hydrogen bond with selection 1 and selection 2.
        selection1_type : {"donor", "acceptor", "both"} (optional)
            Selection 1 can be 'donor', 'acceptor' or 'both'. Note that the
            value for `selection1_type` automatically determines how
            `selection2` handles donors and acceptors: If `selection1` contains
            'both' then `selection2` will also contain 'both'. If `selection1`
            is set to 'donor' then `selection2` is 'acceptor' (and vice versa).
            ['both'].
        update_selection1 : bool (optional)
            Update selection 1 at each frame. Setting to ``True`` if the
            selection is not static. Selection 1 is filtered first to speed up
            performance. Thus, setting to ``True`` is recommended if contact
            surface between selection 1 and selection 2 is constantly
            changing. [``False``]
        update_selection2 : bool (optional)
            Similiar to *update_selection1* but is acted upon selection 2.
            [``False``]
        update_water_selection : bool (optional)
            Update selection of water at each frame. Setting to ``False`` is
            **only** recommended when the total amount of water molecules in the
            simulation are small and when water molecules remain static across
            the simulation.

            However, in normal simulations, only a tiny proportion of water is
            engaged in the formation of water bridge. It is recommended to
            update the water selection and set keyword `filter_first` to
            ``True`` so as to filter out water not residing between the two
            selections. [``True``]
        filter_first : bool (optional)
            Filter the water selection to only include water within 3 *
            `distance` away from `both` selection 1 and selection 2.
            Selection 1 and selection 2 are both filtered to only include atoms
            3 * `distance` away from the other selection. [``True``]
        distance : float (optional)
            Distance cutoff for hydrogen bonds; only interactions with a H-A
            distance <= `distance` (and the appropriate D-H-A angle, see
            `angle`) are recorded. (Note: `distance_type` can change this to
            the D-A distance.) [3.0]
        angle : float (optional)
            Angle cutoff for hydrogen bonds; an ideal H-bond has an angle of
            180º.  A hydrogen bond is only recorded if the D-H-A angle is
            >=  `angle`. The default of 120º also finds fairly non-specific
            hydrogen interactions and a possibly better value is 150º. [120.0]
        forcefield : {"CHARMM27", "GLYCAM06", "other"} (optional)
            Name of the forcefield used. Switches between different
            :attr:`~HydrogenBondAnalysis.DEFAULT_DONORS` and
            :attr:`~HydrogenBondAnalysis.DEFAULT_ACCEPTORS` values.
            ["CHARMM27"]
        donors : sequence (optional)
            Extra H donor atom types (in addition to those in
            :attr:`~HydrogenBondAnalysis.DEFAULT_DONORS`), must be a sequence.
        acceptors : sequence (optional)
            Extra H acceptor atom types (in addition to those in
            :attr:`~HydrogenBondAnalysis.DEFAULT_ACCEPTORS`), must be a
            sequence.
        start : int (optional)
            starting frame-index for analysis, ``None`` is the first one, 0.
            `start` and `stop` are 0-based frame indices and are used to slice
            the trajectory (if supported) [``None``]
        stop : int (optional)
            last trajectory frame for analysis, ``None`` is the last one
            [``None``]
        step : int (optional)
            read every `step` between `start` (included) and `stop` (excluded),
            ``None`` selects 1. [``None``]
        distance_type : {"hydrogen", "heavy"} (optional)
            Measure hydrogen bond lengths between donor and acceptor heavy
            attoms ("heavy") or between donor hydrogen and acceptor heavy
            atom ("hydrogen"). If using "heavy" then one should set the
            *distance* cutoff to a higher value such as 3.5 Å. ["hydrogen"]
        debug : bool (optional)
            If set to ``True`` enables per-frame debug logging. This is disabled
            by default because it generates a very large amount of output in
            the log file. (Note that a logger must have been started to see
            the output, e.g. using :func:`MDAnalysis.start_logging`.)
        verbose : bool (optional)
            Toggle progress output. (Can also be given as keyword argument to
            :meth:`run`.)

        Notes
        -----
        In order to speed up processing, atoms are filtered by a coarse
        distance criterion before a detailed hydrogen bonding analysis is
        performed (`filter_first` = ``True``).

        If selection 1 and selection 2 are very mobile during the simulation
        and the contact surface is constantly changing (i.e. residues are
        moving farther than 3 x `distance`), you might consider setting the
        `update_selection1` and `update_selection2` keywords to ``True`` to
        ensure correctness.
        """
        self.water_selection = water_selection
        self.update_water_selection = update_water_selection
        super(WaterBridgeAnalysis, self).__init__(
            universe=universe, selection1=selection1, selection2=selection2,
            selection1_type=selection1_type,
            update_selection1=update_selection1,
            update_selection2=update_selection2, filter_first=filter_first,
            distance_type=distance_type, distance=distance, angle=angle,
            forcefield=forcefield, donors=donors, acceptors=acceptors,
            start=start, stop=stop, step=step, debug=debug, verbose=verbose)
        self._update_water_selection()

    def _log_parameters(self):
        """Log important parameters to the logfile."""
        logger.info("WBridge analysis: selection1 = %r (update: %r)",
        self.selection1, self.update_selection1)
        logger.info("WBridge analysis: selection2 = %r (update: %r)",
        self.selection2, self.update_selection2)
        logger.info("WBridge analysis: water selection = %r (update: %r)",
        self.water_selection, self.update_water_selection)
        logger.info("WBridge analysis: criterion: donor %s atom and acceptor \
        atom distance <= %.3f A", self.distance_type,
                    self.distance)
        logger.info("WBridge analysis: criterion: angle D-H-A >= %.3f degrees",
        self.angle)
        logger.info("WBridge analysis: force field %s to guess donor and \
        acceptor names", self.forcefield)
        logger.info("WBridge analysis: bonded hydrogen detection algorithm: \
        %r", self.detect_hydrogens)

    def _update_selection_1(self):
        self._s1 = self.u.select_atoms(self.selection1)
        self._s2 = self.u.select_atoms(self.selection2)
        if self.filter_first and self._s1:
            self.logger_debug('Size of selection 1 before filtering:'
                              ' {} atoms'.format(len(self._s1)))
            ns_selection_1 = AtomNeighborSearch(self._s1)
            self._s1 = ns_selection_1.search(self._s2, 3. * self.distance)
        self.logger_debug("Size of selection 1: {0} atoms".format(len(self._s1)))
        self._s1_donors = {}
        self._s1_donors_h = {}
        self._s1_acceptors = {}
        if self.selection1_type in ('donor', 'both'):
            self._s1_donors = self._s1.select_atoms(
                'name {0}'.format(' '.join(self.donors)))
            self._s1_donors_h = {}
            for i, d in enumerate(self._s1_donors):
                tmp = self._get_bonded_hydrogens(d)
                if tmp:
                    self._s1_donors_h[i] = tmp
            self.logger_debug("Selection 1 donors: {0}".format(len(self._s1_donors)))
            self.logger_debug("Selection 1 donor hydrogens: {0}".format(len(self._s1_donors_h)))
        if self.selection1_type in ('acceptor', 'both'):
            self._s1_acceptors = self._s1.select_atoms(
                'name {0}'.format(' '.join(self.acceptors)))
            self.logger_debug("Selection 1 acceptors: {0}".format(len(self._s1_acceptors)))

    def _sanity_check(self, selection, htype):
        """sanity check the selections 1 and 2

        *selection* is 1 or 2, *htype* is "donors" or "acceptors"
        """
        assert selection in (1, 2)
        assert htype in ("donors", "acceptors")
        atoms = getattr(self, "_s{0}_{1}".format(selection, htype))
        update = getattr(self, "update_selection{0}".format(selection))
        if not atoms:
            errmsg = "No {1} found in selection {0}. " \
                "You might have to specify a custom '{1}' keyword.".format(
                selection, htype)
            warnings.warn(errmsg, category=SelectionWarning)
            logger.warning(errmsg)

    def _update_water_selection(self):
        self._water = self.u.select_atoms(self.water_selection)
        self.logger_debug('Size of water selection before filtering:'
                          ' {} atoms'.format(len(self._water)))
        if self.filter_first:
            ns_water_selection = AtomNeighborSearch(self._water)
            self._water = ns_water_selection.search(self._s1,
                                                    3. * self.distance)
            self._water = ns_water_selection.search(self._s2,
                                                    3. * self.distance)

        self.logger_debug("Size of water selection: {0} atoms".format(len(self._water)))
        self._water_donors = {}
        self._water_donors_h = {}
        self._water_acceptors = {}
        if not self._water:
            logger.warning("Water selection '{0}' did not select any atoms."
                           .format(str(self.water_selection)[:80]))
        else:
            self._water_donors = self._water.select_atoms(
                'name {0}'.format(' '.join(self.donors)))
            self._water_donors_h = {}
            for i, d in enumerate(self._water_donors):
                tmp = self._get_bonded_hydrogens(d)
                if tmp:
                    self._water_donors_h[i] = tmp
            self.logger_debug("Water donors: {0}".format(len(self._water_donors)))
            self.logger_debug("Water donor hydrogens: {0}".format(len(self._water_donors_h)))
            self._water_acceptors = self._water.select_atoms(
                'name {0}'.format(' '.join(self.acceptors)))
            self.logger_debug("Water acceptors: {0}".format(len(self._water_acceptors)))

    def run(self, **kwargs):
        """Analyze trajectory and produce timeseries.

        Stores the water bridge data per frame as
        :attr:`WaterBridgeAnalysis.timeseries` (see there for output
        format).

        Parameters
        ----------
        verbose : bool (optional)
             toggle progress meter output
             :class:`~MDAnalysis.lib.log.ProgressMeter` [``True``]
        debug : bool (optional)
             enable detailed logging of debugging information; this can create
             *very big* log files so it is disabled (``False``) by default;
             setting `debug` toggles the debug status for
             :class:`WaterBridgeAnalysis`, namely the value of
             :attr:`WaterBridgeAnalysis.debug`.

        See Also
        --------
        :meth:`WaterBridgeAnalysis.generate_table` :
               processing the data into a different format.
        """
        logger.info("WBridge analysis: starting")
        logger.debug("WBridge analysis: donors    %r", self.donors)
        logger.debug("WBridge analysis: acceptors %r", self.acceptors)
        logger.debug("WBridge analysis: water bridge %r", self.water_selection)

        debug = kwargs.pop('debug', None)
        if debug is not None and debug != self.debug:
            self.debug = debug
            logger.debug("Toggling debug to %r", self.debug)
        if not self.debug:
            logger.debug("WBridge analysis: For full step-by-step debugging output use debug=True")

        self._timeseries = []
        self.timesteps = []
        self._water_network = []

        logger.info("checking trajectory...")  # n_frames can take a while!
        try:
            frames = np.arange(self.u.trajectory.n_frames)[self.traj_slice]
        except:
            logger.error("Problem reading trajectory or trajectory slice incompatible.")
            logger.exception()
            raise
        verbose = _set_verbose(verbose=kwargs.get('verbose', None),
                               quiet=kwargs.get('quiet', None),
                               default=True)
        pm = ProgressMeter(len(frames),
                           format="WBridge frame {current_step:5d}: {step:5d}/{numsteps} [{percentage:5.1f}%]\r",
                           verbose=verbose)

        logger.info("Starting analysis (frame index start=%d stop=%d, step=%d)",
                    (self.traj_slice.start or 0),
                    (self.traj_slice.stop or self.u.trajectory.n_frames),
                    self.traj_slice.step or 1)

        for progress, ts in enumerate(self.u.trajectory[self.traj_slice]):
            # all bonds for this timestep

            # dict of tuples (atom.index, atom.index) for quick check if
            # we already have the bond (to avoid duplicates)

            frame = ts.frame
            timestep = ts.time
            self.timesteps.append(timestep)
            pm.echo(progress, current_step=frame)
            self.logger_debug("Analyzing frame %(frame)d, timestep %(timestep)f ps", vars())
            if self.update_selection1:
                self._update_selection_1()
            if self.update_selection2:
                self._update_selection_2()
            if self.update_water_selection:
                self._update_water_selection()

            s1_frame_results_dict = defaultdict(list)
            if (self.selection1_type in ('donor', 'both') and
                self._water_acceptors):

                self.logger_debug("Selection 1 Donors <-> Water Acceptors")
                ns_acceptors = AtomNeighborSearch(self._water_acceptors)
                for i, donor_h_set in self._s1_donors_h.items():
                    d = self._s1_donors[i]
                    for h in donor_h_set:
                        res = ns_acceptors.search(h, self.distance)
                        for a in res:
                            donor_atom = h if self.distance_type != 'heavy' else d
                            dist = distances.calc_distance(donor_atom.position,
                                                           a.position)
                            if dist <= self.distance:
                                angle = distances.calc_angle(d.position, h.position,
                                                             a.position)
                                if angle >= self.angle:
                                    self.logger_debug(
                                        "S1-D: {0!s} <-> W-A: {1!s} {2:f} A, {3:f} DEG"\
                                        .format(h.index, a.index, dist, angle))
                                    s1_frame_results_dict[(a.resname, a.resid)].append(
                                        (h.index, a.index,
                                        (h.resname, h.resid, h.name),
                                        (a.resname, a.resid, a.name),
                                        dist, angle))

            if (self.selection1_type in ('acceptor', 'both') and
                self._s1_acceptors):

                self.logger_debug("Selection 1 Acceptors <-> Water Donors")
                ns_acceptors = AtomNeighborSearch(self._s1_acceptors)
                for i, donor_h_set in self._water_donors_h.items():
                    d = self._water_donors[i]
                    for h in donor_h_set:
                        res = ns_acceptors.search(h, self.distance)
                        for a in res:
                            donor_atom = h if self.distance_type != 'heavy' else d
                            dist = distances.calc_distance(donor_atom.position,
                                                           a.position)
                            if dist <= self.distance:
                                angle = distances.calc_angle(d.position, h.position,
                                                             a.position)
                                if angle >= self.angle:
                                    self.logger_debug(
                                        "S1-A: {0!s} <-> W-D: {1!s} {2:f} A, {3:f} DEG"\
                                        .format(a.index, h.index, dist, angle))
                                    s1_frame_results_dict[(h.resname, h.resid)].append(
                                        (h.index, a.index,
                                        (h.resname, h.resid, h.name),
                                        (a.resname, a.resid, a.name),
                                        dist, angle))

            # Narrow down the water selection
            selection_resn_id = list(s1_frame_results_dict.keys())
            if not selection_resn_id:
                self._timeseries.append([])
                continue
            selection_resn_id = ['(resname {} and resid {})'.format(
                resname, resid) for resname, resid in selection_resn_id]
            water_bridges = self._water.select_atoms(' or '.join(selection_resn_id))
            self.logger_debug("Size of water bridge selection: {0} atoms".format(len(water_bridges)))
            if not water_bridges:
                logger.warning("No water forming hydrogen bonding with selection 1.")
            water_bridges_donors = water_bridges.select_atoms(
                'name {0}'.format(' '.join(self.donors)))
            water_bridges_donors_h = {}
            for i, d in enumerate(water_bridges_donors):
                tmp = self._get_bonded_hydrogens(d)
                if tmp:
                    water_bridges_donors_h[i] = tmp
            self.logger_debug("water bridge donors: {0}".format(len(water_bridges_donors)))
            self.logger_debug("water bridge donor hydrogens: {0}".format(len(water_bridges_donors_h)))
            water_bridges_acceptors = water_bridges.select_atoms(
                'name {0}'.format(' '.join(self.acceptors)))
            self.logger_debug("water bridge: {0}".format(len(water_bridges_acceptors)))

            # Finding the hydrogen bonds between water bridge and selection 2
            s2_frame_results_dict = defaultdict(list)
            if self._s2_acceptors:
                self.logger_debug("Water bridge Donors <-> Selection 2 Acceptors")
                ns_acceptors = AtomNeighborSearch(self._s2_acceptors)
                for i, donor_h_set in water_bridges_donors_h.items():
                    d = water_bridges_donors[i]
                    for h in donor_h_set:
                        res = ns_acceptors.search(h, self.distance)
                        for a in res:
                            donor_atom = h if self.distance_type != 'heavy'  else d
                            dist = distances.calc_distance(donor_atom.position,
                                                           a.position)
                            if dist <= self.distance:
                                angle = distances.calc_angle(d.position, h.position,
                                                             a.position)
                                if angle >= self.angle:
                                    self.logger_debug(
                                        "WB-D: {0!s} <-> S2-A: {1!s} {2:f} A, {3:f} DEG"\
                                        .format(h.index, a.index, dist, angle))
                                    s2_frame_results_dict[(h.resname, h.resid)].append(
                                        (h.index, a.index,
                                        (h.resname, h.resid, h.name),
                                        (a.resname, a.resid, a.name),
                                        dist, angle))

            if water_bridges_acceptors:
                self.logger_debug("Selection 2 Donors <-> Selection 2 Acceptors")
                ns_acceptors = AtomNeighborSearch(water_bridges_acceptors)
                for i, donor_h_set in self._s2_donors_h.items():
                    d = self._s2_donors[i]
                    for h in donor_h_set:
                        res = ns_acceptors.search(h, self.distance)
                        for a in res:
                            donor_atom = h if self.distance_type != 'heavy' else d
                            dist = distances.calc_distance(donor_atom.position,
                                                           a.position)
                            if dist <= self.distance:
                                angle = distances.calc_angle(d.position, h.position,
                                                             a.position)
                                if angle >= self.angle:
                                    self.logger_debug(
                                        "WB-A: {0!s} <-> S2-D: {1!s} {2:f} A, {3:f} DEG"\
                                        .format(a.index, h.index, dist, angle))
                                    s2_frame_results_dict[(a.resname, a.resid)].append(
                                        (h.index, a.index,
                                        (h.resname, h.resid, h.name),
                                        (a.resname, a.resid, a.name),
                                        dist, angle))

            # Generate the water network
            water_network = {}
            for key in s2_frame_results_dict:
                s1_frame_results = set(s1_frame_results_dict[key])
                s2_frame_results = set(s2_frame_results_dict[key])
                if len(s1_frame_results.union(s2_frame_results)) > 1:
                    # Thus if selection 1 and selection 2 are the same and both
                    # only form a single hydrogen bond with a water, this entry
                    # won't be included.
                    water_network[key] = [s1_frame_results,
                    s2_frame_results.difference(s1_frame_results)]
            # Generate frame_results
            frame_results = []
            for s1_frame_results, s2_frame_results in water_network.values():
                frame_results.extend(list(s1_frame_results))
                frame_results.extend(list(s2_frame_results))

            self._timeseries.append(frame_results)
            self._water_network.append(water_network)


        logger.info("WBridge analysis: complete; timeseries  %s.timeseries",
                    self.__class__.__name__)
    @staticmethod
    def _reformat_hb(hb, atomformat="{0[0]!s}{0[1]!s}:{0[2]!s}"):
        """
        .. deprecated:: 1.0
            This is a compatibility layer so the water bridge
            _timeseries to timeserie conversion can perform as it does in the
            hydrogen bond analysis.

            For more information about the function of _reformat_hb in hydrogen
            bond analysis, please refer to the _reformat_hb in the hydrogen
            bond analysis module.

            As the _timeseries to timeserie conversion will be deprecated in 1.0
            this function will automatically lose its value.
        """

        return (list(hb[:2])
                + [atomformat.format(hb[2]), atomformat.format(hb[3])]
                + list(hb[4:]))

    @property
    def timeseries(self):
        r'''Time series of water bridges.

        Due to the intrinsic complexity of water bridge problem, where several
        atoms :math:`n` can be linked by the same water. A simple ``atom 1 -
        water bridge - atom 2`` mapping will create a very large list
        :math:`n!/(2(n-2)!)` due to the rule of combination.

        Thus, the output is arranged based on each individual bridging water
        allowing the later reconstruction of the water network. The hydrogen
        bonds from selection 1 to the **first** bridging water is followed by
        the hydrogen bonds from selection 2 to the **same** bridging water.
        After that, hydrogen bonds from selection 1 to the **second** bridging
        water is succeeded by hydrogen bonds from selection 2 to the **same
        second** bridging water. An example of the output is given in
        :ref:`wb_Analysis_Output`.

        Note
        ----
        To find an acceptor atom in :attr:`Universe.atoms` by
        *index* one would use ``u.atoms[acceptor_index]``.

        The :attr:`timeseries` is a managed attribute and it is generated
        from the underlying data in :attr:`_timeseries` every time the
        attribute is accessed. It is therefore costly to call and if
        :attr:`timeseries` is needed repeatedly it is recommended that you
        assign to a variable::

           w = WaterBridgeAnalysis(u)
           w.run()
           timeseries = w.timeseries

        See Also
        --------
        :attr:`table` : structured array of the data
        '''
        return super(WaterBridgeAnalysis, self).timeseries

    def generate_table(self):
        """Generate a normalised table of the results.

        The table is stored as a :class:`numpy.recarray` in the
        attribute :attr:`~WaterBridgeAnalysis.table`.

        See Also
        --------
        WaterBridgeAnalysis.table
        """
        super(WaterBridgeAnalysis, self).generate_table()

    def count_by_type(self):
        """Counts the frequency of water bridge of a specific type.

        If one atom *A* from *selection 1* is linked to atom *B* from
        *selection 2* through a bridging water, an entity will be created and
        the proportion of time that this linkage exists in the whole simulation
        will be calculated.

        Returns a :class:`numpy.recarray` containing atom indices for *A* and
        *B*, residue names, residue numbers, atom names (for both A and B) and
        the fraction of the total time during which the water bridge was
        detected. This method returns None if method
        :meth:`WaterBridgeAnalysis.run` was not executed first.

        Returns
        -------
        counts : numpy.recarray
             Each row of the array contains data to define a unique water
             bridge together with the frequency (fraction of the total time)
             that it has been observed.
        """
        if not self._has_timeseries():
            return

        wbridges = defaultdict(int)
        for wframe in self._water_network:
            pairs = set([])
            for water, (selection1, selection2) in wframe.items():
                for (donor_index, acceptor_index, donor, acceptor, distance,
                     angle) in selection1:
                    if donor[:2] == water:
                        sele1 = acceptor
                        sele1_index = acceptor_index
                    else:
                        sele1 = donor
                        sele1_index = donor_index
                    sele1_resnm, sele1_resid, sele1_atom = sele1

                    for (donor_index, acceptor_index, donor, acceptor,
                         distance, angle) in selection2:
                        if donor[:2] == water:
                            sele2 = acceptor
                            sele2_index = acceptor_index
                        else:
                            sele2 = donor
                            sele2_index = donor_index
                        sele2_resnm, sele2_resid, sele2_atom = sele2

                        key = (sele1_index, sele2_index)
                        if key in pairs:
                            continue
                        else:
                            pairs.add(key)
                            wb_key = (sele1_index, sele2_index, sele1_resnm,
                            sele1_resid, sele1_atom, sele2_resnm, sele2_resid,
                            sele2_atom)
                            wbridges[wb_key] += 1

        # build empty output table
        dtype = [
            ("sele1_index", int), ("sele2_index", int), ('sele1_resnm', 'U4'),
            ('sele1_resid', int), ('sele1_atom', 'U4'), ('sele2_resnm', 'U4'),
            ('sele2_resid', int), ('sele2_atom', 'U4'),
            ('frequency', float)
        ]
        out = np.empty((len(wbridges),), dtype=dtype)

        tsteps = float(len(self.timesteps))
        for cursor, (key, count) in enumerate(six.iteritems(wbridges)):
            out[cursor] = key + (count / tsteps,)

        r = out.view(np.recarray)
        return r
