# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2018 The MDAnalysis Development Team and contributors
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
:Year: 2017-2018
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

A higher order water bridge is defined as more than one water bridging
hydrogen bond acceptor and donor. An example of a second order water bridge:

e.g. -CO\ :sub:`2`\ :sup:`-`:···H−O:···H−O:···HN- (where H−O is part of **H−O**\ −H)

The following keyword arguments are important to control the behavior of the
water bridge analysis:

 - *water_selection* (``resname SOL``): the selection string for the bridging
   water
 - *order* the maximum number of water bridging both ends
 - donor-acceptor *distance* (Å): 3.0
 - Angle *cutoff* (degrees): 120.0
 - *forcefield* to switch between default values for different force fields
 - *donors* and *acceptors* atom types (to add additional atom names)

.. _wb_Analysis_Output:

Output
------

The results are a list of hydrogen bonds between the selection 1 or selection 2
and the bridging waters.

Each list is formated similar to the \ :attr:`HydrogenBondAnalysis.timeseries
<MDAnalysis.analysis.hbonds.hbond_analysis.HydrogenBondAnalysis.timeseries>`
and contains

  - the **identities** of donor and acceptor atoms,
  - the **distance** between the heavy atom acceptor atom and the hydrogen atom
  - the **angle** donor-hydrogen-acceptor angle (180º is linear).

Water bridge data are returned per frame, which is stored in \
:attr:`WaterBridgeAnalysis.timeseries` (In the following description, ``#``
indicates comments that are not part of the output.)::

    results = [
        [ # frame 1
           # hbonds linking the selection 1 and selection 2 to the bridging
           # water 1
           [[ # hbond 1 from selection 1 to the bridging water 1
              <donor index (0-based)>,
              <acceptor index (0-based)>, <donor identifer>, <acceptor identifer>,
              <distance>, <angle>
            ],
            [ # hbond 2 from bridging water 1 to the selection 2
              <donor index (0-based)>,
              <acceptor index (0-based)>, <donor identifer>, <acceptor identifer>,
              <distance>, <angle>
            ]],

           # hbonds linking the selection 1 and selection 2 to the bridging
           # water 1 and bridging water 2
           [[ # hbond 3 from selection 1 to the bridging water 1
              <donor index (0-based)>,
              <acceptor index (0-based)>, <donor identifer>, <acceptor identifer>,
              <distance>, <angle>
            ],
            [ # hbond 4 from bridging water 1 to the bridging water 2
              <donor index (0-based)>,
              <acceptor index (0-based)>, <donor identifer>, <acceptor identifer>,
              <distance>, <angle>
            ],
            [ # hbond 5 from bridging water 2 to the selection 2
              <donor index (0-based)>,
              <acceptor index (0-based)>, <donor identifer>, <acceptor identifer>,
              <distance>, <angle>
            ]],
           ....
        ],
        [ # frame 2
          [ ... ], [ ... ], ...
        ],
        ...
    ]


Detection of water bridges
---------------------------
Water bridges are recorded if a bridging water simultaneously forms
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

The maximum number of bridging waters detected can be changed using the order keyword. ::

  w = MDAnalysis.analysis.hbonds.WaterBridgeAnalysis(u, 'resname ARG', 'resname ASP',
                                                     order=3)

Thus, a maximum of three bridging waters will be detected.
The results are stored as the attribute
:attr:`WaterBridgeAnalysis.timeseries`; see :ref:`wb_Analysis_Output` for the
format.

An example of using the :attr:`~WaterBridgeAnalysis` would be
detecting the percentage of time a certain water bridge exits.

Trajectory :code:`u` has two frames, where the first frame contains a water
bridge from the oxygen of the first arginine to one of the oxygens in the carboxylic
group of aspartate (ASP3:OD1). In the second frame, the same water bridge forms but
is between the oxygen of the arginine and the other oxygen in the carboxylic
group (ASP3:OD2). ::

  print(w.timeseries)

prints out. ::

  [ # frame 1
    # A water bridge SOL2 links O from ARG1 to the carboxylic group OD1 of ASP3
   [[[0,1,('ARG',1,'O'),  ('SOL',2,'HW1'),  3.0,180],
     [2,3,('SOL',2,'HW2'),('ASP',3,'OD1'),  3.0,180],],
   ],
    # frame 2
    # Another water bridge SOL2 links O from ARG1 to the other oxygen of the
    # carboxylic group OD2 of ASP3
   [[[0,1,('ARG',1,'O'),  ('SOL',2,'HW1'),  3.0,180],
     [2,4,('SOL',2,'HW2'),('ASP',3,'OD2'),  3.0,180],],
   ],
  ]


.. _wb_count_by_type:

Use count_by_type
---------------------------

To calculate the percentage, we can use the :meth:`~WaterBridgeAnalysis.count_by_type` to
generate the frequence of all water bridges in the simulation. ::

  w.count_by_type()

Returns ::

  [(0, 3, 'ARG', 1, 'O', 'ASP', 3, 'OD1', 0.5),
   (0, 4, 'ARG', 1, 'O', 'ASP', 3, 'OD2', 0.5),]

You might think that the OD1 and OD2 are the same oxygen and the aspartate has just flipped
and thus, they should be counted as the same type of water bridge. The type of the water
bridge can be customised by supplying an analsysis function to
:meth:`~WaterBridgeAnalysis.count_by_type`.  ::

  def analysis(current, output):
      '''This function defines how the type of water bridge should be specified.

        Parameters
        ----------
        current : list
            The current water bridge being analysed is a list of hydrogen bonds from
            selection 1 to selection 2.
        output : dict
            A dictionary where the key is the type of the water bridge and the value
            is the number of this type of water bridge.
      '''

      # decompose the first hydrogen bond.
      s1_index, to_index, (s1_resname, s1_resid, s1_name),
      (to_resname, to_resid, to_name), dist, angle = current[0]
      # decompose the last hydrogen bond.
      from_index, s2_index, (from_resname, from_resid, from_name),
      (s2_resname, s2_resid, s2_name), dist, angle = current[-1]
      # if the residue name is ASP and the atom name is OD2 or OD1,
      # the atom name is changed to OD
      if s2_resname == 'ASP' and (s2_name == 'OD1' or s2_name == 'OD2'):
          s2_name = 'OD'
      # setting up the key which defines this type of water bridge.
      key = (s1_resname, s1_resid, s1_name, s2_resname, s2_resid, s2_name)
      # The number of this type of water bridge is incremented by 1.
      output[key] += 1

  w.count_by_type(analysis_func=analysis)

Returns ::

  [('ARG', 1, 'O', 'ASP', 3, 'OD', 1.0),]

Some people might only interested in contacts between residues and pay no attention
to the details regarding the atom name. This can also be achieved by supplying an analysis
function to :meth:`~WaterBridgeAnalysis.count_by_type`.  ::

  def analysis(current, output):
      '''This function defines how the type of water bridge should be specified.

        Parameters
        ----------
        current : list
            The current water bridge being analysed is a list of hydrogen bonds from
            selection 1 to selection 2.
        output : dict
            A dictionary where the key is the type of the water bridge and the value
            is the number of this type of water bridge.
      '''

      s1_index, to_index, (s1_resname, s1_resid, s1_name),
      (to_resname, to_resid, to_name), dist, angle = current[0]
      from_index, s2_index, (from_resname, from_resid, from_name),
      (s2_resname, s2_resid, s2_name), dist, angle = current[-1]
      # s1_name and s2_name are not included in the key
      key = (s1_resname, s1_resid, s2_resname, s2_resid)
      output[key] += 1

  w.count_by_type(analysis_func=analysis)

Returns ::

  [('ARG', 1, 'ASP', 3, 1.0),]

On the other hand, other people may insist that first order and second order water
bridges shouldn't be mixed together, which can also be achieved by supplying an analysis
function to :meth:`~WaterBridgeAnalysis.count_by_type`.  ::

  def analysis(current, output):
      '''This function defines how the type of water bridge should be specified.

        Parameters
        ----------
        current : list
            The current water bridge being analysed is a list of hydrogen bonds from
            selection 1 to selection 2.
        output : dict
            A dictionary where the key is the type of the water bridge and the value
            is the number of this type of water bridge.
      '''

      s1_index, to_index, (s1_resname, s1_resid, s1_name),
      (to_resname, to_resid, to_name), dist, angle = current[0]
      from_index, s2_index, (from_resname, from_resid, from_name),
      (s2_resname, s2_resid, s2_name), dist, angle = current[-1]
      # order of the current water bridge is computed
      order_of_water_bridge = len(current) - 1
      # and is included in the key
      key = (s1_resname, s1_resid, s2_resname, s2_resid, order_of_water_bridge)
      # The number of this type of water bridge is incremented by 1.
      output[key] += 1

  w.count_by_type(analysis_func=analysis)

The extra number 1 precede the 1.0 indicate that this is a first order water bridge ::

  [('ARG', 1, 'ASP', 3, 1, 1.0),]

Some people might not be interested in the interactions related to arginine. The undesirable
interactions can be discarded by supplying an analysis function to
:meth:`~WaterBridgeAnalysis.count_by_type`.  ::

  def analysis(current, output):
      '''This function defines how the type of water bridge should be specified.

        Parameters
        ----------
        current : list
            The current water bridge being analysed is a list of hydrogen bonds from
            selection 1 to selection 2.
        output : dict
            A dictionary where the key is the type of the water bridge and the value
            is the number of this type of water bridge.
      '''

      s1_index, to_index, (s1_resname, s1_resid, s1_name),
      (to_resname, to_resid, to_name), dist, angle = current[0]
      from_index, s2_index, (from_resname, from_resid, from_name),
      (s2_resname, s2_resid, s2_name), dist, angle = current[-1]
      if not s1_resname == 'ARG':
          key = (s1_resname, s1_resid, s2_resname, s2_resid)
          output[key] += 1

  w.count_by_type(analysis_func=analysis)

Returns nothing in this case ::

  [,]

Additional key words can be supplied to the analysis function by passing through
:meth:`~WaterBridgeAnalysis.count_by_type`.  ::

  def analysis(current, output, **kwargs):
      ...
  w.count_by_type(analysis_func=analysis, **kwargs)


.. _wb_count_by_time:

Use count_by_time
---------------------------

:meth:`~WaterBridgeAnalysis.count_by_type` aggregates data across frames, which
might be desirable in some cases but not the others. :meth:`~WaterBridgeAnalysis.count_by_time`
provides additional functionality for aggregating results for each frame.

The default behaviour of :meth:`~WaterBridgeAnalysis.count_by_time` is counting the number of
water bridges from selection 1 to selection 2 for each frame. Take the previous ASP, ARG salt
bridge for example:  ::

  w.count_by_time()

As one water bridge is found in both frames, the method returns ::

  [1, 1, ]

Similar to :meth:`~WaterBridgeAnalysis.count_by_type`
The behaviour of :meth:`~WaterBridgeAnalysis.count_by_time` can also be modified by supplying
an analysis function.

Suppose we want to count

  - the **number** of water molecules involved in bridging selection 1 to selection 2.
  - only if water bridge terminates in atom name **OD1 of ASP**.
  - only when water bridge is joined by less than **two** water.

The analysis function can be written as::

  def analysis(current, output, **kwargs):
      '''This function defines how the counting of water bridge should be specified.

        Parameters
        ----------
        current : list
            The current water bridge being analysed is a list of hydrogen bonds from
            selection 1 to selection 2.
        output : dict
            A dictionary where the key is the type of the water bridge and the value
            is the number of this type of water bridge.
            The output of this frame is the sum of all the values in this dictionary.
      '''

      # decompose the first hydrogen bond.
      s1_index, to_index, (s1_resname, s1_resid, s1_name),
      (to_resname, to_resid, to_name), dist, angle = current[0]
      # decompose the last hydrogen bond.
      from_index, s2_index, (from_resname, from_resid, from_name),
      (s2_resname, s2_resid, s2_name), dist, angle = current[-1]

      # only the residue name is ASP and the atom name is OD1,
      if s2_resname == 'ASP' and s2_name == 'OD1':
          # only if the order of water bridge is less than 2
          if len(current) -1 < 2:
              # extract all water molecules involved in the water bridge
              # extract the first water from selection 1
              s1_index, to_index, (s1_resname, s1_resid, s1_name),
              (to_resname, to_resid, to_name), dist, angle = current[0]
              key = (to_resname, to_resid)
              output[key] = 1

              # extract all the waters between selection 1 and selection 2
              for hbond in current[1:-1]:
                  # decompose the hydrogen bond.
                  from_index, to_index, (from_resname, from_resid, from_name),
                  (to_resname, to_resid, to_name), dist, angle = hbond
                  # add first water
                  key1 = (from_resname, from_resid)
                  output[key1] = 1
                  # add second water
                  key2 = (to_resname, to_resid)
                  output[key2] = 1

              # extract the last water to selection 2
              from_index, s2_index, (from_resname, from_resid, from_name),
              (s2_resname, s2_resid, s2_name), dist, angle = current[-1]
              key = (from_resname, from_resid)
              output[key] = 1

  w.count_by_time(analysis_func=analysis)

Returns ::

  [1, 0,]

Classes
-------

.. autoclass:: WaterBridgeAnalysis
   :members:

   .. attribute:: timesteps

      List of the times of each timestep. This can be used together with
      :attr:`~WaterBridgeAnalysis.timeseries` to find the specific time point
      of a water bridge existence.

"""
from __future__ import print_function, absolute_import

from collections import defaultdict
import logging
import warnings
import six
import numpy as np

from ..base import AnalysisBase
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
from MDAnalysis import NoDataError, MissingDataWarning, SelectionError
from MDAnalysis.lib import distances

logger = logging.getLogger('MDAnalysis.analysis.wbridges')

class WaterBridgeAnalysis(AnalysisBase):
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
    # use tuple(set()) here so that one can just copy&paste names from the
    # table; set() takes care for removing duplicates. At the end the
    # DEFAULT_DONORS and DEFAULT_ACCEPTORS should simply be tuples.

    #: default heavy atom names whose hydrogens are treated as *donors*
    #: (see :ref:`Default atom names for hydrogen bonding analysis`);
    #: use the keyword `donors` to add a list of additional donor names.
    DEFAULT_DONORS = {
        'CHARMM27': tuple(set([
            'N', 'OH2', 'OW', 'NE', 'NH1', 'NH2', 'ND2', 'SG', 'NE2', 'ND1', 'NZ', 'OG', 'OG1', 'NE1', 'OH'])),
        'GLYCAM06': tuple(set(['N', 'NT', 'N3', 'OH', 'OW'])),
        'other': tuple(set([]))}

    #: default atom names that are treated as hydrogen *acceptors*
    #: (see :ref:`Default atom names for hydrogen bonding analysis`);
    #: use the keyword `acceptors` to add a list of additional acceptor names.
    DEFAULT_ACCEPTORS = {
        'CHARMM27': tuple(set([
            'O', 'OC1', 'OC2', 'OH2', 'OW', 'OD1', 'OD2', 'SG', 'OE1', 'OE1', 'OE2', 'ND1', 'NE2', 'SD', 'OG', 'OG1', 'OH'])),
        'GLYCAM06': tuple(set(['N', 'NT', 'O', 'O2', 'OH', 'OS', 'OW', 'OY', 'SM'])),
        'other': tuple(set([]))}

    #: A :class:`collections.defaultdict` of covalent radii of common donors
    #: (used in :meth`_get_bonded_hydrogens_list` to check if a hydrogen is
    #: sufficiently close to its donor heavy atom). Values are stored for
    #: N, O, P, and S. Any other heavy atoms are assumed to have hydrogens
    #: covalently bound at a maximum distance of 1.5 Å.
    r_cov = defaultdict(lambda: 1.5,  # default value
                        N=1.31, O=1.31, P=1.58, S=1.55)

    def __init__(self, universe, selection1='protein',
                 selection2='not resname SOL', water_selection='resname SOL', order=1,
                 selection1_type='both', update_selection=False, update_water_selection=True,
                 filter_first=True, distance_type='hydrogen', distance=3.0,
                 angle=120.0, forcefield='CHARMM27', donors=None,
                 acceptors=None, debug=None, verbose=False, pbc=False, **kwargs):
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
        order : int (optional)
            The maximum number of water bridges linking both selections.
            if order is set to 3, then all the residues linked with less than
            three water molecules wil be detected. [1]

            Computation of high order water bridges can be very time consuming.
            Think carefully before running the calculation, do you really want
            to compute the 20th order water bridge between domain A and domain B
            or you just want to know the third order water bridge between two residues.
        selection1_type : {"donor", "acceptor", "both"} (optional)
            Selection 1 can be 'donor', 'acceptor' or 'both'. Note that the
            value for `selection1_type` automatically determines how
            `selection2` handles donors and acceptors: If `selection1` contains
            'both' then `selection2` will also contain 'both'. If `selection1`
            is set to 'donor' then `selection2` is 'acceptor' (and vice versa).
            ['both'].
        update_selection : bool (optional)
            Update selection 1 and 2 at each frame. Setting to ``True`` if the
            selection is not static. Selections are filtered first to speed up
            performance. Thus, setting to ``True`` is recommended if contact
            surface between selection 1 and selection 2 is constantly
            changing. [``False``]
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
            Filter the water selection to only include water within 4Å + `order` *
            (2Å + `distance`) away from `both` selection 1 and selection 2.
            Selection 1 and selection 2 are both filtered to only include atoms
            with the same distance away from the other selection. [``True``]
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
        moving farther than 4Å + `order` * (2Å + `distance`)), you might
        consider setting the `update_selection` keywords to ``True``
        to ensure correctness.
        """
        super(WaterBridgeAnalysis, self).__init__(universe.trajectory,
                                          **kwargs)
        self.water_selection = water_selection
        self.update_water_selection = update_water_selection
        # per-frame debugging output?
        self.debug = debug

        self.u = universe
        self.selection1 = selection1
        self.selection2 = selection2
        self.selection1_type = selection1_type

        # if the selection 1 and selection 2 are the same
        if selection1 == selection2:
            # eliminate the duplication
            self.selection1_type = "donor"
        self.update_selection = update_selection
        self.filter_first = filter_first
        self.distance = distance
        self.distance_type = distance_type  # note: everything except 'heavy' will give the default behavior
        self.angle = angle
        self.pbc = pbc and all(self.u.dimensions[:3])
        self.order = order

        # set up the donors/acceptors lists
        if donors is None:
            donors = []
        if acceptors is None:
            acceptors = []
        self.forcefield = forcefield
        self.donors = tuple(set(self.DEFAULT_DONORS[forcefield]).union(donors))
        self.acceptors = tuple(set(self.DEFAULT_ACCEPTORS[forcefield]).union(acceptors))

        if not (self.selection1 and self.selection2):
            raise ValueError('HydrogenBondAnalysis: invalid selections')
        elif self.selection1_type not in ('both', 'donor', 'acceptor'):
            raise ValueError('HydrogenBondAnalysis: Invalid selection type {0!s}'.format(self.selection1_type))

        self._network = []  # final result accessed as self.network
        self.timesteps = None  # time for each frame

        self._log_parameters()

    def _log_parameters(self):
        """Log important parameters to the logfile."""
        logger.info("WBridge analysis: selection = %r (update: %r)",
        self.selection2, self.update_selection)
        logger.info("WBridge analysis: water selection = %r (update: %r)",
        self.water_selection, self.update_water_selection)
        logger.info("WBridge analysis: criterion: donor %s atom and acceptor \
        atom distance <= %.3f A", self.distance_type,
                    self.distance)
        logger.info("WBridge analysis: criterion: angle D-H-A >= %.3f degrees",
        self.angle)
        logger.info("WBridge analysis: force field %s to guess donor and \
        acceptor names", self.forcefield)

    def _update_selection(self):
        self._s1 = self.u.select_atoms(self.selection1)
        self._s2 = self.u.select_atoms(self.selection2)

        if self.filter_first and self._s1:
            self.logger_debug('Size of selection 1 before filtering:'
                              ' {} atoms'.format(len(self._s1)))
            ns_selection_1 = AtomNeighborSearch(self._s1, box=self.box)
            self._s1 = ns_selection_1.search(self._s2, self.selection_distance)
        self.logger_debug("Size of selection 1: {0} atoms".format(len(self._s1)))

        if not self._s1:
            logger.warning('Selection 1 "{0}" did not select any atoms.'
                           .format(str(self.selection1)[:80]))
            return



        if self.filter_first and self._s2:
            self.logger_debug('Size of selection 2 before filtering:'
                              ' {} atoms'.format(len(self._s2)))
            ns_selection_2 = AtomNeighborSearch(self._s2, box=self.box)
            self._s2 = ns_selection_2.search(self._s1, self.selection_distance)
        self.logger_debug('Size of selection 2: {0} atoms'.format(len(self._s2)))

        if not self._s2:
            logger.warning('Selection 2 "{0}" did not select any atoms.'
                           .format(str(self.selection2)[:80]))
            return


        if self.selection1_type in ('donor', 'both'):
            self._s1_donors = self._s1.select_atoms(
                'name {0}'.format(' '.join(self.donors)))
            self._s1_donors_h = {}
            for i, atom in enumerate(self._s1_donors):
                tmp = self._get_bonded_hydrogens(atom)
                if tmp:
                    self._s1_donors_h[i] = tmp
            self.logger_debug("Selection 1 donors: {0}".format(len(self._s1_donors)))
            self.logger_debug("Selection 1 donor hydrogens: {0}".format(len(self._s1_donors_h)))
        if self.selection1_type in ('acceptor', 'both'):
            self._s1_acceptors = self._s1.select_atoms(
                'name {0}'.format(' '.join(self.acceptors)))
            self.logger_debug("Selection 1 acceptors: {0}".format(len(self._s1_acceptors)))

        if not self._s2:
            return None
        if self.selection1_type in ('donor', 'both'):
            self._s2_acceptors = self._s2.select_atoms(
                'name {0}'.format(' '.join(self.acceptors)))
            self.logger_debug("Selection 2 acceptors: {0:d}".format(len(self._s2_acceptors)))
        if self.selection1_type in ('acceptor', 'both'):
            self._s2_donors = self._s2.select_atoms(
                'name {0}'.format(' '.join(self.donors)))
            self._s2_donors_h = {}
            for i, d in enumerate(self._s2_donors):
                tmp = self._get_bonded_hydrogens(d)
                if tmp:
                    self._s2_donors_h[i] = tmp
            self.logger_debug("Selection 2 donors: {0:d}".format(len(self._s2_donors)))
            self.logger_debug("Selection 2 donor hydrogens: {0:d}".format(len(self._s2_donors_h)))

    def _update_water_selection(self):
        self._water = self.u.select_atoms(self.water_selection)
        self.logger_debug('Size of water selection before filtering:'
                          ' {} atoms'.format(len(self._water)))
        if self._water and self.filter_first:
            filtered_s1 = AtomNeighborSearch(self._water, box=self.box).search(self._s1,
                                                    self.selection_distance)
            if filtered_s1:
                self._water = AtomNeighborSearch(filtered_s1, box=self.box).search(self._s2,
                                                        self.selection_distance)

        self.logger_debug("Size of water selection: {0} atoms".format(len(self._water)))

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

    def _sanity_check(self):
        """sanity check the selections 1 and 2
        *selection* is 1 or 2, *htype* is "donors" or "acceptors"
        If selections do not update and the required donor and acceptor
        selections are empty then a :exc:`SelectionError` is immediately
        raised.
        If selections update dynamically then it is possible that the selection
        will yield donors/acceptors at a later step and we only issue a
        warning.
        .. versionadded:: 0.11.0
        """

        for s in (1,2):
            for htype in ('donors', "acceptors"):
                atoms = getattr(self, "_s{0}_{1}".format(s, htype))
                if not atoms and not getattr(self, "update_selection") and not getattr(self, "update_water_selection"):
                    errmsg = '''No {1} found in selection {0}. \
You might have to specify a custom '{1}' keyword.'''.format(s, htype)
                    raise SelectionError(errmsg)

    def _get_bonded_hydrogens(self, atom):
        """Find hydrogens bonded within cutoff to `atom`.

        Hydrogens are detected by either name ("H*", "[123]H*") or type ("H");
        this is not fool-proof as the atom type is not always a character but
        the name pattern should catch most typical occurrences.

        The distance from `atom` is calculated for all hydrogens in the residue
        and only those within a cutoff are kept. The cutoff depends on the
        heavy atom (more precisely, on its element, which is taken as the first
        letter of its name ``atom.name[0]``) and is parameterized in
        :attr:`HydrogenBondAnalysis.r_cov`. If no match is found then the
        default of 1.5 Å is used.


        Parameters
        ----------
        atom : groups.Atom
             heavy atom

        Returns
        -------
        hydrogen_atoms : AtomGroup or []
            list of hydrogens (can be a :class:`~MDAnalysis.core.groups.AtomGroup`)
            or empty list ``[]`` if none were found.
        """
        try:
            return atom.residue.atoms.select_atoms(
                "(name H* 1H* 2H* 3H* or type H) and around {0:f} name {1!s}"
                "".format(self.r_cov[atom.name[0]], atom.name))
        except NoDataError:
            return []

    def logger_debug(self, *args):
        if self.debug:
            logger.debug(*args)


    def _prepare(self):
        # The distance for selection is defined as twice the maximum bond length of an O-H bond (2A)
        # plus order of water bridge times the length of OH bond plus hydrogne bond distance
        self.selection_distance = 1.0 * (2 * 2 + self.order * (2 + self.distance))

        self._s1_donors = {}
        self._s1_donors_h = {}
        self._s1_acceptors = {}

        self._s2_donors = {}
        self._s2_donors_h = {}
        self._s2_acceptors = {}

        self.box = self.u.dimensions if self.pbc else None
        self._update_selection()

        self._water_donors = {}
        self._water_donors_h = {}
        self._water_acceptors = {}

        self.timesteps = []
        if self._s1 and self._s2:
            self._update_water_selection()
            logger.info("WBridge analysis: no atoms found in the selection.")

        logger.info("WBridge analysis: initial checks passed.")

        logger.info("WBridge analysis: starting")
        logger.debug("WBridge analysis: donors    %r", self.donors)
        logger.debug("WBridge analysis: acceptors %r", self.acceptors)
        logger.debug("WBridge analysis: water bridge %r", self.water_selection)

        if self.debug:
            logger.debug("Toggling debug to %r", self.debug)
        else:
            logger.debug("WBridge analysis: For full step-by-step debugging output use debug=True")

        logger.info("Starting analysis (frame index start=%d stop=%d, step=%d)",
                    self.start, self.stop, self.step)

        self._sanity_check()

    def _donor2acceptor(self, donors, donor_hs, acceptor):
        result = []
        ns_acceptors = AtomNeighborSearch(acceptor, self.box)
        for i, donor_h_set in donor_hs.items():
            d = donors[i]
            for h in donor_h_set:
                res = ns_acceptors.search(h, self.distance)
                for a in res:
                    donor_atom = h if self.distance_type != 'heavy' else d
                    dist = distances.calc_bonds(donor_atom.position,
                                                a.position, box=self.box)
                    if dist <= self.distance:
                        angle = distances.calc_angles(d.position, h.position,
                                                      a.position, box=self.box)
                        angle = np.rad2deg(angle)
                        if angle >= self.angle:
                            self.logger_debug(
                                "D: {0!s} <-> A: {1!s} {2:f} A, {3:f} DEG" \
                                    .format(h.index, a.index, dist, angle))
                            result.append((h.index, a.index,
                                     (h.resname, h.resid, h.name),
                                     (a.resname, a.resid, a.name),
                                     dist, angle))
        return result



    def _single_frame(self):
        self.timesteps.append(self._ts.time)
        self.box = self.u.dimensions if self.pbc else None

        if self.update_selection:
            self._update_selection()
        if self._s1 and self._s2:
            if self.update_water_selection:
                self._update_water_selection()
        else:
            self._network.append(defaultdict(dict))
            return


        selection_1 = []
        water_pool = defaultdict(list)
        next_round_water = set([])
        selection_2 = []

        if self.selection1_type in ('donor', 'both'):
            # check for direct hbond from s1 to s2
            if self._s2_acceptors:
                self.logger_debug("Selection 1 Donors <-> Selection 2 Acceptors")
                results = self._donor2acceptor(self._s1_donors, self._s1_donors_h, self._s2_acceptors)
                for line in results:
                    h_index, a_index, (h_resname, h_resid, h_name), (a_resname, a_resid, a_name), dist, angle = line
                    water_pool[(a_resname, a_resid)] = None
                    selection_1.append(line)
                    selection_2.append((a_resname, a_resid))
            if self._water_acceptors:
                self.logger_debug("Selection 1 Donors <-> Water Acceptors")
                results = self._donor2acceptor(self._s1_donors, self._s1_donors_h, self._water_acceptors)
                for line in results:
                    h_index, a_index, (h_resname, h_resid, h_name), (a_resname, a_resid, a_name), dist, angle = line
                    next_round_water.add((a_resname, a_resid))
                    selection_1.append(line)



        if (self.selection1_type in ('acceptor', 'both') and
                self._s1_acceptors):
            self.logger_debug("Selection 1 Acceptors <-> Water Donors")
            results = self._donor2acceptor(self._water_donors, self._water_donors_h, self._s1_acceptors)
            for line in results:
                h_index, a_index, (h_resname, h_resid, h_name), (a_resname, a_resid, a_name), dist, angle = line
                line = (a_index, h_index, (a_resname, a_resid, a_name), (h_resname, h_resid, h_name), dist, angle)
                next_round_water.add((h_resname, h_resid))
                selection_1.append(line)

            self.logger_debug("Selection 2 Donors <-> Selection 1 Acceptors")
            results = self._donor2acceptor(self._s2_donors, self._s2_donors_h, self._s1_acceptors)
            for line in results:
                h_index, a_index, (h_resname, h_resid, h_name), (a_resname, a_resid, a_name), dist, angle = line
                line = (a_index, h_index, (a_resname, a_resid, a_name), (h_resname, h_resid, h_name), dist, angle)
                water_pool[(h_resname, h_resid)] = None
                selection_1.append(line)
                selection_2.append((h_resname, h_resid))

        for i in range(self.order):
            # Narrow down the water selection
            selection_resn_id = list(next_round_water)
            if not selection_resn_id:
                self._network.append(defaultdict(dict))
                logger.warning("No water forming hydrogen bonding with selection 1.")
                return
            selection_resn_id = ['(resname {} and resid {})'.format(
                resname, resid) for resname, resid in selection_resn_id]
            water_bridges = self._water.select_atoms(' or '.join(selection_resn_id))
            self.logger_debug("Size of water bridge selection for the {} order of water bridge: {} atoms".format(i+1,
                                                                                                len(water_bridges)))

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
            self.logger_debug("water bridge acceptors: {0}".format(len(water_bridges_acceptors)))

            if i < self.order:
                next_round_water = set([])
                # first check if the water can touch the other end
                # Finding the hydrogen bonds between water bridge and selection 2
                if self._s2_acceptors:
                    self.logger_debug("Order {} water donor <-> Selection 2 Acceptors".format(i + 1))
                    results = self._donor2acceptor(water_bridges_donors, water_bridges_donors_h, self._s2_acceptors)
                    for line in results:
                        h_index, a_index, (h_resname, h_resid, h_name), (a_resname, a_resid, a_name), dist, angle = line
                        water_pool[(h_resname, h_resid)].append(line)
                        selection_2.append((a_resname, a_resid))

                # Finding the hydrogen bonds between water bridge and selection 2
                if water_bridges_acceptors:
                    self.logger_debug("Selection 2 Donors <-> Order {} water".format(i + 1))
                    results = self._donor2acceptor(self._s2_donors, self._s2_donors_h, water_bridges_acceptors)
                    for line in results:
                        h_index, a_index, (h_resname, h_resid, h_name), (a_resname, a_resid, a_name), dist, angle = line
                        line = (a_index, h_index, (a_resname, a_resid, a_name), (h_resname, h_resid, h_name), dist, angle)
                        water_pool[(a_resname, a_resid)].append(line)
                        selection_2.append((h_resname, h_resid))

                # find the water water hydrogen bond
                if water_bridges_acceptors:
                    self.logger_debug("Order {} water acceptor <-> Order {} water donor".format(i + 1, i + 2))
                    results = self._donor2acceptor(self._water_donors, self._water_donors_h, water_bridges_acceptors)
                    for line in results:
                        h_index, a_index, (h_resname, h_resid, h_name), (a_resname, a_resid, a_name), dist, angle = line
                        line = (a_index, h_index, (a_resname, a_resid, a_name), (h_resname, h_resid, h_name), dist, angle)
                        water_pool[(a_resname, a_resid)].append(line)
                        next_round_water.add((h_resname, h_resid))

                if water_bridges_donors_h:
                    self.logger_debug("Order {} water donor <-> Order {} water acceptor".format(i + 1, i + 2))
                    results = self._donor2acceptor(water_bridges_donors, water_bridges_donors_h, self._water_acceptors)
                    for line in results:
                        h_index, a_index, (h_resname, h_resid, h_name), (a_resname, a_resid, a_name), dist, angle = line
                        water_pool[(h_resname, h_resid)].append(line)
                        next_round_water.add((a_resname, a_resid))

        # solve the connectivity network
        result = {'start': defaultdict(dict), 'water': defaultdict(dict)}

        def add_route(result, route):
            if len(route) == 1:
                result['start'][route[0]] = None
            else:
                # exclude the the selection which goes back to itself
                if (sorted(route[0][:2]) == sorted(route[-1][:2])):
                    return

                # selection 2 to water
                result['water'][route[-1]] = None
                # water to water
                for i in range(1, len(route) - 1):
                    result['water'][route[i]][route[i+1]] = result['water'][route[i+1]]
                # selection 1 to water
                result['start'][route[0]][route[1]] = result['water'][route[1]]

        def traverse_water_network(graph, node, end, route, maxdepth, result):
            if len(route) > self.order + 1:
                return
            else:
                if node in end:
                    # check if any duplication happens
                    if len(route) == len(set(route)):
                        add_route(result, route)
                else:
                    for new_node in graph[node]:
                        new_route = route[:]
                        new_route.append(new_node)
                        new_node = new_node[3][:2]
                        traverse_water_network(graph, new_node, end, new_route, maxdepth, result)

        for s1 in selection_1:

            route = [s1, ]
            next_mol = s1[3][:2]
            traverse_water_network(water_pool, next_mol, selection_2, route[:], self.order, result)

        self._network.append(result['start'])

    def _traverse_water_network(self, graph, current, analysis_func=None, output=None, link_func=None, **kwargs):
        if link_func is None:
            link_func = self._full_link

        if graph is None:
            if not analysis_func is None:
                analysis_func(current, output, **kwargs)
        else:
            # make sure no loop can occur
            if len(current) <= self.order:
                for node in graph:
                    new = link_func(current, node)
                    self._traverse_water_network(graph[node], new, analysis_func, output, link_func, **kwargs)

    @staticmethod
    def _reformat_hb(hb, atomformat="{0[0]!s}{0[1]!s}:{0[2]!s}"):
        """Convert 0.16.1 _timeseries hbond item to 0.16.0 hbond item.
        In 0.16.1, donor and acceptor are stored as a tuple(resname,
        resid, atomid). In 0.16.0 and earlier they were stored as a string.
        .. deprecated:: 1.0
           This is a compatibility layer so that we can provide the same output
           in timeseries as before. However, for 1.0 we should make timeseries
           just return _timeseries, i.e., change the format of timeseries to
           the un-ambiguous representation provided in _timeseries.
        """
        return (list(hb[:2])
                + [atomformat.format(hb[2]), atomformat.format(hb[3])]
                + list(hb[4:]))

    @property
    def timeseries(self):
        r'''Time series of water bridges.

        The output is arranged based on each individual link from selection 1
        to selection 2 allowing the later reconstruction of the water network.
        Each entity is a list of hydrogen bonds from selection 1 to selection 2.
        An example of the output is given in :ref:`wb_Analysis_Output`.

        Note
        ----
        To find an acceptor atom in :attr:`Universe.atoms` by
        *index* one would use ``u.atoms[acceptor_index]``.

        The :attr:`timeseries` is a managed attribute and it is generated
        from the underlying data in :attr:`_network` every time the
        attribute is accessed. It is therefore costly to call and if
        :attr:`timeseries` is needed repeatedly it is recommended that you
        assign to a variable::

           w = WaterBridgeAnalysis(u)
           w.run()
           timeseries = w.timeseries

        '''

        def analysis(current, output):
            output.extend(current)

        timeseries = []
        for frame in self._network:
            new_frame = []
            self._traverse_water_network(frame, [], analysis_func=analysis, output=new_frame, link_func=self._full_link)
            timeseries.append(new_frame)
        self._timeseries = timeseries
        # after 1.0 will be return the unformated _timeseries
        return [[self._reformat_hb(hb) for hb in hframe] for hframe in self._timeseries]

    @classmethod
    def _full_link(self, output, node):
        result = output[:]
        result.append(node)
        return result

    @classmethod
    def _count_by_type_analysis(self, current, output):
        '''
        Generates the key for count_by_type analysis.
        :return:
        '''
        s1_index, to_index, (s1_resname, s1_resid, s1_name), (to_resname, to_resid, to_name), dist, angle = \
            current[0]
        from_index, s2_index, (from_resname, from_resid, from_name), (s2_resname, s2_resid, s2_name), dist, angle = \
            current[-1]
        key = (s1_index, s2_index, s1_resname, s1_resid, s1_name, s2_resname, s2_resid, s2_name)
        output[key] += 1

    def count_by_type(self, analysis_func=None, output='expand', **kwargs):
        """Counts the frequency of water bridge of a specific type.

        If one atom *A* from *selection 1* is linked to atom *B* from
        *selection 2* through one or more bridging waters, an entity will be created and
        the proportion of time that this linkage exists in the whole simulation
        will be calculated.

        The identification of a specific type water bridge can be modified by
        supplying the analysis_func function. See :ref:`wb_count_by_type`
        for detail.

        Returns
        -------
        counts : list
            Returns a :class:`list` containing atom indices for *A* and
            *B*, residue names, residue numbers, atom names (for both A and B) and
            the fraction of the total time during which the water bridge was
            detected. This method returns None if method
            :meth:`WaterBridgeAnalysis.run` was not executed first.


        """
        if analysis_func is None:
            analysis_func = self._count_by_type_analysis

        if self._network:
            length = len(self._network)
            result_dict = defaultdict(int)
            for frame in self._network:
                self._traverse_water_network(frame, [], analysis_func=analysis_func, output=result_dict,
                                             link_func=self._full_link, **kwargs)
            result = [[i for i in key] for key in result_dict]
            [result[i].append(result_dict[key]*1.0/length) for i, key in enumerate(result_dict)]
            return result
        else:
            return None

    @classmethod
    def _count_by_time_analysis(self, current, output):
        s1_index, to_index, (s1_resname, s1_resid, s1_name), (to_resname, to_resid, to_name), dist, angle = \
            current[0]
        from_index, s2_index, (from_resname, from_resid, from_name), (s2_resname, s2_resid, s2_name), dist, angle = \
            current[-1]
        key = (s1_index, s2_index, s1_resname, s1_resid, s1_name, s2_resname, s2_resid, s2_name)
        output[key] += 1

    def count_by_time(self, analysis_func=None, **kwargs):
        """Counts the number of water bridges per timestep.

        The counting behaviour can be adjusted by supplying analysis_func.
        See :ref:`wb_count_by_time` for details.

        Returns
        -------
        counts : list
             Returns a time series ``N(t)`` where ``N`` is the total
             number of observed water bridges at time ``t``.

        """
        if analysis_func is None:
            analysis_func = self._count_by_time_analysis
        if self._network:
            result = []
            for time, frame in zip(self.timesteps, self._network):
                result_dict = defaultdict(int)
                self._traverse_water_network(frame, [], analysis_func=analysis_func, output=result_dict,
                                             link_func=self._full_link, **kwargs)
                result.append((time, sum([result_dict[key] for key in result_dict])))
            return result
        else:
            return None

    @classmethod
    def _timesteps_by_type_analysis(self, current, output, **kwargs):
        s1_index, to_index, (s1_resname, s1_resid, s1_name), (to_resname, to_resid, to_name), dist, angle = \
            current[0]
        from_index, s2_index, (from_resname, from_resid, from_name), (s2_resname, s2_resid, s2_name), dist, angle = \
            current[-1]
        key = (s1_index, s2_index, s1_resname, s1_resid, s1_name, s2_resname, s2_resid, s2_name)
        output[key].append(kwargs.pop('time'))

    def timesteps_by_type(self, analysis_func=None, **kwargs):
        """Frames during which each water bridges existed, sorted by each water bridges.

        Processes :attr:`WaterBridgeAnalysis._network` and returns a
        :class:`numpy.recarray` containing atom indices, residue names, residue
        numbers (for donors and acceptors) and each timestep at which the
        hydrogen bond was detected.

        In principle, this is the same as :attr:`~WaterBridgeAnalysis.table`
        but sorted by hydrogen bond and with additional data for the
        *donor_heavy_atom* and angle and distance omitted.


        Returns
        -------
        data : numpy.recarray



        """
        if analysis_func is None:
            analysis_func = self._timesteps_by_type_analysis

        if self._network:
            result = defaultdict(list)
            for time, frame in zip(self.timesteps, self._network):
                self._traverse_water_network(frame, [], analysis_func=analysis_func, output=result,
                                             link_func=self._full_link, time=time, **kwargs)

            result_list = []
            for key, time_list in six.iteritems(result):
                for time in time_list:
                    key = list(key)
                    key.append(time)
                    result_list.append(key)
            return result_list
        else:
            return None

    def generate_table(self):
        """Generate a normalised table of the results.

        The table is stored as a :class:`numpy.recarray` in the
        attribute :attr:`~HydrogenBondAnalysis.table`.

        See Also
        --------
        HydrogenBondAnalysis.table

        """
        if self._network is None:
            msg = "No data computed, do run() first."
            warnings.warn(msg, category=MissingDataWarning)
            logger.warning(msg)
            return
        if not hasattr(self, '_timeseries'):
            self.timeseries
        timeseries = self._timeseries

        num_records = np.sum([len(hframe) for hframe in timeseries])
        # build empty output table
        dtype = [
            ("time", float),
            ("donor_index", int),  ("acceptor_index", int),
            ("donor_resnm", "|U4"), ("donor_resid", int), ("donor_atom", "|U4"),
            ("acceptor_resnm", "|U4"), ("acceptor_resid", int), ("acceptor_atom", "|U4"),
            ("distance", float), ("angle", float)]
        # according to Lukas' notes below, using a recarray at this stage is ineffective
        # and speedups of ~x10 can be achieved by filling a standard array, like this:
        out = np.empty((num_records,), dtype=dtype)
        cursor = 0  # current row
        for t, hframe in zip(self.timesteps, timeseries):
            for (donor_index, acceptor_index, donor,
                 acceptor, distance, angle) in hframe:
                # donor|acceptor = (resname, resid, atomid)
                out[cursor] = (t, donor_index, acceptor_index) + \
                donor + acceptor + (distance, angle)
                cursor += 1
        assert cursor == num_records, "Internal Error: Not all HB records stored"
        self.table = out.view(np.recarray)
        logger.debug("HBond: Stored results as table with %(num_records)d entries.", vars())