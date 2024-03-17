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

r"""Water Bridge analysis --- :mod:`MDAnalysis.analysis.hydrogenbonds.WaterBridgeAnalysis`
==========================================================================================

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

e.g. -CO\ :sub:`2`\ :sup:`-`:···H−O:···H−O:···HN-
(where H−O is part of **H−O**\ −H)

The following keyword arguments are important to control the behaviour of the
water bridge analysis:

- *water_selection* (``resname SOL``): the selection string for the bridging water
- *order* the maximum number of water bridging both ends
- donor-acceptor *distance* (Å): 3.0
- Angle *cutoff* (degrees): 120.0
- *forcefield* to switch between default values for different force fields
- *donors* and *acceptors* atom types (to add additional atom names)

Theory
------

This module attempts to find multi-order water bridge by an approach similar
to breadth-first search, where the first solvation shell of selection 1 is
selected, followed by the selection of the second solvation shell as well as
any hydrogen bonding partner from selection 1. After that, the third solvation
shell, as well as any binding partners from selection 2, are detected. This
process is repeated until the maximum order of water bridges is reached.

.. _wb_Analysis_Network:

Output as Network
-----------------

Since the waters connecting the two ends of the selections are by nature a
network. We provide a network representation of the water network. Water bridge
data are returned per frame, which is stored in
:attr:`WaterBridgeAnalysis.results.network`. Each frame is represented as a
dictionary, where the keys are the hydrogen bonds originating from selection
1 and the values are new dictionaries representing the hydrogen bonds coming
out of the corresponding molecules making hydrogen bonds with selection 1.

As for the hydrogen bonds which reach the selection 2, the values of the
corresponding keys are None. One example where selection 1 and selection 2 are
joined by one water molecule (A) which also hydrogen bond to another water (B)
which also hydrogen bond to selection 2 would be represented as ::

    # (selection 1)-O:···H-O(water 1):···H-(selection 2)
    #                      |             :
    #                      H·············O-H(water2)
    #                                    H
    {(sele1_acceptor, None, water1_donor, water1_donor_heavy, distance, angle):
         {(water1_acceptor, None, sele2_donor, sele2_donor_heavy,
         distance, angle): None},
         {(water1_donor, water1_donor_heavy, water2_acceptor, None,
         distance, angle):
              {(water2_acceptor, None, sele2_donor, sele2_donor_heavy,
              distance, angle): None}
          },
    }

The atoms are represented by atom index and if the atom is hydrogen bond donor,
it is followed by the index of the corresponding heavy atom
``(donor_proton, donor_heavy_atom)``.
If the atom is a hydrogen bond acceptor, it is followed by none.

.. _wb_Analysis_Timeseries:

Output as Timeseries
--------------------

For lower order water bridges, it might be desirable to represent the
connections as :attr:`WaterBridgeAnalysis.results.timeseries`. The results
are returned per frame and are a list of hydrogen bonds between the selection
1 or selection 2 and the bridging waters. Due to the complexity of the higher
order water bridge and the fact that one hydrogen bond between two waters can
appear in both third and fourth order water bridge, the hydrogen bonds in the
:attr:`WaterBridgeAnalysis.results.timeseries` attribute are generated in a
depth-first search manner to avoid duplication. Example code of how
:attr:`WaterBridgeAnalysis.results.timeseries` is generated::

    def network2timeseries(network, timeseries):
        '''Traverse the network in a depth-first fashion.
        expand_timeseries will expand the compact representation to the
        familiar timeseries representation.'''

        if network is None:
            return
        else:
            for node in network:
                timeseries.append(expand_timeseries(node))
                network2timeseries(network[node], timeseries)

    timeseries = []
    network2timeseries(network, timeseries)

An example would be. ::

    results = [
        [ # frame 1
           [ <donor index>, <acceptor index>,
            (<donor residue name>, <donor residue number>, <donor atom name>),
            (<acceptor residue name>, <acceptor residue number>,
            <acceptor atom name>),
             <distance>, <angle>],
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

Detection of water bridges
--------------------------
Water bridges are recorded if a bridging water simultaneously forms
hydrogen bonds with selection 1 and selection 2.

Hydrogen bonds are detected based on a geometric criterion:

1. The distance between acceptor and hydrogen is less than or equal to
   `distance` (default is 3 Å).

2. The angle between donor-hydrogen-acceptor is greater than or equal to
   `angle` (default is 120º).

The cut-off values `angle` and `distance` can be set as keywords to
:class:`WaterBridgeAnalysis`.

Donor and acceptor heavy atoms are detected from atom names. The current
defaults are appropriate for the CHARMM27 and GLYCAM06 force fields as defined
in Table `Default atom names for water bridge analysis`_.

Hydrogen atoms bonded to a donor are searched based on its distance to the
donor. The algorithm searches for all hydrogens
(name "H*" or name "[123]H" or type "H") in the same residue as the donor atom
within a cut-off distance of 1.2 Å.

.. _Default atom names for water bridge analysis:

.. table:: Default heavy atom names for CHARMM27 force field.

   =========== ==============  =========== ====================================
   group       donor           acceptor    comments
   =========== ==============  =========== ====================================
   main chain  N               O, OC1, OC2 OC1, OC2 from amber99sb-ildn
                                           (Gromacs)
   water       OH2, OW         OH2, OW     SPC, TIP3P, TIP4P (CHARMM27,Gromacs)

   ARG         NE, NH1, NH2
   ASN         ND2             OD1
   ASP                         OD1, OD2
   CYS         SG
   CYH                         SG          possible false positives for CYS
   GLN         NE2             OE1
   GLU                         OE1, OE2
   HIS         ND1, NE2        ND1, NE2    presence of H determines if donor
   HSD         ND1             NE2
   HSE         NE2             ND1
   HSP         ND1, NE2
   LYS         NZ
   MET                         SD          see e.g. :footcite:p:`Gregoret1991`
   SER         OG              OG
   THR         OG1             OG1
   TRP         NE1
   TYR         OH              OH
   =========== ==============  =========== ====================================

.. table:: Heavy atom types for GLYCAM06 force field.

   =========== =========== ==================
   element     donor       acceptor
   =========== =========== ==================
   N           N,NT,N3     N,NT
   O           OH,OW       O,O2,OH,OS,OW,OY
   S                       SM
   =========== =========== ==================

Donor and acceptor names for the CHARMM27 force field will also work for e.g.
OPLS/AA or amber (tested in Gromacs). Residue names in the table are for
information only and are not taken into account when determining acceptors and
donors. This can potentially lead to some ambiguity in the assignment of
donors/acceptors for residues such as histidine or cytosine.

For more information about the naming convention in GLYCAM06 have a look at the
`Carbohydrate Naming Convention in Glycam`_.

.. _`Carbohydrate Naming Convention in Glycam`:
   http://glycam.ccrc.uga.edu/documents/FutureNomenclature.htm

The lists of donor and acceptor names can be extended by providing lists of
atom names in the `donors` and `acceptors` keywords to
:class:`WaterBridgeAnalysis`. If the lists are entirely inappropriate
(e.g. when analysing simulations done with a force field that uses very
different atom names) then one should either use the value "other" for
`forcefield` to set no default values or derive a new class and set the
default list oneself::

 class WaterBridgeAnalysis_OtherFF(WaterBridgeAnalysis):
       DEFAULT_DONORS = {"OtherFF": tuple(set([...]))}
       DEFAULT_ACCEPTORS = {"OtherFF": tuple(set([...]))}

Then simply use the new class instead of the parent class and call it with
```forcefield` = "OtherFF"``. Please also consider contributing the list of
heavy atom names to MDAnalysis.


.. rubric:: References

.. footbibliography::


How to perform ``WaterBridgeAnalysis``
--------------------------------------

All water bridges between arginine and aspartic acid can be analysed with ::

  import MDAnalysis
  from MDAnalysis.analysis.hydrogenbonds import WaterBridgeAnalysis

  u = MDAnalysis.Universe('topology', 'trajectory')
  w = WaterBridgeAnalysis(u, 'resname ARG', 'resname ASP')
  w.run()

The maximum number of bridging waters detected can be changed using the order
keyword. ::

  w = WaterBridgeAnalysis(u, 'resname ARG', 'resname ASP', order=3)

Thus, a maximum of three bridging waters will be detected.

An example of using the :attr:`~WaterBridgeAnalysis` would be
detecting the percentage of time a certain water bridge exits.

Trajectory :code:`u` has two frames, where the first frame contains a water
bridge from the oxygen of the first arginine to one of the oxygens in the
carboxylic group of aspartate (ASP3:OD1). In the second frame, the same water
bridge forms but is between the oxygen of the arginine and the other oxygen in
the carboxylic group (ASP3:OD2). ::

  # index residue id residue name atom name
  #     0          1          ARG         O
  #     1          2          SOL        OW
  #     2          2          SOL       HW1
  #     3          2          SOL       HW2
  #     4          3          ASP       OD1
  #     5          3          ASP       OD2
  print(w.results.timeseries)

prints out ::

  [ # frame 1
    # A water bridge SOL2 links O from ARG1 to the carboxylic group OD1 of ASP3
   [[0,2,('ARG',1,  'O'),('SOL',2,'HW1'),  3.0,180],
    [3,4,('SOL',2,'HW2'),('ASP',3,'OD1'),  3.0,180],
   ],
    # frame 2
    # Another water bridge SOL2 links O from ARG1 to the other oxygen of the
    # carboxylic group OD2 of ASP3
   [[0,2,('ARG',1,  'O'),('SOL',2,'HW1'),  3.0,180],
    [3,5,('SOL',2,'HW2'),('ASP',3,'OD2'),  3.0,180],
   ],
  ]


.. _wb_count_by_type:

Use ``count_by_type``
---------------------

We can use the :meth:`~WaterBridgeAnalysis.count_by_type` to
generate the frequency of all water bridges in the simulation. ::

  w.count_by_type()

Returns ::

  [(0, 3, 'ARG', 1, 'O', 'ASP', 3, 'OD1', 0.5),
   (0, 4, 'ARG', 1, 'O', 'ASP', 3, 'OD2', 0.5),]

You might think that the OD1 and OD2 are the same oxygen and the aspartate has
just flipped and thus, they should be counted as the same type of water bridge.
The type of the water bridge can be customised by supplying an analysis
function to :meth:`~WaterBridgeAnalysis.count_by_type`.

The analysis function has two parameters. The current and the output. The
current is a list of hydrogen bonds from selection 1 to selection 2, formatted
in the same fashion as :attr:`WaterBridgeAnalysis.network`, and an example will
be ::

  [
  # sele1 acceptor idx,   , water donor index, donor heavy atom idx, dist, ang.
   [                 0, None,                 2,                   1, 3.0,180],
  # water donor idx, donor heavy atom idx, sele2 acceptor idx, distance, angle.
   [              3,                    1,                  4, None, 3.0,180],]

Where ``current[0]`` is the first hydrogen bond originating from selection 1
and ``current[-1]`` is the final hydrogen bond ending in selection 2. The
output sums up all the information in the current frame and is a dictionary
with a user-defined key and the value is the weight of the corresponding key.
During the analysis phase, the function analysis iterates through all the water
bridges and modify the output in-place. At the end of the analysis, the keys
from all the frames are collected and the corresponding values will be summed
up and returned. ::

  def analysis(current, output, u):
      r'''This function defines how the type of water bridge should be
      specified.

        Parameters
        ----------
        current : list
            The current water bridge being analysed is a list of hydrogen bonds
            from selection 1 to selection 2.
        output : dict
            A dictionary which is modified in-place where the key is the type
            of the water bridge and the value is the weight of this type of
            water bridge.
        u : MDAnalysis.universe
            The current Universe for looking up atoms.'''

      # decompose the first hydrogen bond.
      sele1_index, sele1_heavy_index, atom2, heavy_atom2, dist, angle =
      current[0]
      # decompose the last hydrogen bond.
      atom1, heavy_atom1, sele2_index, sele2_heavy_index, dist, angle =
      current[-1]
      # expand the atom index to the resname, resid, atom names
      sele1 = u.atoms[sele1_index]
      sele2 = u.atoms[sele2_index]
      (s1_resname, s1_resid, s1_name) = (sele1.resname, sele1.resid,
      sele1.name)
      (s2_resname, s2_resid, s2_name) = (sele2.resname, sele2.resid,
      sele2.name)
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

  [(('ARG', 1, 'O', 'ASP', 3, 'OD'), 1.0),]

Note that the result is arranged in the format of
``(key, the proportion of time)``. When no custom analysis function is supplied
, the key is expanded and is formatted as ::

  [('ARG', 1, 'O', 'ASP', 3, 'OD', 1.0),]

Some people might only interested in contacts between residues and pay no
attention to the details regarding the atom name. However, since multiple water
bridges can exist between two residues, which sometimes can give a result such
that the water bridge between two residues exists 300% of the time. Though this
might be a desirable result for some people, others might want the water bridge
between two residues to be only counted once per frame. This can also be
achieved by supplying an analysis function to
:meth:`~WaterBridgeAnalysis.count_by_type`. ::

  def analysis(current, output, u):
      '''This function defines how the type of water bridge should be specified
      .

        Parameters
        ----------
        current : list
            The current water bridge being analysed is a list of hydrogen bonds
            from selection 1 to selection 2.
        output : dict
            A dictionary which is modified in-place where the key is the type
            of the water bridge and the value is the weight of this type of
            water bridge.
        u : MDAnalysis.universe
            The current Universe for looking up atoms.
      '''

      # decompose the first hydrogen bond.
      sele1_index, sele1_heavy_index, atom2, heavy_atom2, dist, angle =
      current[0]
      # decompose the last hydrogen bond.
      atom1, heavy_atom1, sele2_index, sele2_heavy_index, dist, angle =
      current[-1]
      # expand the atom index to the resname, resid, atom names
      sele1 = u.atoms[sele1_index]
      sele2 = u.atoms[sele2_index]
      (s1_resname, s1_resid, s1_name) = (sele1.resname, sele1.resid,
      sele1.name)
      (s2_resname, s2_resid, s2_name) = (sele2.resname, sele2.resid,
      sele2.name)
      # s1_name and s2_name are not included in the key
      key = (s1_resname, s1_resid, s2_resname, s2_resid)

      # Each residue is only counted once per frame
      output[key] = 1

  w.count_by_type(analysis_func=analysis)

Returns ::

  [(('ARG', 1, 'ASP', 3), 1.0),]

On the other hand, other people may insist that the first order and
second-order water bridges shouldn't be mixed together, which can also be
achieved by supplying an analysis function to
:meth:`~WaterBridgeAnalysis.count_by_type`.  ::

  def analysis(current, output, u):
      '''This function defines how the type of water bridge should be specified
      .

        Parameters
        ----------
        current : list
            The current water bridge being analysed is a list of hydrogen bonds
            from selection 1 to selection 2.
        output : dict
            A dictionary which is modified in-place where the key is the type
            of the water bridge and the value is the weight of this type of
            water bridge.
        u : MDAnalysis.universe
            The current Universe for looking up atoms.
      '''

      # decompose the first hydrogen bond.
      sele1_index, sele1_heavy_index, atom2, heavy_atom2, dist, angle =
      current[0]
      # decompose the last hydrogen bond.
      atom1, heavy_atom1, sele2_index, sele2_heavy_index, dist, angle =
      current[-1]
      # expand the atom index to the resname, resid, atom names
      sele1 = u.atoms[sele1_index]
      sele2 = u.atoms[sele2_index]
      (s1_resname, s1_resid, s1_name) = (sele1.resname, sele1.resid,
      sele1.name)
      (s2_resname, s2_resid, s2_name) = (sele2.resname, sele2.resid,
      sele2.name)
      # order of the current water bridge is computed
      order_of_water_bridge = len(current) - 1
      # and is included in the key
      key = (s1_resname, s1_resid, s2_resname, s2_resid, order_of_water_bridge)
      # The number of this type of water bridge is incremented by 1.
      output[key] += 1

  w.count_by_type(analysis_func=analysis)

The extra number 1 precede the 1.0 indicate that this is a first order water
bridge ::

  [(('ARG', 1, 'ASP', 3, 1), 1.0),]

Some people might not be interested in the interactions related to arginine.
The undesirable interactions can be discarded by supplying an analysis function
to :meth:`~WaterBridgeAnalysis.count_by_type`.  ::

  def analysis(current, output, u):
      '''This function defines how the type of water bridge should be
      specified.

        Parameters
        ----------
        current : list
            The current water bridge being analysed is a list of hydrogen bonds
            from selection 1 to selection 2.
        output : dict
            A dictionary which is modified in-place where the key is the type
            of the water bridge and the value is the number of this type of
            water bridge.
        u : MDAnalysis.universe
            The current Universe for looking up atoms.
      '''

      # decompose the first hydrogen bond.
      sele1_index, sele1_heavy_index, atom2, heavy_atom2, dist, angle =
      current[0]
      # decompose the last hydrogen bond.
      atom1, heavy_atom1, sele2_index, sele2_heavy_index, dist, angle =
      current[-1]
      # expand the atom index to the resname, resid, atom names
      sele1 = u.atoms[sele1_index]
      sele2 = u.atoms[sele2_index]
      (s1_resname, s1_resid, s1_name) = (sele1.resname, sele1.resid,
      sele1.name)
      (s2_resname, s2_resid, s2_name) = (sele2.resname, sele2.resid,
      sele2.name)
      if not s1_resname == 'ARG':
          key = (s1_resname, s1_resid, s2_resname, s2_resid)
          output[key] += 1

  w.count_by_type(analysis_func=analysis)

Returns nothing in this case ::

  [,]

Additional keywords can be supplied to the analysis function by passing through
:meth:`~WaterBridgeAnalysis.count_by_type`.  ::

  def analysis(current, output, **kwargs):
      ...
  w.count_by_type(analysis_func=analysis, **kwargs)


.. _wb_count_by_time:

Use ``count_by_time``
---------------------

:meth:`~WaterBridgeAnalysis.count_by_type` aggregates data across frames, which
might be desirable in some cases but not the others.
:meth:`~WaterBridgeAnalysis.count_by_time` provides additional functionality
for aggregating results for each frame.

The default behaviour of :meth:`~WaterBridgeAnalysis.count_by_time` is counting
the number of water bridges from selection 1 to selection 2 for each frame.
Take the previous ASP, ARG salt bridge for example:  ::

  w.count_by_time()

As one water bridge is found in both frames, the method returns ::

  [(1.0, 1), (2.0, 1), ]

Similar to :meth:`~WaterBridgeAnalysis.count_by_type`
The behaviour of :meth:`~WaterBridgeAnalysis.count_by_time` can also be
modified by supplying an analysis function.

Suppose we want to count

  - the **number** of water molecules involved in bridging selection 1 to
    selection 2.
  - only if water bridge terminates in atom name **OD1 of ASP**.
  - only when water bridge is joined by less than **two** water.

The analysis function can be written as::

  def analysis(current, output, u, **kwargs):
      '''This function defines how the counting of water bridge should be
      specified.

        Parameters
        ----------
        current : list
            The current water bridge being analysed is a list of hydrogen bonds
            from selection 1 to selection 2.
        output : dict
            A dictionary which is modified in-place where the key is the type
            of the water bridge and the value is the number of this type of
            water  bridge.
        u : MDAnalysis.universe
            The current Universe for looking up atoms.
      '''

      # decompose the first hydrogen bond.
      sele1_index, sele1_heavy_index, atom2, heavy_atom2, dist, angle =
      current[0]
      # decompose the last hydrogen bond.
      atom1, heavy_atom1, sele2_index, sele2_heavy_index, dist, angle =
      current[-1]
      # expand the atom index to the resname, resid, atom names
      sele1 = u.atoms[sele1_index]
      sele2 = u.atoms[sele2_index]
      (s1_resname, s1_resid, s1_name) =
      (sele1.resname, sele1.resid, sele1.name)
      (s2_resname, s2_resid, s2_name) =
      (sele2.resname, sele2.resid, sele2.name)

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

  [(1.0, 1), (2.0, 0),]

Classes
-------

.. autoclass:: WaterBridgeAnalysis
   :members:

   .. attribute:: timesteps

      List of the times of each timestep. This can be used together with
      :attr:`~WaterBridgeAnalysis.results.timeseries` to find the specific
      time point of a water bridge existence.

   .. attribute:: results.network

      Network representation of the water network.

      .. versionadded:: 2.0.0

   .. attribute:: network

      Alias to the :attr:`results.network` attribute.

      .. deprecated:: 2.0.0
         Will be removed in MDAnalysis 3.0.0. Please use
         :attr:`results.network` instead.

   .. attribute:: table

      .. deprecated:: 2.0.0
         Will be removed in MDAnalysis 3.0.0. Please generate
         the table with :meth:`generate_table` instead.

   .. attribute:: results.timeseries

      List of hydrogen bonds between the selection 1 or selection 2
      and the bridging waters, for each frame.

      .. versionadded:: 2.0.0

   .. attribute:: timeseries

      Alias to the :attr:`results.timeseries` attribute.

      .. deprecated:: 2.0.0
         Will be removed in MDAnalysis 3.0.0. Please use
         :attr:`results.timeseries` instead.
"""
from collections import defaultdict
import logging
import warnings
import numpy as np

from ..base import AnalysisBase
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
from MDAnalysis.lib.distances import capped_distance, calc_angles
from MDAnalysis import NoDataError, MissingDataWarning

logger = logging.getLogger('MDAnalysis.analysis.WaterBridgeAnalysis')


class WaterBridgeAnalysis(AnalysisBase):
    """Perform a water bridge analysis

    The analysis of the trajectory is performed with the
    :meth:`WaterBridgeAnalysis.run` method. The result is stored in
    :attr:`WaterBridgeAnalysis.results.timeseries`. See
    :meth:`~WaterBridgeAnalysis.run` for the format.

    .. versionadded:: 0.17.0

    """
    # use tuple(set()) here so that one can just copy&paste names from the
    # table; set() takes care for removing duplicates. At the end the
    # DEFAULT_DONORS and DEFAULT_ACCEPTORS should simply be tuples.

    #: default heavy atom names whose hydrogens are treated as *donors*
    #: (see :ref:`Default atom names for water bridge analysis`);
    #: use the keyword `donors` to add a list of additional donor names.
    DEFAULT_DONORS = {
        'CHARMM27': tuple(
            {'N', 'OH2', 'OW', 'NE', 'NH1', 'NH2', 'ND2', 'SG', 'NE2', 'ND1',
             'NZ', 'OG', 'OG1', 'NE1', 'OH'}),
        'GLYCAM06': tuple({'N', 'NT', 'N3', 'OH', 'OW'}),
        'other': tuple(set([]))}

    #: default atom names that are treated as hydrogen *acceptors*
    #: (see :ref:`Default atom names for water bridge analysis`);
    #: use the keyword `acceptors` to add a list of additional acceptor names.
    DEFAULT_ACCEPTORS = {
        'CHARMM27': tuple(
            {'O', 'OC1', 'OC2', 'OH2', 'OW', 'OD1', 'OD2', 'SG', 'OE1', 'OE1',
             'OE2', 'ND1', 'NE2', 'SD', 'OG', 'OG1', 'OH'}),
        'GLYCAM06':
            tuple({'N', 'NT', 'O', 'O2', 'OH', 'OS', 'OW', 'OY', 'SM'}),
        'other': tuple(set([]))}

    #: A :class:`collections.defaultdict` of covalent radii of common donors
    #: (used in :meth`_get_bonded_hydrogens_list` to check if a hydrogen is
    #: sufficiently close to its donor heavy atom). Values are stored for
    #: N, O, P, and S. Any other heavy atoms are assumed to have hydrogens
    #: covalently bound at a maximum distance of 1.5 Å.
    r_cov = defaultdict(lambda: 1.5,  # default value
                        N=1.31, O=1.31, P=1.58, S=1.55)  # noqa: E741

    def __init__(self, universe, selection1='protein',
                 selection2='not resname SOL', water_selection='resname SOL',
                 order=1, selection1_type='both', update_selection=False,
                 update_water_selection=True, filter_first=True,
                 distance_type='hydrogen', distance=3.0, angle=120.0,
                 forcefield='CHARMM27', donors=None, acceptors=None,
                 output_format="sele1_sele2", debug=None,
                 pbc=False, **kwargs):
        """Set up the calculation of water bridges between two selections in a
        universe.

        The timeseries is accessible as the attribute
        :attr:`WaterBridgeAnalysis.results.timeseries`.

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

            However, in theory, this selection can be anything which forms
            a hydrogen bond with selection 1 and selection 2.
        order : int (optional)
            The maximum number of water bridges linking both selections.
            if the order is set to 3, then all the residues linked with less
            than three water molecules will be detected. [1]

            Computation of high order water bridges can be very time-consuming.
            Think carefully before running the calculation, do you really want
            to compute the 20th order water bridge between domain A and domain
            B or you just want to know the third order water bridge between two
            residues.
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
            **only** recommended when the total amount of water molecules in
            the simulation are small and when water molecules remain static
            across the simulation.

            However, in normal simulations, only a tiny proportion of water is
            engaged in the formation of water bridge. It is recommended to
            update the water selection and set keyword `filter_first` to
            ``True`` so as to filter out water not residing between the two
            selections. [``True``]
        filter_first : bool (optional)
            Filter the water selection to only include water within 4 Å +
            `order` * (2 Å + `distance`) away from `both` selection 1 and
            selection 2.
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
            hydrogen interactions and possibly better value is 150º. [120.0]
        forcefield : {"CHARMM27", "GLYCAM06", "other"} (optional)
            Name of the forcefield used. Switches between different
            :attr:`~DEFAULT_DONORS` and
            :attr:`~DEFAULT_ACCEPTORS` values.
            ["CHARMM27"]
        donors : sequence (optional)
            Extra H donor atom types (in addition to those in :attr:`~DEFAULT_DONORS`).
            This shall be the name of the heavy atom that is bonded to the hydrogen.
            For example, the oxygen ('O') in the hydroxyl group. Must be a sequence.
        acceptors : sequence (optional)
            Extra H acceptor atom types (in addition to those in
            :attr:`~DEFAULT_ACCEPTORS`), must be a
            sequence.
        distance_type : {"hydrogen", "heavy"} (optional)
            Measure hydrogen bond lengths between donor and acceptor heavy
            atoms ("heavy") or between donor hydrogen and acceptor heavy
            atom ("hydrogen"). If using "heavy" then one should set the
            *distance* cutoff to a higher value such as 3.5 Å. ["hydrogen"]
        output_format: {"sele1_sele2", "donor_acceptor"} (optional)
            Setting the output format for timeseries and table. If set to
            "sele1_sele2", for each hydrogen bond, the one close to selection 1
            will be placed before selection 2. If set to "donor_acceptor", the
            donor will be placed before acceptor. "sele1_sele2"]
        debug : bool (optional)
            If set to ``True`` enables per-frame debug logging. This is
            disabled by default because it generates a very large amount of
            output in the log file. (Note that a logger must have been started
            to see the output, e.g. using :func:`MDAnalysis.start_logging`.)
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
        moving farther than 4 Å + `order` * (2 Å + `distance`)), you might
        consider setting the `update_selection` keywords to ``True``
        to ensure correctness.

        .. versionchanged 0.20.0
           The :attr:`WaterBridgeAnalysis.timeseries` has been updated
           see :attr:`WaterBridgeAnalysis.timeseries` for detail.
           This class is now based on
           :class:`~MDAnalysis.analysis.base.AnalysisBase`.


        """
        super(WaterBridgeAnalysis, self).__init__(universe.trajectory,
                                                  **kwargs)
        self.water_selection = water_selection
        self.update_water_selection = update_water_selection
        # per-frame debugging output?
        self.debug = debug

        # set the output format
        self.output_format = output_format

        self.u = universe
        self.selection1 = selection1
        self.selection2 = selection2
        self.selection1_type = selection1_type
        if "selection2_type" in kwargs:
            raise ValueError("`selection2_type` is not a keyword argument.")

        # if the selection 1 and selection 2 are the same
        if selection1 == selection2:
            # eliminate the duplication
            self.selection1_type = "donor"
        self.update_selection = update_selection
        self.filter_first = filter_first
        self.distance = distance
        if distance_type not in {"hydrogen", "heavy"}:
            raise ValueError(f"Only 'hydrogen' and 'heavy' are allowed for option `distance_type' ({distance_type}).")
        self.distance_type = distance_type
        # will give the default behavior
        self.angle = angle
        self.pbc = pbc and all(self.u.dimensions[:3])
        self.order = order

        # set up the donors/acceptors lists
        if donors is None:
            donors = ()
        if acceptors is None:
            acceptors = ()
        self.forcefield = forcefield
        self.donors = tuple(set(self.DEFAULT_DONORS[forcefield]).union(donors))
        self.acceptors = tuple(set(self.DEFAULT_ACCEPTORS[forcefield]).union(
            acceptors))

        if self.selection1_type not in ('both', 'donor', 'acceptor'):
            raise ValueError('WaterBridgeAnalysis: '
                             'Invalid selection type {0!s}'.format(
                                self.selection1_type))

        # final result accessed as self.results.network
        self.results.network = []
        self.results.timeseries = None
        self.timesteps = None  # time for each frame

        self._log_parameters()

    def _log_parameters(self):
        """Log important parameters to the logfile."""
        logger.info("WaterBridgeAnalysis: selection = %r (update: %r)",
                    self.selection2, self.update_selection)
        logger.info("WaterBridgeAnalysis: water selection = %r (update: %r)",
                    self.water_selection, self.update_water_selection)
        logger.info("WaterBridgeAnalysis: criterion: donor %s atom and "
                    "acceptor atom distance <= %.3f A", self.distance_type,
                    self.distance)
        logger.info("WaterBridgeAnalysis: criterion: "
                    "angle D-H-A >= %.3f degrees",
                    self.angle)
        logger.info("WaterBridgeAnalysis: force field %s to guess donor and \
        acceptor names", self.forcefield)

    def _build_residue_dict(self, selection):
        # Build the residue_dict where the key is the residue name
        # The content is a dictionary where hydrogen bond donor heavy atom
        # names is the key
        # The content is the hydrogen bond donor hydrogen atom names
        atom_group = self.u.select_atoms(selection)
        for residue in atom_group.residues:
            if residue.resname not in self._residue_dict:
                self._residue_dict[residue.resname] = defaultdict(set)
            for atom in residue.atoms:
                if atom.name in self.donors:
                    self._residue_dict[residue.resname][atom.name].update(
                        self._get_bonded_hydrogens(atom).names)

    def _update_donor_h(self, atom_ix, h_donors, donors_h):
        atom = self.u.atoms[atom_ix]
        residue = atom.residue
        hydrogen_names = self._residue_dict[residue.resname][atom.name]
        if hydrogen_names:
            hydrogens = residue.atoms.select_atoms('name {0}'.format(
                ' '.join(hydrogen_names))).ix
            for atom in hydrogens:
                h_donors[atom] = atom_ix
                donors_h[atom_ix].append(atom)

    def _update_selection(self):
        self._s1_donors = []
        self._s1_h_donors = {}
        self._s1_donors_h = defaultdict(list)
        self._s1_acceptors = []

        self._s2_donors = []
        self._s2_h_donors = {}
        self._s2_donors_h = defaultdict(list)
        self._s2_acceptors = []

        self._s1 = self.u.select_atoms(self.selection1).ix
        self._s2 = self.u.select_atoms(self.selection2).ix

        if self.filter_first and len(self._s1):
            self.logger_debug('Size of selection 1 before filtering:'
                              ' {} atoms'.format(len(self._s1)))
            ns_selection_1 = AtomNeighborSearch(self.u.atoms[self._s1],
                                                box=self.box)
            self._s1 = ns_selection_1.search(self.u.atoms[self._s2],
                                             self.selection_distance).ix
        self.logger_debug("Size of selection 1: {0} atoms".format(
            len(self._s1)))

        if len(self._s1) == 0:
            logger.warning('Selection 1 "{0}" did not select any atoms.'
                           .format(str(self.selection1)[:80]))
            return

        if self.filter_first and len(self._s2):
            self.logger_debug('Size of selection 2 before filtering:'
                              ' {} atoms'.format(len(self._s2)))
            ns_selection_2 = AtomNeighborSearch(self.u.atoms[self._s2],
                                                box=self.box)
            self._s2 = ns_selection_2.search(self.u.atoms[self._s1],
                                             self.selection_distance).ix
        self.logger_debug('Size of selection 2: {0} atoms'.format(
            len(self._s2)))

        if len(self._s2) == 0:
            logger.warning('Selection 2 "{0}" did not select any atoms.'
                           .format(str(self.selection2)[:80]))
            return

        if self.selection1_type in ('donor', 'both'):
            self._s1_donors = self.u.atoms[self._s1].select_atoms(
                'name {0}'.format(' '.join(self.donors))).ix
            for atom_ix in self._s1_donors:
                self._update_donor_h(atom_ix, self._s1_h_donors,
                                     self._s1_donors_h)
            self.logger_debug("Selection 1 donors: {0}".format(
                len(self._s1_donors)))
            self.logger_debug("Selection 1 donor hydrogens: {0}".format(
                len(self._s1_h_donors)))
        if self.selection1_type in ('acceptor', 'both'):
            self._s1_acceptors = self.u.atoms[self._s1].select_atoms(
                'name {0}'.format(' '.join(self.acceptors))).ix
            self.logger_debug("Selection 1 acceptors: {0}".format(
                len(self._s1_acceptors)))

        if len(self._s2) == 0:
            return None
        if self.selection1_type in ('donor', 'both'):
            self._s2_acceptors = self.u.atoms[self._s2].select_atoms(
                'name {0}'.format(' '.join(self.acceptors))).ix
            self.logger_debug("Selection 2 acceptors: {0:d}".format(
                len(self._s2_acceptors)))
        if self.selection1_type in ('acceptor', 'both'):
            self._s2_donors = self.u.atoms[self._s2].select_atoms(
                'name {0}'.format(' '.join(self.donors))).ix
            for atom_ix in self._s2_donors:
                self._update_donor_h(atom_ix, self._s2_h_donors,
                                     self._s2_donors_h)
            self.logger_debug("Selection 2 donors: {0:d}".format(
                len(self._s2_donors)))
            self.logger_debug("Selection 2 donor hydrogens: {0:d}".format(
                len(self._s2_h_donors)))

    def _update_water_selection(self):
        self._water_donors = []
        self._water_h_donors = {}
        self._water_donors_h = defaultdict(list)
        self._water_acceptors = []

        self._water = self.u.select_atoms(self.water_selection).ix
        self.logger_debug('Size of water selection before filtering:'
                          ' {} atoms'.format(len(self._water)))
        if len(self._water) and self.filter_first:
            filtered_s1 = AtomNeighborSearch(self.u.atoms[self._water],
                                             box=self.box).search(
                self.u.atoms[self._s1], self.selection_distance)
            if filtered_s1:
                self._water = AtomNeighborSearch(filtered_s1,
                                                 box=self.box).search(
                    self.u.atoms[self._s2], self.selection_distance).ix

        self.logger_debug("Size of water selection: {0} atoms".format(
            len(self._water)))

        if len(self._water) == 0:
            logger.warning("Water selection '{0}' did not select any atoms."
                           .format(str(self.water_selection)[:80]))
        else:
            self._water_donors = self.u.atoms[self._water].select_atoms(
                'name {0}'.format(' '.join(self.donors))).ix
            for atom_ix in self._water_donors:
                self._update_donor_h(atom_ix, self._water_h_donors,
                                     self._water_donors_h)
            self.logger_debug("Water donors: {0}".format(
                len(self._water_donors)))
            self.logger_debug("Water donor hydrogens: {0}".format(
                len(self._water_h_donors)))
            self._water_acceptors = self.u.atoms[self._water].select_atoms(
                'name {0}'.format(' '.join(self.acceptors))).ix
            self.logger_debug("Water acceptors: {0}".format(
                len(self._water_acceptors)))

    def _get_bonded_hydrogens(self, atom):
        """Find hydrogens bonded within cutoff to `atom`.

        Hydrogens are detected by either name ("H*", "[123]H*") or type ("H");
        this is not fool-proof as the atom type is not always a character but
        the name pattern should catch most typical occurrences.

        The distance from `atom` is calculated for all hydrogens in the residue
        and only those within a cutoff are kept. The cutoff depends on the
        heavy atom (more precisely, on its element, which is taken as the first
        letter of its name ``atom.name[0]``) and is parameterized in
        :attr:`WaterBridgeAnalysis.r_cov`. If no match is found then the
        default of 1.5 Å is used.


        Parameters
        ----------
        atom : groups.Atom
             heavy atom

        Returns
        -------
        hydrogen_atoms : AtomGroup or []
            list of hydrogens (can be a
            :class:`~MDAnalysis.core.groups.AtomGroup`)
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
        # The distance for selection is defined as twice the maximum bond
        # length of an O-H bond (2A) plus order of water bridge times the
        # length of OH bond plus hydrogne bond distance
        self.selection_distance = (2 * 2 + self.order * (2 + self.distance))

        self.box = self.u.dimensions if self.pbc else None
        self._residue_dict = {}
        self._build_residue_dict(self.selection1)
        self._build_residue_dict(self.selection2)
        self._build_residue_dict(self.water_selection)

        self._update_selection()

        self.timesteps = []
        if len(self._s1) and len(self._s2):
            self._update_water_selection()
        else:
            logger.info("WaterBridgeAnalysis: "
                        "no atoms found in the selection.")

        logger.info("WaterBridgeAnalysis: initial checks passed.")

        logger.info("WaterBridgeAnalysis: starting")
        logger.debug("WaterBridgeAnalysis: donors    %r", self.donors)
        logger.debug("WaterBridgeAnalysis: acceptors %r", self.acceptors)
        logger.debug("WaterBridgeAnalysis: water bridge %r",
                     self.water_selection)

        if self.debug:
            logger.debug("Toggling debug to %r", self.debug)
        else:
            logger.debug("WaterBridgeAnalysis: For full step-by-step "
                         "debugging output use debug=True")

        logger.info("Starting analysis "
                    "(frame index start=%d stop=%d, step=%d)",
                    self.start, self.stop, self.step)

    def _donor2acceptor(self, donors, h_donors, acceptor):
        if len(donors) == 0 or len(acceptor) == 0:
            return []
        if self.distance_type != 'heavy':
            donors_idx = list(h_donors.keys())
        else:
            donors_idx = list(donors.keys())
        result = []
        # Code modified from p-j-smith
        pairs, distances = capped_distance(self.u.atoms[donors_idx].positions,
                                           self.u.atoms[acceptor].positions,
                                           max_cutoff=self.distance,
                                           box=self.box,
                                           return_distances=True)
        if self.distance_type == 'hydrogen':
            tmp_distances = distances
            tmp_donors = [h_donors[donors_idx[idx]] for idx in pairs[:, 0]]
            tmp_hydrogens = [donors_idx[idx] for idx in pairs[:, 0]]
            tmp_acceptors = [acceptor[idx] for idx in pairs[:, 1]]
        else:
            # To make sure that for the same index i, the donor (tmp_donors[i]),
            # hydrogen (tmp_hydrogens[i]), acceptor (tmp_acceptors[i]) matches the
            # distance (tmp_distances[i]).
            tmp_donors = []
            tmp_hydrogens = []
            tmp_acceptors = []
            tmp_distances = []
            for idx, distance in enumerate(distances):
                for h in donors[donors_idx[pairs[idx, 0]]]:
                    tmp_donors.append(donors_idx[pairs[idx, 0]])
                    tmp_hydrogens.append(h)
                    tmp_acceptors.append(acceptor[pairs[idx, 1]])
                    tmp_distances.append(distance)

        angles = np.rad2deg(
            calc_angles(
                self.u.atoms[tmp_donors].positions,
                self.u.atoms[tmp_hydrogens].positions,
                self.u.atoms[tmp_acceptors].positions,
                box=self.box
            )
        )
        hbond_indices = np.where(angles > self.angle)[0]
        for index in hbond_indices:
            h = tmp_hydrogens[index]
            d = tmp_donors[index]
            a = tmp_acceptors[index]
            result.append((h, d, a, self._expand_index(h),
                           self._expand_index(a),
                           tmp_distances[index], angles[index]))
        return result

    def _single_frame(self):
        self.timesteps.append(self._ts.time)
        self.box = self.u.dimensions if self.pbc else None

        if self.update_selection:
            self._update_selection()
        if len(self._s1) and len(self._s2):
            if self.update_water_selection:
                self._update_water_selection()
        else:
            self.results.network.append(defaultdict(dict))
            return

        selection_1 = []
        water_pool = defaultdict(list)
        next_round_water = set([])
        selection_2 = []

        if self.selection1_type in ('donor', 'both'):
            # check for direct hbond from s1 to s2
            self.logger_debug("Selection 1 Donors <-> Selection 2 Acceptors")
            results = self._donor2acceptor(
                self._s1_donors_h, self._s1_h_donors, self._s2_acceptors)
            for line in results:
                h_index, d_index, a_index, (h_resname, h_resid, h_name), \
                    (a_resname, a_resid, a_name), dist, angle = line
                water_pool[(a_resname, a_resid)] = None
                selection_1.append(
                    (h_index, d_index, a_index, None, dist, angle))
                selection_2.append((a_resname, a_resid))
            if self.order > 0:
                self.logger_debug("Selection 1 Donors <-> Water Acceptors")
                results = self._donor2acceptor(
                    self._s1_donors_h, self._s1_h_donors,
                    self._water_acceptors)
                for line in results:
                    h_index, d_index, a_index, (h_resname, h_resid, h_name), (
                        a_resname, a_resid, a_name), dist, angle = line
                    selection_1.append(
                        (h_index, d_index, a_index, None, dist, angle))

                self.logger_debug("Water Donors <-> Selection 2 Acceptors")
                results = self._donor2acceptor(
                    self._water_donors_h, self._water_h_donors,
                    self._s2_acceptors)
                for line in results:
                    h_index, d_index, a_index, (h_resname, h_resid, h_name), (
                        a_resname, a_resid, a_name), dist, angle = line
                    water_pool[(h_resname, h_resid)].append(
                        (h_index, d_index, a_index, None, dist, angle))
                    selection_2.append((a_resname, a_resid))

        if self.selection1_type in ('acceptor', 'both'):
            self.logger_debug("Selection 2 Donors <-> Selection 1 Acceptors")
            results = self._donor2acceptor(self._s2_donors_h,
                                           self._s2_h_donors,
                                           self._s1_acceptors)
            for line in results:
                h_index, d_index, a_index, (h_resname, h_resid, h_name), \
                    (a_resname, a_resid, a_name), dist, angle = line
                water_pool[(h_resname, h_resid)] = None
                selection_1.append(
                    (a_index, None, h_index, d_index, dist, angle))
                selection_2.append((h_resname, h_resid))

            if self.order > 0:
                self.logger_debug("Selection 2 Donors <-> Water Acceptors")
                results = self._donor2acceptor(
                    self._s2_donors_h, self._s2_h_donors,
                    self._water_acceptors)
                for line in results:
                    h_index, d_index, a_index, (h_resname, h_resid, h_name), (
                        a_resname, a_resid, a_name), dist, angle = line
                    water_pool[(a_resname, a_resid)].append(
                        (a_index, None, h_index, d_index, dist, angle))
                    selection_2.append((h_resname, h_resid))

                self.logger_debug("Selection 1 Acceptors <-> Water Donors")
                results = self._donor2acceptor(
                    self._water_donors_h, self._water_h_donors,
                    self._s1_acceptors)
                for line in results:
                    h_index, d_index, a_index, (h_resname, h_resid, h_name), (
                        a_resname, a_resid, a_name), dist, angle = line
                    selection_1.append(
                        (a_index, None, h_index, d_index, dist, angle))

        if self.order > 1:
            self.logger_debug("Water donor <-> Water Acceptors")
            results = self._donor2acceptor(self._water_donors_h,
                                           self._water_h_donors,
                                           self._water_acceptors)
            for line in results:
                h_index, d_index, a_index, (h_resname, h_resid, h_name), (
                    a_resname, a_resid, a_name), dist, angle = line
                water_pool[(a_resname, a_resid)].append(
                    (a_index, None, h_index, d_index, dist, angle))
                water_pool[(h_resname, h_resid)].append(
                    (h_index, d_index, a_index, None, dist, angle))

        # solve the connectivity network
        # The following code attempt to generate a water network which is
        # formed by the class dict.
        # Suppose we have a water bridge connection ARG1 to ASP3 via the two
        # hydrogen bonds.
        #     [0,1,('ARG',1,'O'),  ('SOL',2,'HW1'),  3.0,180],
        #     [2,3,('SOL',2,'HW2'),('ASP',3,'OD1'),  3.0,180],
        # The resulting network will be
        # {(0,1,('ARG',1,'O'),  ('SOL',2,'HW1'),  3.0,180):
        # {(2,3,('SOL',2,'HW2'),('ASP',3,'OD1'),  3.0,180): None}}
        # Where the key of the a dict will be all the hydrogen bonds starting
        # from this nodes.
        # The corresponding value of a certain key will be a dictionary whose
        # key will be all the hydrogen bonds from
        # the destination of in the key.
        # If the value of a certain key is None, which means it is reaching
        # selection 2.

        result = {'start': defaultdict(dict), 'water': defaultdict(dict)}

        def add_route(result, route):
            if len(route) == 1:
                result['start'][route[0]] = None
            else:
                # exclude the the selection which goes back to itself
                if (sorted(route[0][0:3:2]) == sorted(route[-1][0:3:2])):
                    return

                # selection 2 to water
                result['water'][route[-1]] = None
                # water to water
                for i in range(1, len(route) - 1):
                    result['water'][route[i]][route[i+1]] = \
                        result['water'][route[i+1]]
                # selection 1 to water
                result['start'][route[0]][route[1]] = result['water'][route[1]]

        def traverse_water_network(graph, node, end, route, maxdepth, result):
            if len(route) > self.order + 1:
                return
            else:
                if node in end:
                    # check if any duplication happens
                    heavy_atom = [line[3] or line[2] for line in route]
                    if len(heavy_atom) == len(set(heavy_atom)):
                        add_route(result, route)
                else:
                    for new_node in graph[node]:
                        new_route = route[:]
                        new_route.append(new_node)
                        new_node = self._expand_timeseries(
                            new_node, 'sele1_sele2')[3][:2]
                        traverse_water_network(graph, new_node, end, new_route,
                                               maxdepth, result)
        for s1 in selection_1:
            route = [s1, ]
            next_mol = self._expand_timeseries(s1, 'sele1_sele2')[3][:2]
            traverse_water_network(water_pool, next_mol, selection_2, route[:],
                                   self.order, result)

        self.results.network.append(result['start'])

    def _traverse_water_network(self, graph, current, analysis_func=None,
                                output=None, link_func=None, **kwargs):
        '''
        This function recursively traverses the water network
        self.results.network and finds the hydrogen bonds which connect the
        current atom to the next atom. The newly found hydrogen bond will be
        appended to the hydrogen bonds connecting the selection 1 to the
        current atom via link_func. When selection 2 is reached, the full list
        of hydrogen bonds connecting the selection 1 to selection 2 will be
        fed into analysis_func, which will then modify the output in place.

        :param graph: The connection network describes the connection between
        the atoms in the water network.
        :param current: The hydrogen bonds from selection 1 until now.
        :param analysis_func: The analysis function which is called to analysis
        the hydrogen bonds.
        :param output: where the result is stored.
        :param link_func: The new hydrogen bonds will be appended to current.
        :param kwargs: the keywords which are passed into the analysis_func.
        :return:
        '''
        if link_func is None:
            # If no link_func is provided, the default link_func will be used
            link_func = self._full_link

        if graph is None:
            # if selection 2 is reached
            if analysis_func is not None:
                # the result is analysed by analysis_func which will change the
                # output
                analysis_func(current, output, self.u, **kwargs)
        else:
            # make sure no loop can occur
            if len(current) <= self.order:
                for node in graph:
                    # the new hydrogen bond will be added to the existing bonds
                    new = link_func(current, node)
                    self._traverse_water_network(graph[node], new,
                                                 analysis_func, output,
                                                 link_func, **kwargs)

    def _expand_index(self, index):
        '''
        Expand the index into (resname, resid, name).
        '''
        atom = self.u.atoms[index]
        return (atom.resname, atom.resid, atom.name)

    def _expand_timeseries(self, entry, output_format=None):
        '''
        Expand the compact data format into the old timeseries form.
        The old is defined as the format for release up to 0.19.2.
        As is discussed in Issue #2177, the internal storage of the hydrogen
        bond information has been changed to the compact format.
        The function takes in the argument `output_format` to see which output
        format will be chosen.
        if `output_format` is not specified, the value will be taken from
        :attr:`output_format`.
        If `output_format` is 'sele1_sele2', the output will be the old water
        bridge analysis format::

          # donor from selection 1 to acceptor in selection 2
          [sele1_index, sele2_index,
           (sele1_resname, sele1_resid, sele1_name),
           (sele2_resname, sele2_resid, sele2_name), dist, angle]

        If `output_format` is 'donor_acceptor', the output will be the old
        hydrogen bond analysis format::

          # From donor to acceptor
          [donor_index, acceptor_index,
           (donor_resname, donor_resid, donor_name),
           (acceptor_resname, acceptor_resid, acceptor_name), dist, angle]
        '''
        output_format = output_format or self.output_format
        # Expand the compact entry into atom1, which is the first index in the
        # output and atom2, which is the second
        # entry.
        atom1, heavy_atom1, atom2, heavy_atom2, dist, angle = entry
        if output_format == 'sele1_sele2':
            # If the output format is the sele1_sele2, no change will be
            # executed
            atom1, atom2 = atom1, atom2
        elif output_format == 'donor_acceptor':
            # If the output format is donor_acceptor, use heavy atom position
            # to check which is donor and which is
            # acceptor
            if heavy_atom1 is None:
                # atom1 is hydrogen bond acceptor and thus, the position of
                # atom1 and atom2 are swapped.
                atom1, atom2 = atom2, atom1
            else:
                # atom1 is hydrogen bond donor, position not swapped.
                atom1, atom2 = atom1, atom2
        else:
            raise KeyError(
                "Only 'sele1_sele2' or 'donor_acceptor' are allowed as output "
                "format")

        return (atom1, atom2, self._expand_index(atom1),
                self._expand_index(atom2), dist, angle)

    def _generate_timeseries(self, output_format=None):
        r'''Time series of water bridges.

        The output is generated per frame as is explained in
        :ref:`wb_Analysis_Timeseries`. The format of output can be changed via
        the output_format selection. If ``output_format="sele1_sele2"``, the
        hydrogen bond forms a directional link from selection 1 to selection 2.
        If ``output_format="donor_acceptor"``, for each hydrogen bond, the
        donor is always written before the acceptor.

        Note
        ----
        To find an acceptor atom in :attr:`Universe.atoms` by
        *index* one would use ``u.atoms[acceptor_index]``.

        .. versionchanged 0.20.0
           The :attr:`WaterBridgeAnalysis.timeseries` has been updated where
           the donor and acceptor string has been changed to tuple
           (resname string, resid, name_string).


        '''
        output_format = output_format or self.output_format

        def analysis(current, output, *args, **kwargs):
            output = current

        timeseries = []
        for frame in self.results.network:
            new_frame = []
            self._traverse_water_network(frame, new_frame,
                                         analysis_func=analysis,
                                         output=new_frame,
                                         link_func=self._compact_link)
            timeseries.append([
                self._expand_timeseries(entry, output_format)
                for entry in new_frame])
        return timeseries


    def set_network(self, network):
        wmsg = ("The `set_network` method was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.network` instead")
        warnings.warn(wmsg, DeprecationWarning)
        self.results.network = network

    @classmethod
    def _full_link(self, output, node):
        '''
        A function used in _traverse_water_network to add the new hydrogen
        bond to the existing bonds.
        :param output: The existing hydrogen bonds from selection 1
        :param node: The new hydrogen bond
        :return: The hydrogen bonds from selection 1 with the new hydrogen
        bond added
        '''
        result = output[:]
        result.append(node)
        return result

    @classmethod
    def _compact_link(self, output, node):
        '''
        A function used in _traverse_water_network to add the new hydrogen
        bond to the existing bonds. In this form no new list is created and
        thus, one bridge will only appear once.
        :param output: The existing hydrogen bonds from selection 1
        :param node: The new hydrogen bond
        :return: The hydrogen bonds from selection 1 with the new hydrogen
        bond added
        '''
        output.append(node)
        return output

    def _count_by_type_analysis(self, current, output, *args, **kwargs):
        '''
        Generates the key for count_by_type analysis.
        :return:
        '''

        s1_index, to_index, s1, to_residue, dist, angle = \
            self._expand_timeseries(current[0])
        s1_resname, s1_resid, s1_name = s1
        from_index, s2_index, from_residue, s2, dist, angle = \
            self._expand_timeseries(current[-1])
        s2_resname, s2_resid, s2_name = s2
        key = (s1_index, s2_index,
               s1_resname, s1_resid, s1_name, s2_resname, s2_resid, s2_name)
        output[key] += 1

    def count_by_type(self, analysis_func=None, **kwargs):
        """Counts the frequency of water bridge of a specific type.

        If one atom *A* from *selection 1* is linked to atom *B* from
        *selection 2* through one or more bridging waters, an entity will be
        created and the proportion of time that this linkage exists in the
        whole simulation will be calculated.

        The identification of a specific type of water bridge can be modified
        by supplying the analysis_func function. See :ref:`wb_count_by_type`
        for detail.

        Returns
        -------
        counts : list
            Returns a :class:`list` containing atom indices for *A* and
            *B*, residue names, residue numbers, atom names (for both A and B)
            and the fraction of the total time during which the water bridge
            was detected. This method returns None if method
            :meth:`WaterBridgeAnalysis.run` was not executed first.


        """
        output = None
        if analysis_func is None:
            analysis_func = self._count_by_type_analysis
            output = 'combined'

        if self.results.network:
            length = len(self.results.network)
            result_dict = defaultdict(int)
            for frame in self.results.network:
                frame_dict = defaultdict(int)
                self._traverse_water_network(frame, [],
                                             analysis_func=analysis_func,
                                             output=frame_dict,
                                             link_func=self._full_link,
                                             **kwargs)
                for key, value in frame_dict.items():
                    result_dict[key] += frame_dict[key]

            if output == 'combined':
                result = [[i for i in key] for key in result_dict]
                [result[i].append(result_dict[key]/length)
                 for i, key in enumerate(result_dict)]
            else:
                result = [(key,
                           result_dict[key]/length) for key in result_dict]
            return result
        else:
            return None

    def _count_by_time_analysis(self, current, output, *args, **kwargs):
        s1_index, to_index, s1, to_residue, dist, angle = \
            self._expand_timeseries(current[0])
        s1_resname, s1_resid, s1_name = s1
        from_index, s2_index, from_residue, s2, dist, angle = \
            self._expand_timeseries(current[-1])
        s2_resname, s2_resid, s2_name = s2
        key = (s1_index, s2_index,
               s1_resname, s1_resid, s1_name, s2_resname, s2_resid, s2_name)
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
        if self.results.network:
            result = []
            for time, frame in zip(self.timesteps, self.results.network):
                result_dict = defaultdict(int)
                self._traverse_water_network(frame, [],
                                             analysis_func=analysis_func,
                                             output=result_dict,
                                             link_func=self._full_link,
                                             **kwargs)
                result.append((time,
                               sum([result_dict[key] for key in result_dict])))
            return result
        else:
            return None

    def _timesteps_by_type_analysis(self, current, output, *args, **kwargs):
        s1_index, to_index, s1, to_residue, dist, angle = \
            self._expand_timeseries(current[0])
        s1_resname, s1_resid, s1_name = s1
        from_index, s2_index, from_residue, s2, dist, angle = \
            self._expand_timeseries(current[-1])
        s2_resname, s2_resid, s2_name = s2
        key = (s1_index, s2_index, s1_resname, s1_resid, s1_name, s2_resname,
               s2_resid, s2_name)
        output[key].append(kwargs.pop('time'))

    def timesteps_by_type(self, analysis_func=None, **kwargs):
        """Frames during which each water bridges existed, sorted by each water
        bridges.

        Processes :attr:`WaterBridgeAnalysis.results.network` and returns a
        :class:`list` containing atom indices, residue names, residue
        numbers (from selection 1 and selection 2) and each timestep at which
        the water bridge was detected.

        Similar to :meth:`~WaterBridgeAnalysis.count_by_type` and
        :meth:`~WaterBridgeAnalysis.count_by_time`, the behavior can be
        adjusted by supplying an analysis_func.

        Returns
        -------
        data : list

        """
        output = None
        if analysis_func is None:
            analysis_func = self._timesteps_by_type_analysis
            output = 'combined'

        if self.results.network:
            result = defaultdict(list)
            if self.timesteps is None:
                timesteps = range(len(self.results.network))
            else:
                timesteps = self.timesteps
            for time, frame in zip(timesteps, self.results.network):
                self._traverse_water_network(frame, [],
                                             analysis_func=analysis_func,
                                             output=result,
                                             link_func=self._full_link,
                                             time=time, **kwargs)

            result_list = []
            for key, time_list in result.items():
                for time in time_list:
                    if output == 'combined':
                        key = list(key)
                        key.append(time)
                        result_list.append(key)
                    else:
                        result_list.append((key, time))
            return result_list
        else:
            return None

    def generate_table(self, output_format=None):
        """Generate a normalised table of the results.

        Parameters
        ----------
        output_format : {'sele1_sele2', 'donor_acceptor'}
            The output format of the `table` can be changed a fashion similar
            to :attr:`WaterBridgeAnalysis.results.timeseries` by changing the
            labels of the columns of the participating atoms.

        Returns
        -------
        table : numpy.rec.recarray
            A "tidy" table with one hydrogen bond per row, labeled according to
            `output_format` and containing information of atom_1, atom_2,
            distance, and angle.

        .. versionchanged:: 2.0.0
           Return the generated table (as well as storing it as :attr:`table`).

        .. deprecated:: 2.0.0
           In release 3.0.0, :meth:`generate_table()` will _only_ return the
           table and no longer store it in :attr:`table`.
        """
        output_format = output_format or self.output_format
        if self.results.network == []:
            msg = "No data computed, do run() first."
            warnings.warn(msg, category=MissingDataWarning)
            logger.warning(msg)
            return None

        if self.results.timeseries is not None \
          and output_format == self.output_format:
            timeseries = self.results.timeseries
        else:
            # Recompute timeseries with correct output format
            timeseries = self._generate_timeseries(output_format)

        num_records = np.sum([len(hframe) for hframe in timeseries])
        # build empty output table
        if output_format == 'sele1_sele2':
            dtype = [
                ("time", float),
                ("sele1_index", int), ("sele2_index", int),
                ("sele1_resnm", "|U4"), ("sele1_resid", int),
                ("sele1_atom", "|U4"),
                ("sele2_resnm", "|U4"), ("sele2_resid", int),
                ("sele2_atom", "|U4"),
                ("distance", float), ("angle", float)]
        elif output_format == 'donor_acceptor':
            dtype = [
                ("time", float),
                ("donor_index", int), ("acceptor_index", int),
                ("donor_resnm", "|U4"), ("donor_resid", int),
                ("donor_atom", "|U4"),
                ("acceptor_resnm", "|U4"), ("acceptor_resid", int),
                ("acceptor_atom", "|U4"),
                ("distance", float), ("angle", float)]

        # according to Lukas' notes below, using a recarray at this stage is
        # ineffective and speedups of ~x10 can be achieved by filling a
        # standard array, like this:
        out = np.empty((num_records,), dtype=dtype)
        cursor = 0  # current row
        for t, hframe in zip(self.timesteps, timeseries):
            for (donor_index, acceptor_index, donor,
                 acceptor, distance, angle) in hframe:
                # donor|acceptor = (resname, resid, atomid)
                out[cursor] = (t, donor_index, acceptor_index) + \
                    donor + acceptor + (distance, angle)
                cursor += 1
        assert cursor == num_records, \
            "Internal Error: Not all wb records stored"
        table = out.view(np.rec.recarray)
        logger.debug(
            "WBridge: Stored results as table with %(num_records)d entries.",
            vars())
        self.table = table

        return table

    def _conclude(self):
        self.results.timeseries = self._generate_timeseries()

    @property
    def network(self):
        wmsg = ("The `network` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.network` instead")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.network

    @property
    def timeseries(self):
        wmsg = ("The `timeseries` attribute was deprecated in MDAnalysis "
                "2.0.0 and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.timeseries` instead")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.timeseries
