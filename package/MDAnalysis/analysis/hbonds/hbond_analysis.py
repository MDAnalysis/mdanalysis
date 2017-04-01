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

# Hydrogen Bonding Analysis
"""
Hydrogen Bond analysis --- :mod:`MDAnalysis.analysis.hbonds.hbond_analysis`
===========================================================================

:Author: David Caplan, Lukas Grossar, Oliver Beckstein
:Year: 2010-2012
:Copyright: GNU Public License v3


Given a :class:`~MDAnalysis.core.universe.Universe` (simulation
trajectory with 1 or more frames) measure all hydrogen bonds for each
frame between selections 1 and 2.

The :class:`HydrogenBondAnalysis` class is modeled after the `VMD
HBONDS plugin`_.

.. _`VMD HBONDS plugin`: http://www.ks.uiuc.edu/Research/vmd/plugins/hbonds/

Options:
  - *update_selections* (``True``): update selections at each frame?
  - *selection_1_type* ("both"): selection 1 is the: "donor", "acceptor", "both"
  - donor-acceptor *distance* (Å): 3.0
  - Angle *cutoff* (degrees): 120.0
  - *forcefield* to switch between default values for different force fields
  - *donors* and *acceptors* atom types (to add additional atom names)

.. _Analysis Output:

Output
------

The results are
  - the **identities** of donor and acceptor heavy-atoms,
  - the **distance** between the heavy atom acceptor atom and the hydrogen atom
    that is bonded to the heavy atom donor,
  - the **angle** donor-hydrogen-acceptor angle (180º is linear).

Hydrogen bond data are returned per frame, which is stored in
:attr:`HydrogenBondAnalysis.timeseries` (In the following description, ``#``
indicates comments that are not part of the output.)::

    results = [
        [ # frame 1
           [ # hbond 1
              <donor index (1-based)>, <acceptor index (1-based)>, <donor index (0-based)>,
              <acceptor index (0-based)>, <donor string>, <acceptor string>,
              <distance>, <angle>
           ],
           [ # hbond 2
              <donor index (1-based)>, <acceptor index (1-based)>, <donor index (0-based)>,
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

.. Note::

   For historic reasons, the *donor index* and *acceptor index* are a 1-based
   indices. To get the :attr:`Atom.index` (the 0-based index typically used in
   MDAnalysis simply subtract 1. For instance, to find an atom in
   :attr:`Universe.atoms` by *index* from the output one would use
   ``u.atoms[index-1]``.
   .. deprecated:: 0.15.0
   This feature is being deprecated in favor of zero-based indices and is targeted
   for removal in 0.16.0.


Using the :meth:`HydrogenBondAnalysis.generate_table` method one can reformat
the results as a flat "normalised" table that is easier to import into a
database for further processing. :meth:`HydrogenBondAnalysis.save_table` saves
the table to a pickled file. The table itself is a :class:`numpy.recarray`.


Detection of hydrogen bonds
---------------------------

Hydrogen bonds are recorded based on a geometric criterion:

1. The distance between acceptor and hydrogen is less than or equal to
   *distance* (default is 3 Å).

2. The angle between donor-hydrogen-acceptor is greater than or equal to
   *angle* (default is 120º).

The cut-off values *angle* and *distance* can be set as keywords to
:class:`HydrogenBondAnalysis`.

Donor and acceptor heavy atoms are detected from atom names. The current
defaults are appropriate for the CHARMM27 and GLYCAM06 force fields as defined
in Table `Default atom names for hydrogen bonding analysis`_.

Hydrogen atoms bonded to a donor are searched with one of two algorithms,
selected with the *detect_hydrogens* keyword.

*distance*

   Searches for all hydrogens (name "H*" or name "[123]H" or type "H") in the
   same residue as the donor atom within a cut-off distance of 1.2 Å.

*heuristic*

   Looks at the next three atoms in the list of atoms following the donor and
   selects any atom whose name matches (name "H*" or name "[123]H"). For

The *distance* search is more rigorous but slower and is set as the
default. Until release 0.7.6, only the heuristic search was implemented.

.. versionchanged:: 0.7.6
   Distance search added (see
   :meth:`HydrogenBondAnalysis._get_bonded_hydrogens_dist`) and heuristic
   search improved (:meth:`HydrogenBondAnalysis._get_bonded_hydrogens_list`)

.. _Default atom names for hydrogen bonding analysis:

.. table:: Default heavy atom names for CHARMM27 force field.

   =========== ==============  =========== ====================================
   group       donor           acceptor    comments
   =========== ==============  =========== ====================================
   main chain  N               O
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
   MET                         SD          see e.g. [Gregoret1991]_
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
OPLS/AA (tested in Gromacs). Residue names in the table are for information
only and are not taken into account when determining acceptors and donors.
This can potentially lead to some ambiguity in the assignment of
donors/acceptors for residues such as histidine or cytosine.

For more information about the naming convention in GLYCAM06 have a look at the
`Carbohydrate Naming Convention in Glycam`.

.. _`Carbohydrate Naming Convention in Glycam`:
   http://glycam.ccrc.uga.edu/documents/FutureNomenclature.htm

The lists of donor and acceptor names can be extended by providing lists of
atom names in the *donors* and *acceptors* keywords to
:class:`HydrogenBondAnalysis`. If the lists are entirely inappropriate
(e.g. when analysing simulations done with a force field that uses very
different atom names) then one should either use the value "other" for *forcefield*
to set no default values, or derive a new class and set the default list oneself::

 class HydrogenBondAnalysis_OtherFF(HydrogenBondAnalysis):
       DEFAULT_DONORS = {"OtherFF": tuple(set([...]))}
       DEFAULT_ACCEPTORS = {"OtherFF": tuple(set([...]))}

Then simply use the new class instead of the parent class and call it with
*forcefield* = "OtherFF". Please also consider to contribute the list of heavy
atom names to MDAnalysis.

.. rubric:: References

.. [Gregoret1991] L.M. Gregoret, S.D. Rader, R.J. Fletterick, and
   F.E. Cohen. Hydrogen bonds involving sulfur atoms in proteins. Proteins,
   9(2):99–107, 1991. `10.1002/prot.340090204`_.

.. _`10.1002/prot.340090204`: http://dx.doi.org/10.1002/prot.340090204


How to perform HydrogenBondAnalysis
-----------------------------------

All protein-water hydrogen bonds can be analysed with ::

  import MDAnalysis
  import MDAnalysis.analysis.hbonds

  u = MDAnalysis.Universe('topology', 'trajectory')
  h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u, 'protein', distance=3.0, angle=120.0)
  h.run()

The results are stored as the attribute
:attr:`HydrogenBondAnalysis.timeseries`; see :ref:`Analysis Output`
for the format and further options.

.. Note::

   Due to the way :class:`HydrogenBondAnalysis` is implemented, it is
   more efficient to have the second selection (*selection2*) be the
   *larger* group, e.g. the water when looking at water-protein
   H-bonds or the whole protein when looking at ligand-protein
   interactions.

.. Note::

   The topology supplied and the trajectory must reflect the same total number
   of atoms.

.. TODO: how to analyse the ouput and notes on selection updating


Classes
-------

.. autoclass:: HydrogenBondAnalysis
   :members:

   .. attribute:: timesteps

      List of the times of each timestep. This can be used together with
      :attr:`~HydrogenBondAnalysis.timeseries` to find the specific time point
      of a hydrogen bond existence, or see :attr:`~HydrogenBondAnalysis.table`.

   .. attribute:: timeseries

      Results of the hydrogen bond analysis, stored for each frame. In
      the following description, # indicates comments that are not
      part of the output::

        results = [
            [ # frame 1
               [ # hbond 1
                  <donor index (1-based)>, <acceptor index (1-based)>, <donor index (0-based)>,
                  <acceptor index (0-based)>, <donor string>, <acceptor string>,
                  <distance>, <angle>
               ],
               [ # hbond 2
                  <donor index (1-based)>, <acceptor index (1-based)>, <donor index (0-based)>,
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
      The time of each step is not stored with each hydrogen bond frame but in
      :attr:`~HydrogenBondAnalysis.timesteps`.

      .. Note::

         The *index* is a 1-based index. To get the :attr:`Atom.index` (the
         0-based index typically used in MDAnalysis simply subtract 1. For
         instance, to find an atom in :attr:`Universe.atoms` by *index* one
         would use ``u.atoms[index-1]``.



   .. attribute:: table

      A normalised table of the data in
      :attr:`HydrogenBondAnalysis.timeseries`, generated by
      :meth:`HydrogenBondAnalysis.generate_table`. It is a
      :class:`numpy.recarray` with the following columns:

          0. "time"
          1. "donor_idx"
          2. "acceptor_idx"
          3. "donor_index"
          4. "acceptor_index"
          5. "donor_resnm"
          6. "donor_resid"
          7. "donor_atom"
          8. "acceptor_resnm"
          9. "acceptor_resid"
          10. "acceptor_atom"
          11. "distance"
          12. "angle"

      It takes up more space than
      :attr:`~HydrogenBondAnalysis.timeseries` but it is easier to
      analyze and to import into databases (e.g. using recsql_).

      .. Note::

         The *index* is a 1-based index. To get the :attr:`Atom.index` (the
         0-based index typically used in MDAnalysis simply subtract 1. For
         instance, to find an atom in :attr:`Universe.atoms` by *index* one
         would use ``u.atoms[idx_zero]``. The 1-based index is deprecated and
         targeted for removal in 0.16.0



   .. automethod:: _get_bonded_hydrogens

   .. automethod:: _get_bonded_hydrogens_dist

   .. automethod:: _get_bonded_hydrogens_list

   .. deprecated:: 0.15.0
        The donor and acceptor indices being 1-based is deprecated in favor of
        a zero-based index. This can be accessed by "donor_index" or
        "acceptor_index" removal of the 1-based indices is targeted
        for version 0.16.0

"""
import six
from six.moves import range, zip, map, cPickle

from collections import defaultdict
import numpy as np
import warnings
import logging

from MDAnalysis import MissingDataWarning, NoDataError, SelectionError, SelectionWarning
from MDAnalysis.lib.util import parse_residue
from MDAnalysis.lib.mdamath import norm, angle
from MDAnalysis.lib.log import ProgressMeter, _set_verbose
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch


logger = logging.getLogger('MDAnalysis.analysis.hbonds')


class HydrogenBondAnalysis(object):
    """Perform a hydrogen bond analysis

    The analysis of the trajectory is performed with the
    :meth:`HydrogenBondAnalysis.run` method. The result is stored in
    :attr:`HydrogenBondAnalysis.timeseries`. See
    :meth:`~HydrogenBondAnalysis.run` for the format.

    The default atom names are taken from the CHARMM 27 force field files, which
    will also work for e.g. OPLS/AA in Gromacs, and GLYCAM06.

    *Donors* (associated hydrogens are deduced from topology)
      *CHARMM 27*
        N of the main chain, water OH2/OW, ARG NE/NH1/NH2, ASN ND2, HIS ND1/NE2,
        SER OG, TYR OH, CYS SG, THR OG1, GLN NE2, LYS NZ, TRP NE1
      *GLYCAM06*
        N,NT,N3,OH,OW

    *Acceptors*
      *CHARMM 27*
        O of the main chain, water OH2/OW, ASN OD1, ASP OD1/OD2, CYH SG, GLN OE1,
        GLU OE1/OE2, HIS ND1/NE2, MET SD, SER OG, THR OG1, TYR OH
      *GLYCAM06*
        N,NT,O,O2,OH,OS,OW,OY,P,S,SM

    .. SeeAlso:: Table :ref:`Default atom names for hydrogen bonding analysis`

    .. versionchanged:: 0.7.6
       DEFAULT_DONORS/ACCEPTORS is now embedded in a dict to switch between
       default values for different force fields.
    """

    # use tuple(set()) here so that one can just copy&paste names from the
    # table; set() takes care for removing duplicates. At the end the
    # DEFAULT_DONORS and DEFAULT_ACCEPTORS should simply be tuples.

    # : default heavy atom names whose hydrogens are treated as *donors*
    # : (see :ref:`Default atom names for hydrogen bonding analysis`)
    #: Use the keyword *donors* to add a list of additional donor names.
    DEFAULT_DONORS = {
        'CHARMM27': tuple(set([
            'N', 'OH2', 'OW', 'NE', 'NH1', 'NH2', 'ND2', 'SG', 'NE2', 'ND1', 'NZ', 'OG', 'OG1', 'NE1', 'OH'])),
        'GLYCAM06': tuple(set(['N', 'NT', 'N3', 'OH', 'OW'])),
        'other': tuple(set([]))}

    #: default atom names that are treated as hydrogen *acceptors*
    #: (see :ref:`Default atom names for hydrogen bonding analysis`)
    #: Use the keyword *acceptors* to add a list of additional acceptor names.
    DEFAULT_ACCEPTORS = {
        'CHARMM27': tuple(set([
            'O', 'OH2', 'OW', 'OD1', 'OD2', 'SG', 'OE1', 'OE1', 'OE2', 'ND1', 'NE2', 'SD', 'OG', 'OG1', 'OH'])),
        'GLYCAM06': tuple(set(['N', 'NT', 'O', 'O2', 'OH', 'OS', 'OW', 'OY', 'SM'])),
        'other': tuple(set([]))}

    #: A :class:`collections.defaultdict` of covalent radii of common donors
    #: (used in :meth`_get_bonded_hydrogens_list` to check if a hydrogen is
    #: sufficiently close to its donor heavy atom). Values are stored for
    #: N, O, P, and S. Any other heavy atoms are assumed to have hydrogens
    #: covalently bound at a maximum distance of 1.5 Å.
    r_cov = defaultdict(lambda: 1.5,  # default value
                        N=1.31, O=1.31, P=1.58, S=1.55)

    def __init__(self, universe, selection1='protein', selection2='all', selection1_type='both',
                 update_selection1=True, update_selection2=True, filter_first=True, distance_type='hydrogen',
                 distance=3.0, angle=120.0,
                 forcefield='CHARMM27', donors=None, acceptors=None,
                 start=None, stop=None, step=None,
                 debug=None, detect_hydrogens='distance', verbose=None):
        """Set up calculation of hydrogen bonds between two selections in a universe.

        The timeseries is accessible as the attribute :attr:`HydrogenBondAnalysis.timeseries`.

        Some initial checks are performed. If there are no atoms selected by
        *selection1* or *selection2* or if no donor hydrogens or acceptor atoms
        are found then a :exc:`SelectionError` is raised for any selection that
        does *not* update (*update_selection1* and *update_selection2*
        keywords). For selections that are set to update, only a warning is
        logged because it is assumed that the selection might contain atoms at
        a later frame (e.g. for distance based selections).

        If no hydrogen bonds are detected or if the initial check fails, look
        at the log output (enable with :func:`MDAnalysis.start_logging` and set
        *verbose* = ``True``). It is likely that the default names for donors
        and acceptors are not suitable (especially for non-standard
        ligands). In this case, either change the *forcefield* or use
        customized *donors* and/or *acceptors*.

        .. Note::

           In order to speed up processing, atoms are filtered by a coarse
           distance criterion before a detailed hydrogen bonding analysis is
           performed (*filter_first* = ``True``). If one of your selections is
           e.g. the solvent then *update_selection1* (or *update_selection2*)
           must also be ``True`` so that the list of candidate atoms is updated
           at each step: this is now the default.

           If your selections will essentially remain the same for all time
           steps (i.e. residues are not moving farther than 3 x *distance*), for
           instance, if no water or large conformational changes are involved
           or if the optimization is disabled (*filter_first* = ``False``) then
           you can improve performance by setting the *update_selection*
           keywords to ``False``.

        Parameters
        ----------
        universe : universe
            Universe object
        selection1 : str (optional)
            Selection string for first selection ['protein']
        selection2 : str (optional)
            Selection string for second selection ['all']
        selection1_type : str (optional)
            Selection 1 can be 'donor', 'acceptor' or 'both'. Note that the
            value for *selection1_type* automatically determines how
            *selection2* handles donors and acceptors: If *selection1* contains
            'both' then *selection2* will also contain *both*. If *selection1*
            is set to 'donor' then *selection2* is 'acceptor' (and vice versa).
            ['both'].
        update_selection1 : bool (optional)
            Update selection 1 at each frame? [``False``]
        update_selection2 : bool (optional)
            Update selection 2 at each frame? [``False``]
        filter_first : bool (optional)
            Filter selection 2 first to only atoms 3*distance away [``True``]
        distance : float (optional)
            Distance cutoff for hydrogen bonds; only interactions with a H-A distance
            <= *distance* (and the appropriate D-H-A angle, see *angle*) are
            recorded. (Note: *distance_type* can change this to the D-A distance.) [3.0]
        angle : float (optional)
            Angle cutoff for hydrogen bonds; an ideal H-bond has an angle of
            180º.  A hydrogen bond is only recorded if the D-H-A angle is
            >=  *angle*. The default of 120º also finds fairly non-specific
            hydrogen interactions and a possibly better value is 150º. [120.0]
        forcefield : {"CHARMM27", "GLYCAM06", "other"} (optional)
            Name of the forcefield used. Switches between different
            :attr:`~HydrogenBondAnalysis.DEFAULT_DONORS` and
            :attr:`~HydrogenBondAnalysis.DEFAULT_ACCEPTORS` values.
            Available values: "CHARMM27", "GLYCAM06", "other" ["CHARMM27"]
        donors : sequence (optional)
            Extra H donor atom types (in addition to those in
            :attr:`~HydrogenBondAnalysis.DEFAULT_DONORS`), must be a sequence.
        acceptors : sequence (optional)
            Extra H acceptor atom types (in addition to those in
            :attr:`~HydrogenBondAnalysis.DEFAULT_ACCEPTORS`), must be a sequence.
        start : int (optional)
            starting frame-index for analysis, ``None`` is the first one, 0.
            *start* and *stop* are 0-based frame indices and are used to slice
            the trajectory (if supported) [``None``]
        stop : int (optional)
            last trajectory frame for analysis, ``None`` is the last one [``None``]
        step : int (optional)
            read every *step* between *start* and *stop*, ``None`` selects 1.
            Note that not all trajectory reader from 1 [``None``]
        debug : bool (optional)
            If set to ``True`` enables per-frame debug logging. This is disabled
            by default because it generates a very large amount of output in
            the log file. (Note that a logger must have been started to see
            the output, e.g. using :func:`MDAnalysis.start_logging`.)
        detect_hydrogens : {"distance", "heuristic"} (optional)
            Determine the algorithm to find hydrogens connected to donor
            atoms. Can be "distance" (default; finds all hydrogens in the
            donor's residue within a cutoff of the donor) or "heuristic"
            (looks for the next few atoms in the atom list). "distance" should
            always give the correct answer but "heuristic" is faster,
            especially when the donor list is updated each
            for each frame. ["distance"]
        distance_type : {"hydrogen", "heavy"} (optional)
            Measure hydrogen bond lengths between donor and acceptor heavy
            attoms ("heavy") or between donor hydrogen and acceptor heavy
            atom ("hydrogen"). If using "heavy" then one should set the *distance*
            cutoff to a higher value such as 3.5 Å. ["hydrogen"]

        Raises
        ------
        :exc:`SelectionError`
            is raised for each static selection without the required
            donors and/or acceptors.

        .. versionchanged:: 0.7.6
           New *verbose* keyword (and per-frame debug logging disabled by
           default).

           New *detect_hydrogens* keyword to switch between two different
           algorithms to detect hydrogens bonded to donor. "distance" is a new,
           rigorous distance search within the residue of the donor atom,
           "heuristic" is the previous list scan (improved with an additional
           distance check).

           New *forcefield* keyword to switch between different values of
           DEFAULT_DONORS/ACCEPTORS to accomodate different force fields.
           Also has an option "other" for no default values.

        .. versionchanged:: 0.8
           The new default for *update_selection1* and *update_selection2* is now
           ``True`` (see `Issue 138`_). Set to ``False`` if your selections only
           need to be determined once (will increase performance).

        .. versionchanged:: 0.9.0
           New keyword *distance_type* to select between calculation between
           heavy atoms or hydrogen-acceptor. It defaults to the previous
           behavior (i.e. "hydrogen").

        .. versionchanged:: 0.11.0
           Initial checks for selections that potentially raise :exc:`SelectionError`.

        .. deprecated:: 0.16
           The *verbose* keyword argument is replaced by *debug*. Note that the
           *verbose* keyword argument is now comsistently used to toggle
           progress meters throuthout the library.

        .. _`Issue 138`: https://github.com/MDAnalysis/mdanalysis/issues/138
        """
        warnings.warn(
             "The donor and acceptor indices being 1-based is deprecated in favor"
             " of a zero-based index. These can be accessed by 'donor_index' or"
             " 'acceptor_index', removal of the 1-based indices is targeted for"
             " version 0.16.0", category=DeprecationWarning)

        self._get_bonded_hydrogens_algorithms = {
            "distance": self._get_bonded_hydrogens_dist,  # 0.7.6 default
            "heuristic": self._get_bonded_hydrogens_list,  # pre 0.7.6
        }
        if not detect_hydrogens in self._get_bonded_hydrogens_algorithms:
            raise ValueError("detect_hydrogens must be one of {0!r}".format(
                             self._get_bonded_hydrogens_algorithms.keys()))
        self.detect_hydrogens = detect_hydrogens

        self.u = universe
        self.selection1 = selection1
        self.selection2 = selection2
        self.selection1_type = selection1_type
        self.update_selection1 = update_selection1
        self.update_selection2 = update_selection2
        self.filter_first = filter_first
        self.distance = distance
        self.distance_type = distance_type  # note: everything except 'heavy' will give the default behavior
        self.angle = angle
        self.traj_slice = slice(start, stop, step)

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

        self.timeseries = None  # final result
        self.timesteps = None  # time for each frame

        self.table = None  # placeholder for output table

        self.debug = True  # always enable debug output for initial selection update
        self._update_selection_1()
        self._update_selection_2()
        # per-frame debugging output?
        # This line must be changed at the end of the deprecation period for
        # the *quiet* keyword argument. Then it must become:
        # self.debug = debug
        # In the signature, *verbose* must be removed and the default value
        # for *debug* must be set to False.
        # See the docstring for lib.log._set_verbose, the pull request #1150,
        # and the issue #903.
        self.debug = _set_verbose(debug, verbose, default=False,
                                  was='verbose', now='debug')

        self._log_parameters()

        if self.selection1_type == 'donor':
            self._sanity_check(1, 'donors')
            self._sanity_check(2, 'acceptors')
        elif self.selection1_type == 'acceptor':
            self._sanity_check(1, 'acceptors')
            self._sanity_check(2, 'donors')
        else:  # both
            self._sanity_check(1, 'donors')
            self._sanity_check(1, 'acceptors')
            self._sanity_check(2, 'acceptors')
            self._sanity_check(2, 'donors')
        logger.info("HBond analysis: initial checks passed.")


    def _sanity_check(self, selection, htype):
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
        assert selection in (1, 2)
        assert htype in ("donors", "acceptors")
        # horrible data organization:  _s1_donors, _s2_acceptors, etc, update_selection1, ...
        atoms = getattr(self, "_s{0}_{1}".format(selection, htype))
        update = getattr(self, "update_selection{0}".format(selection))
        if not atoms:
            errmsg = "No {1} found in selection {0}. " \
                "You might have to specify a custom '{1}' keyword.".format(
                selection, htype)
            if not update:
                logger.error(errmsg)
                raise SelectionError(errmsg)
            else:
                errmsg += " Selection will update so continuing with fingers crossed."
                warnings.warn(errmsg, category=SelectionWarning)
                logger.warn(errmsg)

    def _log_parameters(self):
        """Log important parameters to the logfile."""
        logger.info("HBond analysis: selection1 = %r (update: %r)", self.selection1, self.update_selection1)
        logger.info("HBond analysis: selection2 = %r (update: %r)", self.selection2, self.update_selection2)
        logger.info("HBond analysis: criterion: donor %s atom and acceptor atom distance <= %.3f A", self.distance_type,
                    self.distance)
        logger.info("HBond analysis: criterion: angle D-H-A >= %.3f degrees", self.angle)
        logger.info("HBond analysis: force field %s to guess donor and acceptor names", self.forcefield)
        logger.info("HBond analysis: bonded hydrogen detection algorithm: %r", self.detect_hydrogens)

    def _get_bonded_hydrogens(self, atom, **kwargs):
        """Find hydrogens bonded to *atom*.

        This method is typically not called by a user but it is documented to
        facilitate understanding of the internals of
        :class:`HydrogenBondAnalysis`.

        :Returns: list of hydrogens (can be a
                  :class:`~MDAnalysis.core.groups.AtomGroup`) or empty list
                  ``[]`` if none were found.

        .. SeeAlso::

           :meth:`_get_bonded_hydrogens_dist` and :meth:`_get_bonded_hydrogens_list`

        .. versionchanged:: 0.7.6
           Can switch algorithm by using the *detect_hydrogens* keyword to the
           constructor. *kwargs* can be used to supply arguments for algorithm.
        """
        return self._get_bonded_hydrogens_algorithms[self.detect_hydrogens](atom, **kwargs)

    def _get_bonded_hydrogens_dist(self, atom):
        """Find hydrogens bonded within *cutoff* to *atom*.

        * hydrogens are detected by either name ("H*", "[123]H*") or type
          ("H"); this is not fool-proof as the atom type is not always a
          character but the name pattern should catch most typical occurrences.

        * The distance from *atom* is calculated for all hydrogens in the
          residue and only those within a cutoff are kept. The cutoff depends
          on the heavy atom (more precisely, on its element, which is taken as
          the first letter of its name ``atom.name[0]``) and is parameterized
          in :attr:`HydrogenBondAnalysis.r_cov`. If no match is found then the
          default of 1.5 Å is used.

        The performance of this implementation could be improved once the
        topology always contains bonded information; it currently uses the
        selection parser with an "around" selection.

        .. versionadded:: 0.7.6
        """
        try:
            return atom.residue.atoms.select_atoms(
                "(name H* 1H* 2H* 3H* or type H) and around {0:f} name {1!s}"
                "".format(self.r_cov[atom.name[0]], atom.name))
        except NoDataError:
            return []

    def _get_bonded_hydrogens_list(self, atom, **kwargs):
        """Find "bonded" hydrogens to the donor *atom*.

        At the moment this relies on the **assumption** that the
        hydrogens are listed directly after the heavy atom in the
        topology. If this is not the case then this function will
        fail.

        Hydrogens are detected by name ``H*``, ``[123]H*`` and they have to be
        within a maximum distance from the heavy atom. The cutoff distance
        depends on the heavy atom and is parameterized in
        :attr:`HydrogenBondAnalysis.r_cov`.

        .. versionchanged:: 0.7.6

           Added detection of ``[123]H`` and additional check that a
           selected hydrogen is bonded to the donor atom (i.e. its
           distance to the donor is less than the covalent radius
           stored in :attr:`HydrogenBondAnalysis.r_cov` or the default
           1.5 Å).

           Changed name to
           :meth:`~HydrogenBondAnalysis._get_bonded_hydrogens_list`
           and added *kwargs* so that it can be used instead of
           :meth:`~HydrogenBondAnalysis._get_bonded_hydrogens_dist`.

        """
        warnings.warn("_get_bonded_hydrogens_list() does not always find "
                      "all hydrogens; detect_hydrogens='distance' is safer.",
                      category=DeprecationWarning)
        try:
            hydrogens = [
                a for a in self.u.atoms[atom.index + 1:atom.index + 4]
                if a.name.startswith(('H', '1H', '2H', '3H')) \
                    and self.calc_eucl_distance(atom, a) < self.r_cov[atom.name[0]]]
        except IndexError:
            hydrogens = []  # weird corner case that atom is the last one in universe
        return hydrogens

    def _update_selection_1(self):
        self._s1 = self.u.select_atoms(self.selection1)
        self.logger_debug("Size of selection 1: {0} atoms".format(len(self._s1)))
        if not self._s1:
            logger.warn("Selection 1 '{0}' did not select any atoms.".format(str(self.selection1)[:80]))
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

    def _update_selection_2(self):
        self._s2 = self.u.select_atoms(self.selection2)
        if self.filter_first and self._s2:
            self.logger_debug('Size of selection 2 before filtering:'
                              ' {} atoms'.format(len(self._s2)))
            ns_selection_2 = AtomNeighborSearch(self._s2)
            self._s2 = ns_selection_2.search(self._s1, 3. * self.distance)
        self.logger_debug('Size of selection 2: {0} atoms'.format(len(self._s2)))
        if not self._s2:
            logger.warn('Selection 2 "{0}" did not select any atoms.'.format(
                str(self.selection2)[:80]))
        self._s2_donors = {}
        self._s2_donors_h = {}
        self._s2_acceptors = {}
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

    def logger_debug(self, *args):
        if self.debug:
            logger.debug(*args)

    def run(self, **kwargs):
        """Analyze trajectory and produce timeseries.

        Stores the hydrogen bond data per frame as
        :attr:`HydrogenBondAnalysis.timeseries` (see there for output
        format).

        The method accepts a number of keywords, amongst them *verbose*
        (default ``True``), which toggles the porgress output (see
        :class:`~MDAnalysis.lib.log.ProgressMeter`) and *debug* which can
        be used to change the debug value provided to the class constructor.

        .. SeeAlso:: :meth:`HydrogenBondAnalysis.generate_table` for processing
                     the data into a different format.

        .. versionchanged:: 0.7.6
           Results are not returned, only stored in
           :attr:`~HydrogenBondAnalysis.timeseries` and duplicate hydrogen bonds
           are removed from output (can be suppressed with *remove_duplicates* =
           ``False``)

        .. versionchanged:: 0.11.0
           Accept *quiet* keyword. Analysis will now proceed through frames even if
           no donors or acceptors were found in a particular frame.

        .. deprecated:: 0.15.0
           The donor and acceptor indices being 1-based is deprecated in favor
           of a zero-based index. This can be accessed by "donor_index" or
           "acceptor_index" removal of the 1-based indices is targeted
           for version 0.16.0

        .. deprecated:: 0.16
           The *quiet* keyword argument is deprecated in favor of the *verbose*
           one. Previous use of *verbose* now corresponds to the new keyword
           argument *debug*.

        """
        logger.info("HBond analysis: starting")
        logger.debug("HBond analysis: donors    %r", self.donors)
        logger.debug("HBond analysis: acceptors %r", self.acceptors)

        remove_duplicates = kwargs.pop('remove_duplicates', True)  # False: old behaviour
        if not remove_duplicates:
            logger.warn("Hidden feature remove_duplicates=False activated: you will probably get duplicate H-bonds.")

        debug = kwargs.pop('debug', None)
        if debug is not None and debug != self.debug:
            self.debug = debug
            logger.debug("Toggling debug to %r", self.debug)
        if not self.debug:
            logger.debug("HBond analysis: For full step-by-step debugging output use debug=True")

        self.timeseries = []
        self.timesteps = []

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
                           format="HBonds frame {current_step:5d}: {step:5d}/{numsteps} [{percentage:5.1f}%]\r",
                           verbose=verbose)

        try:
            self.u.trajectory.time

            def _get_timestep():
                return self.u.trajectory.time

            logger.debug("HBond analysis is recording time step")
        except NotImplementedError:
            # chained reader or xyz(?) cannot do time yet
            def _get_timestep():
                return self.u.trajectory.frame

            logger.warn("HBond analysis is recording frame number instead of time step")

        logger.info("Starting analysis (frame index start=%d stop=%d, step=%d)",
                    (self.traj_slice.start or 0),
                    (self.traj_slice.stop or self.u.trajectory.n_frames), self.traj_slice.step or 1)

        for progress, ts in enumerate(self.u.trajectory[self.traj_slice]):
            # all bonds for this timestep
            frame_results = []
            # dict of tuples (atomid, atomid) for quick check if
            # we already have the bond (to avoid duplicates)
            already_found = {}

            frame = ts.frame
            timestep = _get_timestep()
            self.timesteps.append(timestep)

            pm.echo(progress, current_step=frame)
            self.logger_debug("Analyzing frame %(frame)d, timestep %(timestep)f ps", vars())
            if self.update_selection1:
                self._update_selection_1()
            if self.update_selection2:
                self._update_selection_2()

            if self.selection1_type in ('donor', 'both') and self._s2_acceptors:
                self.logger_debug("Selection 1 Donors <-> Acceptors")
                ns_acceptors = AtomNeighborSearch(self._s2_acceptors)
                for i, donor_h_set in self._s1_donors_h.items():
                    d = self._s1_donors[i]
                    for h in donor_h_set:
                        res = ns_acceptors.search(h, self.distance)
                        for a in res:
                            angle = self.calc_angle(d, h, a)
                            donor_atom = h if self.distance_type != 'heavy' else d
                            dist = self.calc_eucl_distance(donor_atom, a)
                            if angle >= self.angle and dist <= self.distance:
                                self.logger_debug(
                                    "S1-D: {0!s} <-> S2-A: {1!s} {2:f} A, {3:f} DEG".format(h.index + 1, a.index + 1, dist, angle))
                                #self.logger_debug("S1-D: %r <-> S2-A: %r %f A, %f DEG" % (h, a, dist, angle))
                                frame_results.append(
                                    [h.index + 1, a.index + 1, h.index, a.index,
                                    '{0!s}{1!s}:{2!s}'.format(h.resname, repr(h.resid), h.name),
                                    '{0!s}{1!s}:{2!s}'.format(a.resname, repr(a.resid), a.name),
                                    dist, angle])

                                already_found[(h.index + 1, a.index + 1)] = True
            if self.selection1_type in ('acceptor', 'both') and self._s1_acceptors:
                self.logger_debug("Selection 1 Acceptors <-> Donors")
                ns_acceptors = AtomNeighborSearch(self._s1_acceptors)
                for i, donor_h_set in self._s2_donors_h.items():
                    d = self._s2_donors[i]
                    for h in donor_h_set:
                        res = ns_acceptors.search(h, self.distance)
                        for a in res:
                            if remove_duplicates and (
                                    (h.index + 1, a.index + 1) in already_found
                                    or (a.index + 1, h.index + 1) in already_found):
                                continue
                            angle = self.calc_angle(d, h, a)
                            donor_atom = h if self.distance_type != 'heavy' else d
                            dist = self.calc_eucl_distance(donor_atom, a)
                            if angle >= self.angle and dist <= self.distance:
                                self.logger_debug(
                                    "S1-A: {0!s} <-> S2-D: {1!s} {2:f} A, {3:f} DEG".format(a.index + 1, h.index + 1, dist, angle))
                                #self.logger_debug("S1-A: %r <-> S2-D: %r %f A, %f DEG" % (a, h, dist, angle))
                                frame_results.append(
                                    [h.index + 1, a.index + 1, h.index, a.index,
                                    '{0!s}{1!s}:{2!s}'.format(h.resname, repr(h.resid), h.name),
                                    '{0!s}{1!s}:{2!s}'.format(a.resname, repr(a.resid), a.name),
                                    dist, angle])

            self.timeseries.append(frame_results)

        logger.info("HBond analysis: complete; timeseries with %d hbonds in %s.timeseries",
                    self.count_by_time().count.sum(), self.__class__.__name__)

    @staticmethod
    def calc_angle(d, h, a):
        """Calculate the angle (in degrees) between two atoms with H at apex."""
        v1 = h.position - d.position
        v2 = h.position - a.position
        if np.all(v1 == v2):
            return 0.0
        return np.rad2deg(angle(v1, v2))

    @staticmethod
    def calc_eucl_distance(a1, a2):
        """Calculate the Euclidean distance between two atoms. """
        return norm(a2.position - a1.position)

    def generate_table(self):
        """Generate a normalised table of the results.

        The table is stored as a :class:`numpy.recarray` in the
        attribute :attr:`~HydrogenBondAnalysis.table` and can be used
        with e.g. `recsql`_.

        Columns:
          0. "time"
          1. "donor_idx"
          2. "acceptor_idx"
          3. "donor_index"
          4. "acceptor_index"
          4. "donor_resnm"
          5. "donor_resid"
          6. "donor_atom"
          7. "acceptor_resnm"
          8. "acceptor_resid"
          9. "acceptor_atom"
          10. "distance"
          11. "angle"

        .. _recsql: http://pypi.python.org/pypi/RecSQL
        """
        if self.timeseries is None:
            msg = "No timeseries computed, do run() first."
            warnings.warn(msg, category=MissingDataWarning)
            logger.warn(msg)
            return

        num_records = np.sum([len(hframe) for hframe in self.timeseries])
        # build empty output table
        dtype = [
            ("time", float), ("donor_idx", int), ("acceptor_idx", int),
            ("donor_index", int),  ("acceptor_index", int),
            ("donor_resnm", "|U4"), ("donor_resid", int), ("donor_atom", "|U4"),
            ("acceptor_resnm", "|U4"), ("acceptor_resid", int), ("acceptor_atom", "|U4"),
            ("distance", float), ("angle", float)]
        # according to Lukas' notes below, using a recarray at this stage is ineffective
        # and speedups of ~x10 can be achieved by filling a standard array, like this:
        out = np.empty((num_records,), dtype=dtype)
        cursor = 0  # current row
        for t, hframe in zip(self.timesteps, self.timeseries):
            for (donor_idx, acceptor_idx, donor_index, acceptor_index, donor,
            acceptor, distance, angle) in hframe:
                out[cursor] = (t, donor_idx, acceptor_idx, donor_index, acceptor_index) + \
                parse_residue(donor) + parse_residue(acceptor) + (distance, angle)
                cursor += 1
        assert cursor == num_records, "Internal Error: Not all HB records stored"
        self.table = out.view(np.recarray)
        logger.debug("HBond: Stored results as table with %(num_records)d entries.", vars())

    def save_table(self, filename="hbond_table.pickle"):
        """Saves :attr:`~HydrogenBondAnalysis.table` to a pickled file.

        Load with ::

           import cPickle
           table = cPickle.load(open(filename))

        .. SeeAlso:: :mod:`cPickle` module and :class:`numpy.recarray`
        """
        if self.table is None:
            self.generate_table()
        cPickle.dump(self.table, open(filename, 'wb'), protocol=cPickle.HIGHEST_PROTOCOL)

    def count_by_time(self):
        """Counts the number of hydrogen bonds per timestep.

        :Returns: a class:`numpy.recarray`
        """

        if self.timeseries is None:
            msg = "No timeseries computed, do run() first."
            warnings.warn(msg, category=MissingDataWarning)
            logger.warn(msg)
            return

        out = np.empty((len(self.timesteps),), dtype=[('time', float), ('count', int)])
        for cursor, time_count in enumerate(zip(self.timesteps,
                                                map(len, self.timeseries))):
            out[cursor] = time_count
        return out.view(np.recarray)

    def count_by_type(self):
        """Counts the frequency of hydrogen bonds of a specific type.

        Processes :attr:`HydrogenBondAnalysis.timeseries` and returns
        a :class:`numpy.recarray` containing atom indices, residue
        names, residue numbers (for donors and acceptors) and the
        fraction of the total time during which the hydrogen bond was
        detected.

        :Returns: a class:`numpy.recarray`
        """
        if self.timeseries is None:
            msg = "No timeseries computed, do run() first."
            warnings.warn(msg, category=MissingDataWarning)
            logger.warn(msg)
            return

        hbonds = defaultdict(int)
        for hframe in self.timeseries:
            for (donor_idx, acceptor_idx, donor_index, acceptor_index, donor,
                acceptor, distance, angle) in hframe:
                donor_resnm, donor_resid, donor_atom = parse_residue(donor)
                acceptor_resnm, acceptor_resid, acceptor_atom = parse_residue(acceptor)
                # generate unambigous key for current hbond \
                # (the donor_heavy_atom placeholder '?' is added later)
                # idx_zero is redundant for an unambigous key, but included for
                # consistency.
                hb_key = (
                    donor_idx, acceptor_idx, donor_index, acceptor_index,
                    donor_resnm, donor_resid, "?", donor_atom,
                    acceptor_resnm, acceptor_resid, acceptor_atom)

                hbonds[hb_key] += 1

        # build empty output table
        dtype = [
            ('donor_idx', int), ('acceptor_idx', int),
            ("donor_index", int), ("acceptor_index", int), ('donor_resnm', 'U4'),
            ('donor_resid', int), ('donor_heavy_atom', 'U4'), ('donor_atom', 'U4'),
            ('acceptor_resnm', 'U4'), ('acceptor_resid', int), ('acceptor_atom', 'U4'),
            ('frequency', float)
        ]
        out = np.empty((len(hbonds),), dtype=dtype)

        # float because of division later
        tsteps = float(len(self.timesteps))
        for cursor, (key, count) in enumerate(six.iteritems(hbonds)):
            out[cursor] = key + (count / tsteps,)

        # return array as recarray
        # The recarray has not been used within the function, because accessing the
        # the elements of a recarray (3.65 us) is much slower then accessing those
        # of a ndarray (287 ns).
        r = out.view(np.recarray)

        # patch in donor heavy atom names (replaces '?' in the key)
        h2donor = self._donor_lookup_table_byindex()
        r.donor_heavy_atom[:] = [h2donor[idx] for idx in r.donor_index]

        return r

    def timesteps_by_type(self):
        """Frames during which each hydrogen bond existed, sorted by hydrogen bond.

        Processes :attr:`HydrogenBondAnalysis.timeseries` and returns
        a :class:`numpy.recarray` containing atom indices, residue
        names, residue numbers (for donors and acceptors) and a list
        of timesteps at which the hydrogen bond was detected.

        :Returns: a class:`numpy.recarray`
        """

        if self.timeseries is None:
            msg = "No timeseries computed, do run() first."
            warnings.warn(msg, category=MissingDataWarning)
            logger.warn(msg)
            return

        hbonds = defaultdict(list)
        for (t, hframe) in zip(self.timesteps, self.timeseries):
            for (donor_idx, acceptor_idx, donor_index, acceptor_index, donor,
            acceptor, distance, angle) in hframe:
                donor_resnm, donor_resid, donor_atom = parse_residue(donor)
                acceptor_resnm, acceptor_resid, acceptor_atom = parse_residue(acceptor)
                # generate unambigous key for current hbond
                # (the donor_heavy_atom placeholder '?' is added later)
                # idx_zero is redundant for key but added for consistency
                hb_key = (
                    donor_idx, acceptor_idx,
                    donor_index, acceptor_index,
                    donor_resnm, donor_resid, "?", donor_atom,
                    acceptor_resnm, acceptor_resid, acceptor_atom)
                hbonds[hb_key].append(t)

        out_nrows = 0
        # count number of timesteps per key to get length of output table
        for ts_list in six.itervalues(hbonds):
            out_nrows += len(ts_list)

        # build empty output table
        dtype = [
            ('donor_idx', int), ('acceptor_idx', int),('donor_index', int),
            ('acceptor_index', int), ('donor_resnm', 'U4'), ('donor_resid', int),
            ('donor_heavy_atom', 'U4'), ('donor_atom', 'U4'),('acceptor_resnm', 'U4'),
            ('acceptor_resid', int), ('acceptor_atom', 'U4'), ('time', float)]
        out = np.empty((out_nrows,), dtype=dtype)

        out_row = 0
        for (key, times) in six.iteritems(hbonds):
            for tstep in times:
                out[out_row] = key + (tstep,)
                out_row += 1

        # return array as recarray
        # The recarray has not been used within the function, because accessing the
        # the elements of a recarray (3.65 us) is much slower then accessing those
        # of a ndarray (287 ns).
        r = out.view(np.recarray)

        # patch in donor heavy atom names (replaces '?' in the key)
        h2donor = self._donor_lookup_table_byindex()
        r.donor_heavy_atom[:] = [h2donor[idx - 1] for idx in r.donor_idx]

        return r

    def _donor_lookup_table_byres(self):
        """Look-up table to identify the donor heavy atom from resid and hydrogen name.

        Assumptions:
        * resids are unique
        * hydrogen atom names are unique within a residue
        * selections have not changed (because we are simply looking at the last content
          of the donors and donor hydrogen lists)

        Donors from *selection1* and *selection2* are merged.

        Output dictionary ``h2donor`` can be used as::

           heavy_atom_name = h2donor[resid][hydrogen_name]

        """
        s1d = self._s1_donors  # list of donor Atom instances
        s1h = self._s1_donors_h  # dict indexed by donor position in donor list, containg AtomGroups of H
        s2d = self._s2_donors
        s2h = self._s2_donors_h

        def _make_dict(donors, hydrogens):
            # two steps so that entry for one residue can be UPDATED for multiple donors
            d = dict((donors[k].resid, {}) for k in range(len(donors)) if k in hydrogens)
            for k in range(len(donors)):
                if k in hydrogens:
                    d[donors[k].resid].update(dict((atom.name, donors[k].name) for atom in hydrogens[k]))
            return d

        h2donor = _make_dict(s2d, s2h)  # 2 is typically the larger group
        # merge (in principle h2donor.update(_make_dict(s1d, s1h) should be sufficient
        # with our assumptions but the following should be really safe)
        for resid, names in _make_dict(s1d, s1h).items():
            if resid in h2donor:
                h2donor[resid].update(names)
            else:
                h2donor[resid] = names

        return h2donor

    def _donor_lookup_table_byindex(self):
        """Look-up table to identify the donor heavy atom from hydrogen atom index.

        Assumptions:
        * selections have not changed (because we are simply looking at the last content
          of the donors and donor hydrogen lists)

        Donors from *selection1* and *selection2* are merged.

        Output dictionary ``h2donor`` can be used as::

           heavy_atom_name = h2donor[index]

        .. Note::

           *index* is the 0-based MDAnalysis index
           (:attr:`MDAnalysis.core.groups.Atom.index`).  The
           tables generated by :class:`HydrogenBondAnalysis` contain
           1-based indices and zero-based indices.

        .. deprecated:: 0.15.0
            The 1-based indices are deprecated in favor of the zero-based indices
            given by "idx_zero".
        """
        s1d = self._s1_donors  # list of donor Atom instances
        s1h = self._s1_donors_h  # dict indexed by donor position in donor list, containg AtomGroups of H
        s2d = self._s2_donors
        s2h = self._s2_donors_h

        def _make_dict(donors, hydrogens):
            #return dict(flatten_1([(atom.id, donors[k].name) for atom in hydrogens[k]] for k in range(len(donors))
            # if k in hydrogens))
            x = []
            for k in range(len(donors)):
                if k in hydrogens:
                    x.extend([(atom.index, donors[k].name) for atom in hydrogens[k]])
            return dict(x)

        h2donor = _make_dict(s2d, s2h)  # 2 is typically the larger group
        h2donor.update(_make_dict(s1d, s1h))

        return h2donor
