# -*- encoding: utf-8 -*-
# Hydrogen Bonding Analysis
"""
Hydrogen Bond analysis --- :mod:`MDAnalysis.analysis.hbonds`
===================================================================

:Author: David Caplan
:Year: 2010-2011
:Copyright: GNU Public License v3


Given a :class:`~MDAnalysis.core.AtomGroup.Universe` (simulation
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
  - *donors* and *acceptors* atom types (to add additional atom names)


.. _Analysis Output:

Output
------

The results are hydrogen bond data per frame (# indicates comments that are not part of the output)::

    results = [
        [ # frame 1
           [ # hbond 1
              <donor index>, <acceptor index>, <donor string>, <acceptor string>, <distance>, <angle>
           ],
           [ # hbond 2
              <donor index>, <acceptor index>, <donor string>, <acceptor string>, <distance>, <angle>
           ],
           ....
        ],
        [ # frame 2
          [ ... ], [ ... ], ...
        ],
        ...
    ]

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
defaults are appropriate for the CHARMM27 force field as defined in Table
`Default atom names for hydrogen bonding analysis`_. Hydrogen atoms are
searched for as atom names following in sequence and starting with 'H'.


.. _Default atom names for hydrogen bonding analysis:

.. table:: Default heavy atom names for hydrogen bonding analysis.

   =========== ==============  =========== ====================================
   group       donor           acceptor    comments
   =========== ==============  =========== ====================================
   main chain  N               C
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


Donor and acceptor names are based on the CHARMM27 force field but will also
work for e.g. OPLS/AA (tested in Gromacs) . Residue names in the table are for
information only and are not taken into account when determining acceptors and
donors. This can potentially lead to some ambiguity in the assignment of
donors/acceptors for residues such as histidine or cytosine.

The lists of donor and acceptor names can be extended by providing lists of
atom names in the *donors* and *acceptors* keywords to
:class:`HydrogenBondAnalysis`. If the lists are entirely inapprpriate
(e.g. when analysing simulations done with a force field that uses very
different atom names or when analysing non-protein system) then one can derive
a new class and set the default lists oneself::

 class HydrogenBondAnalysis_OtherFF(HydrogenBondAnalysis):
       DEFAULT_DONORS = (....)      # add donor heavy atoms here
       DEFAULT_ACCEPTORS = (....)   # add acceptor heavy atoms here

Then simply use the new class instead of the parent class.


.. rubric:: References

.. [Gregoret1991] L.M. Gregoret, S.D. Rader, R.J. Fletterick, and
   F.E. Cohen. Hydrogen bonds involving sulfur atoms in proteins. Proteins,
   9(2):99–107, 1991. `10.1002/prot.340090204`_.

.. _`10.1002/prot.340090204`: http://dx.doi.org/10.1002/prot.340090204


Example
-------

All protein-water hydrogen bonds can be analysed with ::

  import MDAnalysis.analysis.hbonds

  u = MDAnalysis.Universe(PSF, PDB, permissive=True)
  h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u, 'protein', 'resname TIP3', distance=3.0, angle=120.0)
  results = h.run()

The results are also stored as the attribute :attr:`h.timeseries`; see
:ref:`Analysis Output` for the format and further options.

.. Note:: Due to the way :class:`HydrogenBondAnalysis` is implemented, it is
          more efficient to have the second selection (*selection2*) be the
          *larger* group, e.g. the water when looking at water-protein H-bonds
          or the whole protein when looking at ligand-protein interactions.

.. TODO: how to analyse the ouput and notes on selection updating


Classes
-------

.. autoclass:: HydrogenBondAnalysis
   :members:

"""

from __future__ import with_statement

import numpy

from MDAnalysis import MissingDataWarning
from MDAnalysis.core.AtomGroup import AtomGroup
import MDAnalysis.KDTree.NeighborSearch as NS
from MDAnalysis.core.util import norm, angle, parse_residue

import warnings
import logging
logger = logging.getLogger('MDAnalysis.analysis.hbonds')

class HydrogenBondAnalysis(object):
    """Perform a hydrogen bond analysis

    The analysis of the trajectory is performed with the
    :meth:`HydrogenBondAnalysis.run` method. The result is stored in
    :attr:`HydrogenBondAnalysis.timeseries`. See
    :meth:`~HydrogenBondAnalysis.run` for the format.

    The default atom names are taken from the CHARMM 27 force field files but
    also work for e.g. OPLS/AA in Gromacs.

    *Donors* (associated hydrogens are deduced from topology)
      N of the main chain, water OH2/OW, ARG NE/NH1/NH2, ASN ND2, HIS ND1/NE2,
      SER OG, TYR OH, CYS SG, THR OG1, GLN NE2, LYS NZ, TRP NE1

    *Acceptors*
      O of the main chain, water OH2/OW, ASN OD1, ASP OD1/OD2, CYH SG, GLN OE1,
      GLU OE1/OE2, HIS ND1/NE2, MET SD, SER OG, THR OG1, TYR OH

    .. SeeAlso:: Table :ref:`Default atom names for hydrogen bonding analysis`

    """

    # use tuple(set()) here so that one can just copy&paste names from the
    # table; set() takes care for removing duplicates. At the end the
    # DEFAULT_DONORS and DEFAULT_ACCEPTORS should simply be tuples.

    #: default heavy atom names whose hydrogens are treated as *donors*
    #: (see :ref:`Default atom names for hydrogen bonding analysis`)
    #: Use the keyword *donors* to add a list of additional donor names.
    DEFAULT_DONORS = tuple(set(['N', 'OH2', 'OW', 'NE', 'NH1' 'NH2', 'ND2', 'SG', 'NE2',
                               'ND1', 'NE2', 'ND1', 'NE2', 'ND1', 'NE2', 'NZ', 'OG', 'OG1',
                               'NE1', 'OH',]))

    #: default atom names that are treated as hydrogen *acceptors*
    #: (see :ref:`Default atom names for hydrogen bonding analysis`)
    #: Use the keyword *acceptors* to add a list of additional acceptor names.
    DEFAULT_ACCEPTORS = tuple(set(['C', 'OH2', 'OW', 'OD1', 'OD2', 'SG', 'OE1', 'OE1', 'OE2',
                                   'ND1', 'NE2', 'NE2', 'ND1', 'SD', 'OG', 'OG1', 'OH',]))

    def __init__(self, universe, selection1='protein', selection2='all', selection1_type='both',
                update_selection1=False, update_selection2=False, filter_first=True, distance=3.0, angle=120.0,
                donors=None, acceptors=None):
        """Set up calculation of hydrogen bonds between two selections in a universe.

        :Arguments:
          *universe*
            Universe object
          *selection1*
            Selection string for first selection
          *selection2*
            Selection string for second selection
          *selection1_type*
            Selection 1 can be 'donor', 'acceptor' or 'both'
          *update_selection1*
            Update selection 1 at each frame?
          *update_selection2*
            Update selection 2 at each frame?
          *filter_first*
            Filter selection 2 first to only atoms 3*distance away
          *distance*
            Distance cutoff for hydrogen bonds
          *angle*
            Angle cutoff for hydrogen bonds
          *donors*
            Extra H donor atom types (in addition to those in
            :attr:`~HydrogenBondAnalysis.DEFAULT_DONORS`), must be a sequence.
          *acceptors*
            Extra H acceptor atom types (in addition to those in
            :attr:`~HydrogenBondAnalysis.DEFAULT_ACCEPTORS`), must be a sequence.

        The timeseries is accessible as the attribute :attr:`HydrogenBondAnalysis.timeseries`.
        """

        self.u = universe
        self.selection1 = selection1
        self.selection2 = selection2
        self.selection1_type = selection1_type
        self.update_selection1 = update_selection1
        self.update_selection2 = update_selection2
        self.filter_first = filter_first
        self.distance = distance
        self.angle = angle

        # set up the donors/acceptors lists
        if donors is None:
            donors = []
        if acceptors is None:
            acceptors = []
        self.donors = self.DEFAULT_DONORS + tuple(donors)
        self.acceptors = self.DEFAULT_ACCEPTORS + tuple(acceptors)

        if not (self.selection1 and self.selection2):
            raise ValueError('HydrogenBondAnalysis: invalid selections')
        elif self.selection1_type not in ('both', 'donor', 'acceptor'):
            raise ValueError('HydrogenBondAnalysis: Invalid selection type %s' % self.selection1_type)

        self.timeseries = None  # final result
        self._timesteps = None  # time for each frame

        self._update_selection_1()
        self._update_selection_2()

    def _get_bonded_hydrogens(self, atom):
        hydrogens = []
        for i in range(3):
            try:
                next_atom = self.u.atoms[atom.number+1+i]
            except:
                break
            else:
                if next_atom.name.startswith('H'):
                    hydrogens.append(next_atom)
        return hydrogens

    def _update_selection_1(self):
        self._s1 = self.u.selectAtoms(self.selection1)
        logger.debug("Size of selection 1: %d atoms" % len(self._s1))
        if self.selection1_type in ('donor', 'both'):
            self._s1_donors = self._s1.selectAtoms(' or '.join([ 'name %s' % i for i in self.donors ]))
            self._s1_donors_h = {}
            for i,d in enumerate(self._s1_donors):
                tmp = self._get_bonded_hydrogens(d)
                if tmp:
                    self._s1_donors_h[i] = tmp
            logger.debug("Selection 1 donors: %d" % len(self._s1_donors))
            logger.debug("Selection 1 donor hydrogens: %d" % len(self._s1_donors_h.keys()))
        if self.selection1_type in ('acceptor', 'both'):
            self._s1_acceptors = self._s1.selectAtoms(' or '.join([ 'name %s' % i for i in self.acceptors ]))
            logger.debug("Selection 1 acceptors: %d" % len(self._s1_acceptors))

    def _update_selection_2(self):
        self._s2 = self.u.selectAtoms(self.selection2)
        if self.filter_first:
            logger.debug("Size of selection 2 before filtering: %d atoms" % len(self._s2))
            ns_selection_2 = NS.AtomNeighborSearch(self._s2)
            self._s2 = ns_selection_2.search_list(self._s1, 3.*self.distance)
        logger.debug("Size of selection 2: %d atoms" % len(self._s2))

        if self.selection1_type in ('donor', 'both'):
            self._s2_acceptors = self._s2.selectAtoms(' or '.join([ 'name %s' % i for i in self.acceptors ]))
            logger.debug("Selection 2 acceptors: %d" % len(self._s2_acceptors))
        if self.selection1_type in ('acceptor', 'both'):
            self._s2_donors = self._s2.selectAtoms(' or '.join([ 'name %s' % i for i in self.donors ]))
            self._s2_donors_h = {}
            for i,d in enumerate(self._s2_donors):
                tmp = self._get_bonded_hydrogens(d)
                if tmp:
                    self._s2_donors_h[i] = tmp
            logger.debug("Selection 2 donors: %d" % len(self._s2_donors))
            logger.debug("Selection 2 donor hydrogens: %d" % len(self._s2_donors_h.keys()))

    def run(self):
        """Analyze trajectory and produce timeseries.

        Returns hydrogen bond data per frame (# indicates comments that are not part of the output)::

          results = [
              [ # frame 1
                 [ # hbond 1
                    <donor index>, <acceptor index>, <donor string>, <acceptor string>, <distance>, <angle>
                 ],
                 [ # hbond 2
                    <donor index>, <acceptor index>, <donor string>, <acceptor string>, <distance>, <angle>
                 ],
                 ....
              ],
              [ # frame 2
                [ ... ], [ ... ], ...
              ],
              ...
          ]

        The data are also stored as :attr:`HydrogenBondAnalysis.timeseries`.

        """
        logger.info("HBond analysis: starting")

        self.timeseries = []
        self._timesteps = []

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

        for ts in self.u.trajectory:
            frame_results = []

            frame = ts.frame
            timestep = _get_timestep()
            self._timesteps.append(timestep)

            logger.debug("Analyzing frame %(frame)d, timestep %(timestep)f ps", vars())
            if self.update_selection1:
                self._update_selection_1()
            if self.update_selection2:
                self._update_selection_2()

            if self.selection1_type in ('donor', 'both'):
                logger.debug("Selection 1 Donors <-> Acceptors")
                ns_acceptors = NS.AtomNeighborSearch(self._s2_acceptors)
                for i,donor_h_set in self._s1_donors_h.items():
                    d = self._s1_donors[i]
                    for h in donor_h_set:
                        res = ns_acceptors.search_list(AtomGroup([h]), self.distance)
                        for a in res:
                            angle = self.calc_angle(d,h,a)
                            dist = self.calc_eucl_distance(h,a)
                            if angle >= self.angle and dist <= self.distance:
                                logger.debug("S1-D: %s <-> S2-A: %s %fA, %f DEG" % (h.number, a.number, dist, angle))
                                frame_results.append([h.number+1, a.number+1, '%s%s:%s' % (h.resname, repr(h.resid), h.name), '%s%s:%s' % (a.resname, repr(a.resid), a.name), dist, angle])
            if self.selection1_type in ('acceptor', 'both'):
                logger.debug("Selection 1 Acceptors <-> Donors")
                ns_acceptors = NS.AtomNeighborSearch(self._s1_acceptors)
                for i,donor_h_set in self._s2_donors_h.items():
                    d = self._s2_donors[i]
                    for h in donor_h_set:
                        res = ns_acceptors.search_list(AtomGroup([h]), self.distance)
                        for a in res:
                            angle = self.calc_angle(d,h,a)
                            dist = self.calc_eucl_distance(h,a)
                            if angle >= self.angle and dist <= self.distance:
                                logger.debug("S1-A: %s <-> S2-D: %s %fA, %f DEG" % (a.number+1, h.number, dist, angle))
                                frame_results.append([h.number+1, a.number+1, '%s%s:%s' % (h.resname, repr(h.resid), h.name), '%s%s:%s' % (a.resname, repr(a.resid), a.name), dist, angle])
            self.timeseries.append(frame_results)

        logger.info("HBond analysis: complete; timeseries in %s.timeseries", self.__class__.__name__)
        return self.timeseries  # can we omitt this?

    def calc_angle(self, d, h, a):
        """Calculate the angle (in degrees) between two atoms with H at apex.
        """
        v1 = h.pos-d.pos
        v2 = h.pos-a.pos
        if numpy.all(v1 == v2):
            return 0.0
        return numpy.rad2deg(angle(v1, v2))

    def calc_eucl_distance(self, a1, a2):
        """Calculate the Euclidean distance between two atoms.
        """
        return norm(a2.pos - a1.pos)


    def generate_table(self):
        """Generate a normalised table of the results.

        The table is stored as a :class:`numpy.recarray` in the
        attribute :attr:`~HydrogenBondAnalysis.table` and can be used
        with e.g. `recsql`_.

        Columns:
          0. "time"
          1. "donor_idx"
          2. "acceptor_idx"
          3. "donor_resnm"
          4. "donor_resid"
          5. "donor_atom"
          6. "acceptor_resnm"
          7. "acceptor_resid"
          8. "acceptor_atom"
          9. "distance"
          10. "angle"

        .. _recsql: http://sbcb.bioch.ox.ac.uk/oliver/software/RecSQL/html/index.html
        """
        from itertools import izip
        if self.timeseries is None:
            msg = "No timeseries computed, do run() first."
            warnings.warn(msg)
            logger.warn(msg, category=MissingDataWarning)
            return

        num_records = numpy.sum([len(f) for f in self.timeseries])
        dtype = [("time",float), ("donor_idx",int), ("acceptor_idx",int),
                 ("donor_resnm","|S4"), ("donor_resid",int), ("donor_atom","|S4"),
                 ("acceptor_resnm","|S4"), ("acceptor_resid",int), ("acceptor_atom","|S4"),
                 ("distance",float), ("angle",float)]
        self.table = numpy.recarray((num_records,), dtype=dtype)
        row = 0
        for t,hframe in izip(self._timesteps, self.timeseries):
            self.table[row:row+len(hframe)].time = t
            for hrow, (donor_idx, acceptor_idx, donor, acceptor, distance, angle) in enumerate(hframe):
                cursor = row + hrow
                r = self.table[cursor]
                r.donor_idx = donor_idx
                r.donor_resnm, r.donor_resid, r.donor_atom = parse_residue(donor)
                r.acceptor_idx = acceptor_idx
                r.acceptor_resnm, r.acceptor_resid, r.acceptor_atom = parse_residue(acceptor)
                r.distance = distance
                r.angle = angle
            row = cursor + 1

    def save_table(self, filename="hbond_table.pickle"):
        """Saves the table to a pickled file.

        Load with ::

           import cPickle
           table = cPickle.load(open(filename))

        .. SeeAlso:: :mod:`cPickle` module and :class:`numpy.recarray`
        """
        import cPickle
        if self.table is None:
            self.generate_table()
        cPickle.dump(self.table, open(filename, 'wb'), protocol=cPickle.HIGHEST_PROTOCOL)
