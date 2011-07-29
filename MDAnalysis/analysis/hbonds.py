# Hydrogen Bonding Analysis
"""
Hydrogen Bond analysis --- :mod:`MDAnalysis.analysis.hbonds`
===================================================================

:Author: David Caplan
:Year: 2010
:Copyright: GNU Public License v3


This is modeled after the VMD HBONDS plugin (http://www.ks.uiuc.edu/Research/vmd/plugins/hbonds/)

Given a Universe (simulation trajectory with 1 or more frames) measure all hydrogen bonds for each frame between selections 1 and 2.

Options:
  - update_selections (True): update selections at each frame?
  - selection_1_type ('both'): selection 1 is the: donor, acceptor, both
  - donor-acceptor distance (A): 3.0
  - Angle cutoff (degrees): 120.0
  - donor and acceptor atom types

Returns hydrogen bond data per frame:
    results[ [ <donor index>, <acceptor index>, <donor string>, <acceptor string>, <distance>, <angle> ], [frame 1], [frame 2] ... ]


Example
-------

TODO

  import MDAnalysis
  import hbonds
  
  u = MDAnalysis.Universe(PSF, PDB, permissive=True)
  h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u, 'protein', 'resname TIP3', distance=3.0, angle=120.0)
  results = h.run()

Classes
-------

.. autoclass:: HydrogenBondAnalysis
   :members:

"""

from __future__ import with_statement

import numpy

from MDAnalysis.core.AtomGroup import AtomGroup
import MDAnalysis.KDTree.NeighborSearch as NS

import logging
logger = logging.getLogger('MDAnalysis.analysis.hbonds')

DEFAULT_DONORS = ('NH', 'OH2', 'OW', 'NE', 'ND2', 'NE2', 'OG', 'OH', 'NH1', 'SG', 'ND1', 'OG1', 'NH2', 'NE2', 'NZ', 'NE1', )
DEFAULT_ACCEPTORS = ('CO', 'OH2', 'OW', 'OD1', 'OE1', 'SD', 'OD1', 'OE1', 'OG', 'OD2', 'OE2', 'OG1', 'SG', 'ND1', 'OH', )

class HydrogenBondAnalysis(object):
    """Perform a hydrogen bond analysis

    The analysis of the trajectory is performed with the
    :meth:`HydrogenBondAnalysis.run` method. The result is stored in
    :attr:`HydrogenBondAnalysis.timeseries`. It is a numpy array which
    contains the frame number at index 0, selection1 and selection2,
    and the total number of hydrogen bonds ::

        frame  selection1 selection2  num_hbonds


    Donors: NH of the main chain, water-H1/H2, ARG NE, ASN ND2, HIS NE2, SER OG, TYR OH, ARG NH1, CYS SG, HIS ND1, THR OG1, ARG NH2, GLN NE2, LYS NZ, TRP NE1
    Acceptors: CO main chain, water-OH2, water-OW, ASN OD1, GLN OE1, MET SD, ASP OD1, GLU OE1, SER OG, ASP OD2, GLU OE2, THR OG1,CYH SG, HIS ND1, TYR OH.

    """
    
    def __init__(self, universe, selection1='protein', selection2='all', selection1_type='both',
                update_selection1=False, update_selection2=False, filter_first=True, distance=3.0, angle=120.0,
                donors=[], acceptors=[]):
        """Calculate hydrogen bonds between two selections in a universe.

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
            Extra H donor atom types (in addition to those in DEFAULT_DONORS)
          *acceptors*
            Extra H acceptor atom types (in addition to those in DEFAULT_ACCEPTORS)
            
        The timeseries accessible as the attribute :attr:`HydrogenBondAnalysis.timeseries`.
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
        self.donors = DEFAULT_DONORS + tuple(donors)
        self.acceptors = DEFAULT_ACCEPTORS + tuple(acceptors)
        
        if not (self.selection1 and self.selection2):
            raise Exception('HydrogenBondAnalysis: invalid selections')
        elif self.selection1_type not in ('both', 'donor', 'acceptor'):
            raise Exception('HydrogenBondAnalysis: Invalid selection type %s' % self.selection1_type)

        self.timeseries = None  # final result
        
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
        """
        self.timeseries = []
        
        for ts in self.u.trajectory:
            frame_results = []
            
            frame = ts.frame
            logger.debug("Analyzing frame %d" % frame)
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
        return self.timeseries
        
    def calc_angle(self, d, h, a):
        """Calculate the angle (in degrees) between two atoms
        """
        v1 = h.pos-d.pos
        v2 = h.pos-a.pos
        v1 /= numpy.linalg.norm(v1)
        v2 /= numpy.linalg.norm(v2)

        if v1.tolist() == v2.tolist():
            return 0.0
        return numpy.arccos(numpy.dot(v1, v2) / (numpy.linalg.norm(v1)*numpy.linalg.norm(v2))) * 180 / numpy.pi


    def calc_eucl_distance(self, a1, a2):
        """Calculate the Euclidean distance between two atoms.
        """
        v1 = a1.pos
        v2 = a2.pos
        return numpy.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)

