# $Id$
"""Hydrogen bond analysis in MDAnalysis

Hydrogen bond ... :      x-A...H-D
   A: acceptor
   x: acceptor antecendent
   H: hydrogen
   D: donor (hydrogen antecendent)

Requires a PSF file with DONOR and ACCEPTOR entries. PSF entries are
pair-tuples of A/D and x/H:

   acceptors:     (A,x)      [if x = -1 then no antecedent]
   donors:        (D,H)

(Note that conceptually the order of 'atom of interest' (A,H) and
antecedent (x,D) is switched between acceptors and donors.)
  
"""
import numpy
import AtomGroup
#from sets import Set as set   # works for python 2.3, py2.4 has set builtin

# organize topology as AtomGroup objects

def init_hbonds(u):
    """Topology data for hydrogen bonds: Store donor and acceptor data with
    individual atoms.

    ACCEPTOR:
    atom.acceptor = Atom(antecedent)      # reference to antecedent atom x (not 
                                          # really used at the moment, just make 
                                          # sure its not None
    DONOR:
    atom.donor = AtomGroup([H1, H2, ...]) # AtomGroup of all HB hydrogens from this atom

    HYDROGEN:
    atom.donor = Atom(donor)              # donor of the hydrogen

    all other atoms:
    atom.acceptor = atom.donor = None
    """
    # u = universe
    atoms = u.atoms
        
    u._HBonds = {}
    u._HBonds['acceptors'] = AtomGroup([atoms[iacc] for (iacc,iantecedent) in u._psf['_acceptors']])
    donors = numpy.unique([idonor for (idonor,ihydrogen) in u._psf['_donors']])
    u._HBonds['donors'] = AtomGroup([atoms[i] for i in donors])
    u._HBonds['hydrogens'] = AtomGroup([atoms[ihydrogen] for (idonor,ihydrogen) in u._psf['_donors']])
    
    # set donor/acceptor only for specific ones; note that antecedent = -1 is valid
    # To test for existence one should check 'if atom.donor is not None'
    for (iacc,iantecedent) in u._psf['_acceptors']:
        if iantecedent == -1:
            iantecedent = iacc    # point to itself if no explicit antecedent given
        atoms[iacc].acceptor = atoms[iantecedent]
        
    for (idonor,ihydrogen) in u._psf['_donors']:
        atoms[ihydrogen].donor = atoms[idonor]
        try:
            atoms[idonor].donor += atoms[ihydrogen]
        except TypeError:
            atoms[idonor].donor = AtomGroup([atoms[ihydrogen]])

class HBonds:
    def __init__(self,universe):
        init_hbonds(universe)
        self.universe = universe
        HBonds = dict( [(x,u._HBonds[x]) for x in 'acceptors','donors','hydrogens'] )
        self.__dict__.update(HBonds)
        # coordinate functions
        self.h_coords = self.hydrogens.coordinates
        self.d_coords = self.donors.coordinates
        self.a_coords = self.acceptors.coordinates

def dist(i,j):
    global d
    print "dist(%d,%d) = %f" % (i,j,d[i,j])
    print "acceptor: "+str(HB.acceptors[i])
    print "hydrogen: "+str(HB.hydrogens[j])
    print "donor:    "+str(HB.hydrogens[j].donor)

    
if __name__ == '__main__':
    import MDAnalysis
    u = MDAnalysis.Universe(psffilename="inp/ifabp_water_hbond.psf",dcdfilename="trj/rmsfit_ifabp_water_1.dcd")

    HB = HBonds(u)
    print "HB -- contains hydrogen bonds of "+str(u)+" for easy access"
    import MDAnalysis.distances

    #d = MDAnalysis.distances.distance_array(A_coords(),h_coords(),box=None)

