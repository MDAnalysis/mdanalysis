from base import AnalysisBase
import numpy as np

class TotalDipole(AnalysisBase):
    
    def __init__(self, universe=None, filename='order.dat', selection=None, start=None, stop=None, step=None):

        if selection is None:
            raise RuntimeError('In class TotalDipole: constroctur requires a selection')
        else:
            self.selection_string        = selection

        self._universe = universe
        #self._trajectory = self._universe.trajectory
        self.filename         = filename
        self.time             = []

        self._setup_frames(self._universe, start, stop, step)
        self.dipoles          = []

    def _single_frame(self):
        selection = self._universe.select_atoms(self.selection_string)

        dipole = np.zeros(3)

        for residue in selection.residues:
            dipole += MolecularDipole(residue)

        self.dipoles.append(dipole)

    def _conclude(self):
        total_dipole = np.sum(self.dipoles, axis=0)
        print "Average dipole:", total_dipole/self.nframes
        return total_dipole/self.nframes

    def __iadd__(self,other):
        self.dipoles += other.dipoles
        self.nframes += other.nframes
        return self


#--- Functions ---#
def MolecularDipole(residue):
    charges    = residue.charges
    abscharges = np.absolute(charges)
    charge_sum = np.sum(abscharges)
    positions  = residue.positions

    charge_center = []

    for coord in [0,1,2]:
        charge_center.append(np.sum(np.multiply(abscharges,positions[:,coord]))/charge_sum)

    dipole = []
    # 4.803 converts to debyes
    for coord in [0,1,2]:
        dipole.append(np.sum(np.multiply(charges,positions[:,coord]-charge_center[coord]))*4.803)
                      
    return dipole 
