from six.moves import zip
from six import string_types
import numpy as np
import logging
import warnings


import MDAnalysis.lib.qcprot as qcp
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.exceptions import SelectionError, NoDataError
from MDAnalysis.lib.log import ProgressMeter, _set_verbose
from MDAnalysis.lib.util import asiterable, iterable, get_weights


import numpy as np

def dihedral_calc(atomgroups):
    """Calculates phi and psi angles for a list of AtomGroups over trajectory.

    Parameters
    ----------
    atomgroups : list of AtomGroups
        must be a list of one or more AtomGroups containing 5 atoms in the
        correct order (i.e. C-N-CA-C-N)

    Returns
    -------
    angles : numpy.ndarray
        An array of time steps which contain (phi,psi) for all AtomGroups.
    """

    dihedrals = [atomgroup.dihedral for atomgroup in atomgroups]
    angles = [dih.value() for dih in dihedrals]

    return angles


class Ramachandran(AnalysisBase):
    """
    """
    def __init__(self, atomgroup, **kwargs):
        r"""
        """
        super(Ramachandran, self).__init__(atomgroup.universe.trajectory, **kwargs)
        self.atomgroup = atomgroup

    def run(self, start=None, stop=None, step=None, verbose=None, quiet=None):
        """Perform the analysis."""

        if any([el is not None for el in (start, stop, step, quiet)]):
            warnings.warn("run arguments are deprecated. Please pass them at "
                          "class construction. These options will be removed in 0.17.0",
                          category=DeprecationWarning)
            verbose = _set_verbose(verbose, quiet, default=False)
            # regenerate class with correct args
            super(Ramachandran, self).__init__(self.atomgroup.universe.trajectory,
                                       start=start, stop=stop, step=step,
                                       verbose=verbose)
        return super(Ramachandran, self).run()

    def _prepare(self):
        self.protein = self.atomgroup.universe.atoms.select_atoms("protein")
        self.resids = self.atomgroup.residues.resids
        self.phi_atoms = [self.protein.residues[resid-1].phi_selection()
                          for resid in self.resids
                          if 1 < resid < len(self.protein.residues)]
        self.psi_atoms = [self.protein.residues[resid-1].psi_selection()
                          for resid in self.resids
                          if 1 < resid < len(self.protein.residues)]
        self.angles = []

    def _single_frame(self):
        self.phi_angles = dihedral_calc(self.phi_atoms)
        self.psi_angles = dihedral_calc(self.psi_atoms)
        self.angles.append((self.phi_angles,self.psi_angles))

    def _conclude(self):
        self.angles = np.array(self.angles)

    def plot(self):
        fig = plt.figure(figsize=(10,10))
        ax1 = plt.subplot(111)
        ax1.axis([-180,180,-180,180])
        ax1.axhline(0, color='k', lw=1)
        ax1.axvline(0, color='k', lw=1)
        plt.xticks(np.arange(-180,181,60))
        plt.yticks(np.arange(-180,181,60))
        plt.xlabel(r"$\phi$ (deg)")
        plt.ylabel(r"$\psi$ (deg)")
        for angles in self.angles:
            ax1.plot(angles[0],angles[1],'ks')
