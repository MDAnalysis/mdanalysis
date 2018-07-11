from MDAnalysis.analysis.base import AnalysisBase
import numpy as np
import matplotlib.pyplot as plt


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
    """Calculate phi and psi dihedral angles of specified residues.

    Note
    ----
    Run the analysis with :meth:`Ramachandran.run()`, which stores the results
    in the array :attr:`Ramachandran.angles`. A basic plot can be obtained
    with :meth: `Ramachandran.run().plot()`.   

    """
    def __init__(self, atomgroup, **kwargs):
        r"""Parameters
        ----------
        atomgroup : Atomgroup
            atoms for residues for which phi and psi are calculated
        start : int (optional)
            starting frame, default None becomes 0.
        stop : int (optional)
            Frame index to stop analysis. Default: None becomes
            n_frames. Iteration stops *before* this frame number,
            which means that the trajectory would be read until the end.
        step : int (optional)
            step between frames, default None becomes 1.
        verbose : bool (optional)
             Show detailed progress of the calculation if set to ``True``; the
             default is ``False``.

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
