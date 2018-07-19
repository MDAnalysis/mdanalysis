import numpy as np
import matplotlib.pyplot as plt
import warnings
import six 
from six.moves import zip

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import calc_dihedrals


class Ramachandran(AnalysisBase):
    """Calculate phi and psi dihedral angles of specified residues.

    Note
    ----
    Run the analysis with :meth:`Ramachandran.run()`, which stores the results
    in the array :attr:`Ramachandran.angles`. A axes object can be obtained
    with :meth: `Ramachandran.run().plot()`.

    If the residue selection is beyond the scope of the protein, then an error
    will be raised. If the residue selection includes the first or last residue
    then a warning will be raised, and the final array of angles will not
    include those residues.

    """
    def __init__(self, atomgroup, **kwargs):
        r"""Parameters
        ----------
        atomgroup : Atomgroup
            atoms for residues for which phi and psi are calculated
        start : int, optional
            starting frame, default None becomes 0.
        stop : int, optional
            Frame index to stop analysis. Default: None becomes
            n_frames. Iteration stops *before* this frame number,
            which means that the trajectory would be read until the end.
        step : int, optional
            step between frames, default None becomes 1.

        """
        super(Ramachandran, self).__init__(atomgroup.universe.trajectory, **kwargs)
        self.atomgroup = atomgroup
        residues = self.atomgroup.residues
        protein = self.atomgroup.universe.select_atoms("protein").residues

        if not residues.issubset(protein):
            raise IndexError("Found atoms outside of protein. Only atoms "
                             "inside of a 'protein' selection can be used to "
                             "calculate dihedrals.")
        elif not residues.isdisjoint(protein[[0, -1]]):
            warnings.warn("Cannot determine phi and psi angles for the first "
                          "or last residues")
            residues = residues.difference(protein[[0, -1]])

        phi_sel = [res.phi_selection() for res in residues]
        psi_sel = [res.psi_selection() for res in residues]
        self.ag1 = mda.AtomGroup([atoms[0] for atoms in phi_sel])
        self.ag2 = mda.AtomGroup([atoms[1] for atoms in phi_sel])
        self.ag3 = mda.AtomGroup([atoms[2] for atoms in phi_sel])
        self.ag4 = mda.AtomGroup([atoms[3] for atoms in phi_sel])
        self.ag5 = mda.AtomGroup([atoms[3] for atoms in psi_sel])

    def _prepare(self):
        self.angles = []

    def _single_frame(self):
        phi_angles = calc_dihedrals(self.ag1.positions, self.ag2.positions,
                                    self.ag3.positions, self.ag4.positions,
                                    box=self.ag1.dimensions)
        psi_angles = calc_dihedrals(self.ag2.positions, self.ag3.positions,
                                    self.ag4.positions, self.ag5.positions,
                                    box=self.ag1.dimensions)
        phi_psi = zip(phi_angles, psi_angles)
        self.angles.append(phi_psi)

    def _conclude(self):
        self.angles = np.rad2deg(np.array(self.angles))


    def plot(self, ax=None, **kwargs):
        """Plots data into standard ramachandran plot. Each time step in
        self.angles is plotted onto the same graph.

        Parameters
        ----------
        ax : :class:`matplotlib.axes.Axes`
              If no `ax` is supplied or set to ``None`` then the plot will
              be added to the current active axes.

        Returns
        -------
        ax : :class:`~matplotlib.axes.Axes`
             Axes with the plot, either `ax` or the current axes.

        """
        if ax is None:
            ax = plt.gca()
        ax.axis([-180,180,-180,180])
        ax.axhline(0, color='k', lw=1)
        ax.axvline(0, color='k', lw=1)
        ax.set(xticks=range(-180,181,60), yticks=range(-180,181,60),
               xlabel=r"$\phi$ (deg)", ylabel=r"$\psi$ (deg)")
        a = self.angles.reshape(np.prod(self.angles.shape[:2]), 2)
        ax.scatter(a[:,0], a[:,1])
