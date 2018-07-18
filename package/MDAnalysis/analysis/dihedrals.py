import numpy as np
import matplotlib.pyplot as plt
import warnings

from MDAnalysis.analysis.base import AnalysisBase

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

    def _prepare(self):
        self.residues = self.atomgroup.residues
        res_min = np.min(self.atomgroup.universe.select_atoms("protein").residues)
        res_max = np.max(self.atomgroup.universe.select_atoms("protein").residues)
        if any([(residue < res_min or residue > res_max) for residue in self.residues]):
            raise IndexError("Selection exceeds protein length")
        elif any([residue == (res_min or res_max) for residue in self.residues]):
            warnings.warn("Cannot determine phi and psi angles for the first or last residues")

        self.phi_atoms = [residue.phi_selection() for residue in self.residues
                          if res_min < residue < res_max]
        self.psi_atoms = [residue.psi_selection() for residue in self.residues
                          if res_min < residue < res_max]

        self.angles = []

    def _single_frame(self):
        self.phi_angles = dihedral_calc(self.phi_atoms)
        self.psi_angles = dihedral_calc(self.psi_atoms)
        self.angles.append((self.phi_angles,self.psi_angles))

    def _conclude(self):
        self.angles = np.array(self.angles)

    def plot(self, ax=None, color='k', marker='s', title=None):
        """Plots data into standard ramachandran plot. Each time step in
        self.angles is plotted onto the same graph.

        Parameters
        ----------
        ax : :class:`matplotlib.axes.Axes`
              If no `ax` is supplied or set to ``None`` then the plot will
              be added to the current active axes.
        color : string, optional
              Color used for markers in the plot; the default color is 'black'.
        marker : string, optional
               Marker used in plot; the default marker is 'square'.
        title : string, optional
              Title of axes object; the default `None` leaves plot without a
              title.

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
        ax.set_xticks(range(-180,181,60))
        ax.set_yticks(range(-180,181,60))
        ax.set_xlabel(r"$\phi$ (deg)")
        ax.set_ylabel(r"$\psi$ (deg)")
        if title is not None:
            ax.set_title(title)
        for angles in self.angles:
            ax.plot(angles[0],angles[1], color=color,
                    marker=marker, linestyle='')
