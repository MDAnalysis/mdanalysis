import numpy as np

from scipy.spatial import KDTree
from numpy import arange, pi, sin, cos, arccos
from .base import AnalysisBase


ATOMIC_RADI = {
    "H": 1.200,
    "HE": 1.400,
    "C": 1.700,
    "N": 1.550,
    "O": 1.520,
    "F": 1.470,
    "NA": 2.270,
    "MG": 1.730,
    "P": 1.800,
    "S": 1.800,
    "CL": 1.750,
    "K": 2.750,
    "CA": 2.310,
    "NI": 1.630,
    "CU": 1.400,
    "ZN": 1.390,
    "SE": 1.900,
    "BR": 1.850,
    "CD": 1.580,
    "I": 1.980,
    "HG": 1.550,
}


def atomic_radi(e):
    """van der Waals radii""" 
    DEFAULT_ATOMIC_RADI = 2. 
    return ATOMIC_RADI[e]  if e in ATOMIC_RADI  else DEFAULT_ATOMIC_RADI


class SASA(AnalysisBase):
    r""" Calculation of solvent accesible surface areas using Shrake-Rupley "rolling ball'
         algorithm 
    """
    def __init__(self, ag, n_points=100, probe_radius=1.40, **kwargs):
        """Parameters
        ----------
        ag : AtomGroup
            Atom group used for the calculation.
        n_points: int (optional)
            Number of points used for the estimation of the atom surface area
        probe_radius: double ( optional )
            Radius of the probe atom, by default watter radious
        """
        super(SASA, self).__init__(ag.universe.trajectory, **kwargs)
        self._ag = ag
        self._atom_neighbours = KDTree(ag.positions, 10)
        self._n_points = n_points
        self._sphere = self._compute_sphere()
        self._twice_max_radi = (max(list(ATOMIC_RADI.values())) + probe_radius) * 2
        self._probe_radius = probe_radius


    def _compute_sphere(self):
        """Fibonacci lattice to evently distribute points in the sphere."""
        golden_ratio = (1 + 5**0.5) / 2
        i = arange(0, self._n_points)
        theta = 2 * pi * i / golden_ratio
        phi = arccos(1 - 2 * (i + 0.5) / self._n_points)
        x_points, y_points, z_points = (
            cos(theta) * sin(phi),
            sin(theta) * sin(phi),
            cos(phi),
        )
        return np.transpose(np.array([x_points, y_points, z_points]))


    def _compute(self):
        asas = np.zeros(len(self._ag.atoms), dtype=np.double)

        for atom in self._ag.atoms:
            r_w = atomic_radi(atom.type)
            raddi = r_w + self._probe_radius
            # translate and transform the sphere to the centroid
            sphere = np.array(self._sphere, copy=True) * raddi + atom.position
            neighbours = self._atom_neighbours.query_ball_point(
                atom.position, self._twice_max_radi
            )
            kdt_sphere = KDTree(sphere, 10)
            intersect = set()
            for n in neighbours:
                if n == atom.index:
                    continue
                n_raddi = atomic_radi(self._ag.atoms[n].type) + self._probe_radius
                cut = kdt_sphere.query_ball_point(self._ag.atoms[n].position, n_raddi)
                intersect |= set(cut)

            points = self._n_points - len(intersect)
            total_surface_area = raddi * raddi * 4 * np.pi
            area_per_point = total_surface_area / self._n_points
            asas[atom.index] = points * area_per_point
        return asas


    def atoms(self):
        """ returns accessible surface area in A**2"""
        return self._compute()
