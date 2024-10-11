# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Nucleic acid analysis --- :mod:`MDAnalysis.analysis.nuclinfo`
=============================================================

:Author: Elizabeth Denning
:Year: 2011
:Copyright: GNU Public License v3

The module provides functions to analyze nucleic acid structures, in
particular

- backbone dihedrals,
- chi dihedrals,
- AS or CP phase angles,
- Watson-Crick N1-N3 distances, C2-O2 distances, N6-O4 distances, O6-N4 distances.

For applications of this kind of analysis see
:footcite:p:`Denning2011,Denning2012`.

All functions take a :class:`~MDAnalysis.core.universe.Universe` as an
argument together with further parameters that specify the base or bases in
question. Angles are in degrees. The functions use standard CHARMM names for
nucleic acids and atom names.


.. rubric:: References

.. footbibliography::

Distances
---------

.. autofunction:: wc_pair

.. autofunction:: minor_pair

.. autofunction:: major_pair


Phases
------

.. autofunction:: phase_cp

.. autofunction:: phase_as


Dihedral angles
---------------

.. autofunction:: tors

.. autofunction:: tors_alpha

.. autofunction:: tors_beta

.. autofunction:: tors_gamma

.. autofunction:: tors_delta

.. autofunction:: tors_eps

.. autofunction:: tors_zeta

.. autofunction:: tors_chi

.. autofunction:: hydroxyl

.. autofunction:: pseudo_dihe_baseflip

"""
import numpy as np
from math import pi, sin, cos, atan2, sqrt, pow

from MDAnalysis.lib import mdamath


def wc_pair(universe, i, bp, seg1="SYSTEM", seg2="SYSTEM"):
    """Watson-Crick basepair distance for residue `i` with residue `bp`.

    The distance of the nitrogen atoms in a Watson-Crick hydrogen bond is
    computed.

    Parameters
    ----------
    universe : Universe
         :class:`~MDAnalysis.core.universe.Universe` containing the trajectory
    i : int
        resid of the first base
    bp : int
        resid of the second base
    seg1 : str (optional)
        segment id for first base ["SYSTEM"]
    seg2 : str (optional)
        segment id for second base ["SYSTEM"]

    Returns
    -------
    float
        Watson-Crick base pair distance


    Notes
    -----
    If failure occurs be sure to check the segment identification.


    .. versionadded:: 0.7.6
    """
    if universe.select_atoms(" resid {0!s} ".format(i)).resnames[0] in ["DC", "DT", "U", "C", "T", "CYT", "THY", "URA"]:
        a1, a2 = "N3", "N1"
    if universe.select_atoms(" resid {0!s} ".format(i)).resnames[0] in ["DG", "DA", "A", "G", "ADE", "GUA"]:
        a1, a2 = "N1", "N3"
    wc_dist = universe.select_atoms("(segid {0!s} and resid {1!s} and name {2!s}) "
                                    "or (segid {3!s} and resid {4!s} and name {5!s}) "
                                    .format(seg1, i, a1, seg2, bp, a2))
    wc = mdamath.norm(wc_dist[0].position - wc_dist[1].position)
    return wc


def minor_pair(universe, i, bp, seg1="SYSTEM", seg2="SYSTEM"):
    """Minor-Groove basepair distance for residue `i` with residue `bp`.

    The distance of the nitrogen and oxygen atoms in a Minor-groove hydrogen
    bond is computed.

    Parameters
    ----------
    universe : Universe
         :class:`~MDAnalysis.core.universe.Universe` containing the trajectory
    i : int
        resid of the first base
    bp : int
        resid of the second base
    seg1 : str (optional)
        segment id for first base ["SYSTEM"]
    seg2 : str (optional)
        segment id for second base ["SYSTEM"]

    Returns
    -------
    float
        Minor groove base pair distance

    Notes
    -----
    If failure occurs be sure to check the segment identification.


    .. versionadded:: 0.7.6
    """
    if universe.select_atoms(" resid {0!s} ".format(i)).resnames[0] in ["DC", "DT", "U", "C", "T", "CYT", "THY", "URA"]:
        a1, a2 = "O2", "C2"
    if universe.select_atoms(" resid {0!s} ".format(i)).resnames[0] in ["DG", "DA", "A", "G", "ADE", "GUA"]:
        a1, a2 = "C2", "O2"
    c2o2_dist = universe.select_atoms("(segid {0!s} and resid {1!s} and name {2!s}) "
                                      "or (segid {3!s} and resid {4!s} and name {5!s})"
                                      .format(seg1, i, a1, seg2, bp, a2))
    c2o2 = mdamath.norm(c2o2_dist[0].position - c2o2_dist[1].position)
    return c2o2


def major_pair(universe, i, bp, seg1="SYSTEM", seg2="SYSTEM"):
    """Major-Groove basepair distance for residue `i` with residue `bp`.

    The distance of the nitrogen and oxygen atoms in a Major-groove hydrogen
    bond is computed.


    Parameters
    ----------
    universe : Universe
         :class:`~MDAnalysis.core.universe.Universe` containing the trajectory
    i : int
        resid of the first base
    bp : int
        resid of the second base
    seg1 : str (optional)
        segment id for first base ["SYSTEM"]
    seg2 : str (optional)
        segment id for second base ["SYSTEM"]

    Returns
    -------
    float
        Major groove base pair distance

    Notes
    -----
    If failure occurs be sure to check the segment identification.


    .. versionadded:: 0.7.6
    """
    if universe.select_atoms(" resid {0!s} ".format(i)).resnames[0] in ["DC", "DG", "C", "G", "CYT", "GUA"]:
        if universe.select_atoms(" resid {0!s} ".format(i)).resnames[0] in ["DC", "C", "CYT"]:
            a1, a2 = "N4", "O6"
        else:
            a1, a2 = "O6", "N4"
    if universe.select_atoms(" resid {0!s} ".format(i)).resnames[0] in ["DT", "DA", "A", "T", "U", "ADE", "THY", "URA"]:
        if universe.select_atoms(" resid {0!s} ".format(i)).resnames[0] in ["DT", "T", "THY", "U", "URA"]:
            a1, a2 = "O4", "N6"
        else:
            a1, a2 = "N6", "O4"
    no_dist = universe.select_atoms("(segid {0!s} and resid {1!s} and name {2!s}) "
                                    "or (segid {3!s} and resid {4!s} and name {5!s}) "
                                    .format(seg1, i, a1, seg2, bp, a2))
    major = mdamath.norm(no_dist[0].position - no_dist[1].position)
    return major


def phase_cp(universe, seg, i):
    """
    Pseudo-angle describing the phase of the ribose pucker for residue `i` using the CP method.
    The angle is computed by the positions of atoms in the ribose ring.

    Parameters
    ----------
    universe : Universe
        :class:`~MDAnalysis.core.universe.Universe` containing the trajectory
    seg : str
        Segment id for base
    i : int
        Resid of the first base

    Returns
    -------
    float
        Phase angle in degrees
    """
    
    # Select atoms
    atom1 = universe.select_atoms(f"atom {seg} {i} O4\'")
    atom2 = universe.select_atoms(f"atom {seg} {i} C1\'")
    atom3 = universe.select_atoms(f"atom {seg} {i} C2\'")
    atom4 = universe.select_atoms(f"atom {seg} {i} C3\'")
    atom5 = universe.select_atoms(f"atom {seg} {i} C4\'")

    # Get positions
    data1, data2, data3, data4, data5 = [atom.positions for atom in [atom1, atom2, atom3, atom4, atom5]]

    # Calculate mean position (r0)
    r0 = (data1 + data2 + data3 + data4 + data5) / 5.0

    # Calculate relative positions (r1, r2, r3, r4, r5)
    r1, r2, r3, r4, r5 = [data - r0 for data in [data1, data2, data3, data4, data5]]

    # Calculate R1 and R2
    R1 = (r1 * sin(2 * pi * 0.0 / 5.0)) + (r2 * sin(2 * pi * 1.0 / 5.0)) + \
         (r3 * sin(2 * pi * 2.0 / 5.0)) + (r4 * sin(2 * pi * 3.0 / 5.0)) + \
         (r5 * sin(2 * pi * 4.0 / 5.0))

    R2 = (r1 * cos(2 * pi * 0.0 / 5.0)) + (r2 * cos(2 * pi * 1.0 / 5.0)) + \
         (r3 * cos(2 * pi * 2.0 / 5.0)) + (r4 * cos(2 * pi * 3.0 / 5.0)) + \
         (r5 * cos(2 * pi * 4.0 / 5.0))

    # Normalize vector x and calculate dot products
    x = np.cross(R1[0], R2[0])
    n = x / np.linalg.norm(x)
    
    r_d = [np.dot(r, n) for r in [r1, r2, r3, r4, r5]]

    # Calculate D and C components
    D = sum(r_d[j] * sin(4 * pi * j / 5.0) for j in range(5)) * -1 * sqrt(2.0 / 5.0)
    C = sum(r_d[j] * cos(4 * pi * j / 5.0) for j in range(5)) * sqrt(2.0 / 5.0)

    # Calculate phase angle
    phase_ang = (np.arctan2(D, C) + (np.pi / 2.)) * 180. / np.pi
    return phase_ang % 360

def phase_as(universe, seg, i):
    """Pseudo-angle describing the phase of the ribose pucker for residue `i` using the AS method

    The angle is computed by the position vector of atoms in the ribose ring.

    Parameters
    ----------
    universe : Universe
         :class:`~MDAnalysis.core.universe.Universe` containing the trajectory
    seg : str
        segment id for base
    i : int
        resid of the first base

    Returns
    -------
    float
        phase angle in degrees


    .. versionadded:: 0.7.6
    """
    angle1 = universe.select_atoms(" atom {0!s} {1!s} C1\' ".format(seg, i),
                                   " atom {0!s} {1!s} C2\' ".format(seg, i),
                                   " atom {0!s} {1!s} C3\' ".format(seg, i),
                                   " atom {0!s} {1!s} C4\' ".format(seg, i))

    angle2 = universe.select_atoms(" atom {0!s} {1!s} C2\' ".format(seg, i),
                                   " atom {0!s} {1!s} C3\' ".format(seg, i),
                                   " atom {0!s} {1!s} C4\' ".format(seg, i),
                                   " atom {0!s} {1!s} O4\' ".format(seg, i))

    angle3 = universe.select_atoms(" atom {0!s} {1!s} C3\' ".format(seg, i),
                                   " atom {0!s} {1!s} C4\' ".format(seg, i),
                                   " atom {0!s} {1!s} O4\' ".format(seg, i),
                                   " atom {0!s} {1!s} C1\' ".format(seg, i))

    angle4 = universe.select_atoms(" atom {0!s} {1!s} C4\' ".format(seg, i),
                                   " atom {0!s} {1!s} O4\' ".format(seg, i),
                                   " atom {0!s} {1!s} C1\' ".format(seg, i),
                                   " atom {0!s} {1!s} C2\' ".format(seg, i))

    angle5 = universe.select_atoms(" atom {0!s} {1!s} O4\' ".format(seg, i),
                                   " atom {0!s} {1!s} C1\' ".format(seg, i),
                                   " atom {0!s} {1!s} C2\' ".format(seg, i),
                                   " atom {0!s} {1!s} C3\' ".format(seg, i))

    data1 = angle1.dihedral.value()
    data2 = angle2.dihedral.value()
    data3 = angle3.dihedral.value()
    data4 = angle4.dihedral.value()
    data5 = angle5.dihedral.value()

    B = ((data1 * sin(2 * 2 * pi * (1 - 1.) / 5.))
         + (data2 * sin(2 * 2 * pi * (2 - 1.) / 5.))
         + (data3 * sin(2 * 2 * pi * (3 - 1.) / 5.))
         + (data4 * sin(2 * 2 * pi * (4 - 1.) / 5.))
         + (data5 * sin(2 * 2 * pi * (5 - 1.) / 5.))) * -2. / 5.

    A = ((data1 * cos(2 * 2 * pi * (1 - 1.) / 5.))
         + (data2 * cos(2 * 2 * pi * (2 - 1.) / 5.))
         + (data3 * cos(2 * 2 * pi * (3 - 1.) / 5.))
         + (data4 * cos(2 * 2 * pi * (4 - 1.) / 5.))
         + (data5 * cos(2 * 2 * pi * (5 - 1.) / 5.))) * 2. / 5.

    phase_ang = (np.arctan2(D, C) + (np.pi / 2.)) * 180. / np.pi
    return phase_ang % 360


def tors(universe, seg, i):
    """Calculation of nucleic backbone dihedral angles.

    The dihedral angles are alpha, beta, gamma, delta, epsilon, zeta, chi.

    The dihedral is computed based on position of atoms for resid `i`.

    Parameters
    ----------
    universe : Universe
         :class:`~MDAnalysis.core.universe.Universe` containing the trajectory
    seg : str
        segment id for base
    i : int
        resid of the first base

    Returns
    -------
    [alpha, beta, gamma, delta, epsilon, zeta, chi] : list of floats
        torsion angles in degrees

    Notes
    -----
    If failure occurs be sure to check the segment identification.


    .. versionadded:: 0.7.6

    """
    a = universe.select_atoms(" atom {0!s} {1!s} O3\' ".format(seg, i - 1),
                              " atom {0!s} {1!s} P  ".format(seg, i),
                              " atom {0!s} {1!s} O5\' ".format(seg, i),
                              " atom {0!s} {1!s} C5\' ".format(seg, i))

    b = universe.select_atoms(" atom {0!s} {1!s} P    ".format(seg, i),
                              " atom {0!s} {1!s} O5\' ".format(seg, i),
                              " atom {0!s} {1!s} C5\' ".format(seg, i),
                              " atom {0!s} {1!s} C4\' ".format(seg, i))

    g = universe.select_atoms(" atom {0!s} {1!s} O5\' ".format(seg, i),
                              " atom {0!s} {1!s} C5\' ".format(seg, i),
                              " atom {0!s} {1!s} C4\' ".format(seg, i),
                              " atom {0!s} {1!s} C3\' ".format(seg, i))

    d = universe.select_atoms(" atom {0!s} {1!s} C5\' ".format(seg, i),
                              " atom {0!s} {1!s} C4\' ".format(seg, i),
                              " atom {0!s} {1!s} C3\' ".format(seg, i),
                              " atom {0!s} {1!s} O3\' ".format(seg, i))

    e = universe.select_atoms(" atom {0!s} {1!s} C4\' ".format(seg, i),
                              " atom {0!s} {1!s} C3\' ".format(seg, i),
                              " atom {0!s} {1!s} O3\' ".format(seg, i),
                              " atom {0!s} {1!s} P    ".format(seg, i + 1))

    z = universe.select_atoms(" atom {0!s} {1!s} C3\' ".format(seg, i),
                              " atom {0!s} {1!s} O3\' ".format(seg, i),
                              " atom {0!s} {1!s} P    ".format(seg, i + 1),
                              " atom {0!s} {1!s} O5\' ".format(seg, i + 1))
    c = universe.select_atoms(" atom {0!s} {1!s} O4\' ".format(seg, i),
                              " atom {0!s} {1!s} C1\' ".format(seg, i),
                              " atom {0!s} {1!s} N9 ".format(seg, i),
                              " atom {0!s} {1!s} C4  ".format(seg, i))
    if len(c) < 4:
        c = universe.select_atoms(" atom {0!s} {1!s} O4\' ".format(seg, i),
                                  " atom {0!s} {1!s} C1\' ".format(seg, i),
                                  " atom {0!s} {1!s} N1 ".format(seg, i),
                                  " atom {0!s} {1!s} C2  ".format(seg, i))

    alpha = a.dihedral.value() % 360
    beta = b.dihedral.value() % 360
    gamma = g.dihedral.value() % 360
    delta = d.dihedral.value() % 360
    epsilon = e.dihedral.value() % 360
    zeta = z.dihedral.value() % 360
    chi = c.dihedral.value() % 360

    return [alpha, beta, gamma, delta, epsilon, zeta, chi]


def tors_alpha(universe, seg, i):
    """alpha backbone dihedral

    The dihedral is computed based on position atoms for resid `i`.

    Parameters
    ----------
    universe : Universe
         :class:`~MDAnalysis.core.universe.Universe` containing the trajectory
    seg : str
        segment id for base
    i : int
        resid of the first base

    Returns
    -------
    alpha : float
        torsion angle in degrees


    .. versionadded:: 0.7.6
    """
    a = universe.select_atoms(" atom {0!s} {1!s} O3\' ".format(seg, i - 1),
                              " atom {0!s} {1!s} P  ".format(seg, i),
                              " atom {0!s} {1!s} O5\' ".format(seg, i),
                              " atom {0!s} {1!s} C5\' ".format(seg, i))
    alpha = a.dihedral.value() % 360
    return alpha


def tors_beta(universe, seg, i):
    """beta backbone dihedral

    The dihedral is computed based on position atoms for resid `i`.

    Parameters
    ----------
    universe : Universe
         :class:`~MDAnalysis.core.universe.Universe` containing the trajectory
    seg : str
        segment id for base
    i : int
        resid of the first base

    Returns
    -------
    beta : float
        torsion angle in degrees

    .. versionadded:: 0.7.6
    """
    b = universe.select_atoms(" atom {0!s} {1!s} P    ".format(seg, i),
                              " atom {0!s} {1!s} O5\' ".format(seg, i),
                              " atom {0!s} {1!s} C5\' ".format(seg, i),
                              " atom {0!s} {1!s} C4\' ".format(seg, i))
    beta = b.dihedral.value() % 360
    return beta


def tors_gamma(universe, seg, i):
    """ Gamma backbone dihedral

    The dihedral is computed based on position atoms for resid `i`.

    Parameters
    ----------
    universe : Universe
         :class:`~MDAnalysis.core.universe.Universe` containing the trajectory
    seg : str
        segment id for base
    i : int
        resid of the first base

    Returns
    -------
    gamma : float
        torsion angle in degrees


    .. versionadded:: 0.7.6
    """
    g = universe.select_atoms(" atom {0!s} {1!s} O5\' ".format(seg, i),
                              " atom {0!s} {1!s} C5\' ".format(seg, i),
                              " atom {0!s} {1!s} C4\' ".format(seg, i),
                              " atom {0!s} {1!s} C3\' ".format(seg, i))
    gamma = g.dihedral.value() % 360
    return gamma


def tors_delta(universe, seg, i):
    """delta backbone dihedral

    The dihedral is computed based on position atoms for resid `i`.

    Parameters
    ----------
    universe : Universe
         :class:`~MDAnalysis.core.universe.Universe` containing the trajectory
    seg : str
        segment id for base
    i : int
        resid of the first base

    Returns
    -------
    delta : float
        torsion angle in degrees


    .. versionadded:: 0.7.6
    """
    d = universe.select_atoms(" atom {0!s} {1!s} C5\' ".format(seg, i),
                              " atom {0!s} {1!s} C4\' ".format(seg, i),
                              " atom {0!s} {1!s} C3\' ".format(seg, i),
                              " atom {0!s} {1!s} O3\' ".format(seg, i))
    delta = d.dihedral.value() % 360
    return delta


def tors_eps(universe, seg, i):
    """Epsilon backbone dihedral

    The dihedral is computed based on position atoms for resid `i`.

    Parameters
    ----------
    universe : Universe
         :class:`~MDAnalysis.core.universe.Universe` containing the trajectory
    seg : str
        segment id for base
    i : int
        resid of the first base

    Returns
    -------
    epsilon : float
        torsion angle in degrees


    .. versionadded:: 0.7.6
    """
    e = universe.select_atoms(" atom {0!s} {1!s} C4\' ".format(seg, i),
                              " atom {0!s} {1!s} C3\' ".format(seg, i),
                              " atom {0!s} {1!s} O3\' ".format(seg, i),
                              " atom {0!s} {1!s} P    ".format(seg, i + 1))
    epsilon = e.dihedral.value() % 360
    return epsilon


def tors_zeta(universe, seg, i):
    """Zeta backbone dihedral

    The dihedral is computed based on position atoms for resid `i`.

    Parameters
    ----------
    universe : Universe
         :class:`~MDAnalysis.core.universe.Universe` containing the trajectory
    seg : str
        segment id for base
    i : int
        resid of the first base

    Returns
    -------
    zeta : float
        torsion angle in degrees


    .. versionadded:: 0.7.6
    """
    z = universe.select_atoms(" atom {0!s} {1!s} C3\' ".format(seg, i),
                              " atom {0!s} {1!s} O3\' ".format(seg, i),
                              " atom {0!s} {1!s} P    ".format(seg, i + 1),
                              " atom {0!s} {1!s} O5\' ".format(seg, i + 1))
    zeta = z.dihedral.value() % 360
    return zeta


def tors_chi(universe, seg, i):
    """chi nucleic acid dihedral

     The dihedral is computed based on position atoms for resid `i`.

    Parameters
    ----------
    universe : Universe
         :class:`~MDAnalysis.core.universe.Universe` containing the trajectory
    seg : str
        segment id for base
    i : int
        resid of the first base

    Returns
    -------
    chi : float
        torsion angle in degrees


    .. versionadded:: 0.7.6
    """
    c = universe.select_atoms(" atom {0!s} {1!s} O4\' ".format(seg, i),
                              " atom {0!s} {1!s} C1\' ".format(seg, i),
                              " atom {0!s} {1!s} N9 ".format(seg, i),
                              " atom {0!s} {1!s} C4  ".format(seg, i))
    if len(c) < 4:
        c = universe.select_atoms(" atom {0!s} {1!s} O4\' ".format(seg, i),
                                  " atom {0!s} {1!s} C1\' ".format(seg, i),
                                  " atom {0!s} {1!s} N1 ".format(seg, i),
                                  " atom {0!s} {1!s} C2  ".format(seg, i))
    chi = c.dihedral.value() % 360
    return chi


def hydroxyl(universe, seg, i):
    """2-hydroxyl dihedral. Useful only for RNA calculations.

     .. Note:: This dihedral calculation will only work if using atom
               names as documented by charmm force field parameters,
               namely "C1', C2', O2', H2'".

    Parameters
    ----------
    universe : Universe
         :class:`~MDAnalysis.core.universe.Universe` containing the trajectory
    seg : str
        segment id for base
    i : int
        resid of the first base

    Returns
    -------
    hydroxyl_angle : float
        torsion angle in degrees


    .. versionadded:: 0.7.6

    """
    h = universe.select_atoms("atom {0!s} {1!s} C1'".format(seg, i),
                              "atom {0!s} {1!s} C2'".format(seg, i),
                              "atom {0!s} {1!s} O2'".format(seg, i),
                              "atom {0!s} {1!s} H2'".format(seg, i))
    try:
        hydr = h.dihedral.value() % 360
    except ValueError:
        errmsg = (f"Resid {i} does not contain atoms C1', C2', O2', H2' but "
                  f"atoms {list(h.atoms)}")
        raise ValueError(errmsg) from None

    return hydr


def pseudo_dihe_baseflip(universe, bp1, bp2, i,
                         seg1="SYSTEM", seg2="SYSTEM", seg3="SYSTEM"):
    """pseudo dihedral for flipped bases. Useful only for nucleic acid base flipping

    The dihedral is computed based on position atoms for resid `i`

    .. Note:: This dihedral calculation will only work if using atom names as
              documented by charmm force field parameters.

    Parameters
    ----------
    universe : Universe
        :class:`~MDAnalysis.core.universe.Universe` containing the
        trajectory
    bp1 : int
        resid that base pairs with `bp2`
    bp2 : int
        resid below the base that flips
    i : int
        resid of the base that flips
    segid1 : str (optional)
        segid of resid base pairing with `bp2`
    segid2 : str (optional)
        segid, same as that of segid of flipping resid `i`
    segid3 : str (optional)
        segid of resid `i` that flips

    Returns
    -------
    float
          pseudo dihedral angle in degrees


    .. versionadded:: 0.8.0
    """
    bf1 = universe.select_atoms(
        " ( segid {0!s} and resid {1!s} and nucleicbase ) "
        "or ( segid {2!s} and resid {3!s} and nucleicbase ) "
        .format( seg1, bp1, seg2, bp2))
    bf4 = universe.select_atoms("(segid {0!s} and resid {1!s} and nucleicbase) ".format(seg3, i))
    bf2 = universe.select_atoms("(segid {0!s} and resid {1!s} and nucleicsugar) ".format(seg2, bp2))
    bf3 = universe.select_atoms("(segid {0!s} and resid {1!s} and nucleicsugar) ".format(seg3, i))
    x = [bf1.center_of_mass(), bf2.center_of_mass(),
         bf3.center_of_mass(), bf4.center_of_mass()]
    pseudo = mdamath.dihedral(x[0] - x[1], x[1] - x[2], x[2] - x[3])
    pseudo = np.rad2deg(pseudo) % 360
    return pseudo
