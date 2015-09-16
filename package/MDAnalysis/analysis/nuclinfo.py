# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

# MDAnalysis -- nucleic acid analysis
# Copyright (c) 2011 Elizabeth Denning <denniej0@gmail.com>

"""
Nucleic acid analysis --- :mod:`MDAnalysis.analysis.nuclinfo`
=============================================================

:Author: Elizabeth Denning
:Year: 2011
:Copyright: GNU Public License v3

The module provides functions to analyze nuclic acid structures, in
particular

- backbone dihedrals,
- chi dihedrals,
- AS or CP phase angles,
- Watson-Crick N1-N3 distances, C2-O2 distances, N6-O4 distances, O6-N4 distances.

For applications of this kind of analysis see [Denning2011]_ and [Denning2012]_.

All functions take a :class:`~MDAnalysis.core.AtomGroup.Universe` as an argument
together with further parameters that specify the base or bases in question. Angles are
in degrees. The functions use standard CHARMM names for nucleic acids and atom names.


.. rubric:: References

.. [Denning2011] E.J. Denning, U.D. Priyakumar, L. Nilsson, and A.D. Mackerell, Jr. Impact of
              2'-hydroxyl sampling on the conformational properties of RNA: update of the
              CHARMM all-atom additive force field for RNA. *J. Comput. Chem.* 32 (2011),
              1929--1943. doi: `10.1002/jcc.21777`_

.. [Denning2012] E.J. Denning and A.D. MacKerell, Jr. Intrinsic Contribution of the 2'-Hydroxyl to
              RNA Conformational Heterogeneity. *J. Am. Chem. Soc.* 134 (2012), 2800--2806.
              doi: `10.1021/ja211328g`_


.. _`10.1002/jcc.21777`: http://dx.doi.org/10.1002/jcc.21777
.. _`10.1021/ja211328g`: http://dx.doi.org/10.1021/ja211328g


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
    """Watson-Crick basepair distance for residue *i* with residue *bp*.

    The distance of the nitrogen atoms in a Watson-Crick hydrogen bond is
    computed.

    :Arguments:
      *universe*
          :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory
      *seg1*
          segment id for first base
      *i*
          resid of the first base
      *seg2*
          segment id for second base
      *bp*
          resid of the second base

    .. NOTE:: If failure occurs be sure to check the segment identification.

    .. versionadded:: 0.7.6
    """
    if universe.select_atoms(" resid %s " % (i,)).resnames[0] in ["DC", "DT", "U", "C", "T", "CYT", "THY", "URA"]:
        a1, a2 = "N3", "N1"
    if universe.select_atoms(" resid %s " % (i,)).resnames[0] in ["DG", "DA", "A", "G", "ADE", "GUA"]:
        a1, a2 = "N1", "N3"
    wc_dist = universe.select_atoms(
        " (segid %s and resid %s and name %s)  or (segid %s and resid %s and name %s) " % (seg1, i, a1, seg2, bp, a2))
    wc = mdamath.norm(wc_dist[0].pos - wc_dist[1].pos)
    return wc


def minor_pair(universe, i, bp, seg1="SYSTEM", seg2="SYSTEM"):
    """Minor-Groove basepair distance for residue *i* with residue *bp*.

    The distance of the nitrogen and oxygen atoms in a Minor-groove hydrogen bond is
    computed.

    :Arguments:
      *universe*
          :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory
      *seg1*
          segment id for first base
      *i*
          resid of the first base
      *seg2*
          segment id for second base
      *bp*
          resid of the second base

    .. NOTE:: If failure occurs be sure to check the segment identification.

    .. versionadded:: 0.7.6
    """
    if universe.select_atoms(" resid %s " % (i,)).resnames[0] in ["DC", "DT", "U", "C", "T", "CYT", "THY", "URA"]:
        a1, a2 = "O2", "C2"
    if universe.select_atoms(" resid %s " % (i,)).resnames[0] in ["DG", "DA", "A", "G", "ADE", "GUA"]:
        a1, a2 = "C2", "O2"
    c2o2_dist = universe.select_atoms(
        " (segid %s and resid %s and name %s)  or (segid %s and resid %s and name %s) " % (seg1, i, a1, seg2, bp, a2))
    c2o2 = mdamath.norm(c2o2_dist[0].pos - c2o2_dist[1].pos)
    return c2o2


def major_pair(universe, i, bp, seg1="SYSTEM", seg2="SYSTEM"):
    """Major-Groove basepair distance for residue *i* with residue *bp*.

    The distance of the nitrogen and oxygen atoms in a Major-groove hydrogen bond is
    computed.

    :Arguments:
      *universe*
          :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory       *i*
      *seg1*
          segment id for first base
      *i*
          resid of the first base
      *seg2*
          segment id for second base
      *bp*
          resid of the second base

    .. NOTE:: If failure occurs be sure to check the segment identification

    .. versionadded:: 0.7.6
    """
    if universe.select_atoms(" resid %s " % (i,)).resnames[0] in ["DC", "DG", "C", "G", "CYT", "GUA"]:
        if universe.select_atoms(" resid %s " % (i,)).resnames[0] in ["DC", "C", "CYT"]:
            a1, a2 = "N4", "O6"
        else:
            a1, a2 = "O6", "N4"
    if universe.select_atoms(" resid %s " % (i,)).resnames[0] in ["DT", "DA", "A", "T", "U", "ADE", "THY", "URA"]:
        if universe.select_atoms(" resid %s " % (i,)).resnames[0] in ["DT", "T", "THY", "U", "URA"]:
            a1, a2 = "O4", "N6"
        else:
            a1, a2 = "N6", "O4"
    no_dist = universe.select_atoms(
        " (segid %s and resid %s and name %s)  or (segid %s and resid %s and name %s) " % (seg1, i, a1, seg2, bp, a2))
    major = mdamath.norm(no_dist[0].pos - no_dist[1].pos)
    return major


def phase_cp(universe, seg, i):
    """Pseudo-angle describing the phase of the ribose pucker for residue *i* using the CP method.

    The angle is computed by the positions of atoms in the ribose ring.

    :Arguments:
      *universe*
         :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory

       *segid*
         segment identity of resid
      *i*
         resid of the base

    .. versionadded:: 0.7.6
    """
    atom1 = universe.select_atoms(" atom %s %s O4\' " % (seg, i))
    atom2 = universe.select_atoms(" atom %s %s C1\' " % (seg, i))
    atom3 = universe.select_atoms(" atom %s %s C2\' " % (seg, i))
    atom4 = universe.select_atoms(" atom %s %s C3\' " % (seg, i))
    atom5 = universe.select_atoms(" atom %s %s C4\' " % (seg, i))

    data1 = atom1.coordinates()
    data2 = atom2.coordinates()
    data3 = atom3.coordinates()
    data4 = atom4.coordinates()
    data5 = atom5.coordinates()

    r0 = (data1 + data2 + data3 + data4 + data5) * (1.0 / 5.0)
    r1 = data1 - r0
    r2 = data2 - r0
    r3 = data3 - r0
    r4 = data4 - r0
    r5 = data5 - r0

    R1 = ((r1 * sin(2 * pi * 0.0 / 5.0)) + (r2 * sin(2 * pi * 1.0 / 5.0)) +
          (r3 * sin(2 * pi * 2.0 / 5.0)) + (r4 * sin(2 * pi * 3.0 / 5.0)) + (r5 * sin(2 * pi * 4.0 / 5.0)))

    R2 = ((r1 * cos(2 * pi * 0.0 / 5.0)) + (r2 * cos(2 * pi * 1.0 / 5.0)) +
          (r3 * cos(2 * pi * 2.0 / 5.0)) + (r4 * cos(2 * pi * 3.0 / 5.0)) + (r5 * cos(2 * pi * 4.0 / 5.0)))

    x = np.cross(R1[0], R2[0])
    n = x / sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2))

    r1_d = np.dot(r1, n)
    r2_d = np.dot(r2, n)
    r3_d = np.dot(r3, n)
    r4_d = np.dot(r4, n)
    r5_d = np.dot(r5, n)

    D = ((r1_d * sin(4 * pi * 0.0 / 5.0)) + (r2_d * sin(4 * pi * 1.0 / 5.0)) +
         (r3_d * sin(4 * pi * 2.0 / 5.0)) + (r4_d * sin(4 * pi * 3.0 / 5.0)) + (r5_d * sin(4 * pi * 4.0 / 5.0))) \
        * -1 * sqrt(2.0 / 5.0)

    C = ((r1_d * cos(4 * pi * 0.0 / 5.0)) + (r2_d * cos(4 * pi * 1.0 / 5.0)) + (r3_d * cos(4 * pi * 2.0 / 5.0)) +
         (r4_d * cos(4 * pi * 3.0 / 5.0)) + (r5_d * cos(4 * pi * 4.0 / 5.0))) * sqrt(2.0 / 5.0)

    phase_ang = (atan2(D, C) + (pi / 2.)) * 180. / pi
    if phase_ang < 0:
        phase_ang = phase_ang + 360
    else:
        phase_ang
    return phase_ang


def phase_as(universe, seg, i):
    """Pseudo-angle describing the phase of the ribose pucker for residue *i* using the AS method

    The angle is computed by the position vector of atoms in the ribose ring.

    :Arguments:
      *universe*
         :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory
      *segid*
         segment identity of resid
      *i*
         resid of the base

    .. versionadded:: 0.7.6
    """
    angle1 = universe.select_atoms(" atom %s %s C1\' " % (seg, i), " atom %s %s C2\' " % (seg, i),
                                  " atom %s %s C3\' " % (seg, i), " atom %s %s C4\' " % (seg, i))
    angle2 = universe.select_atoms(" atom %s %s C2\' " % (seg, i), " atom %s %s C3\' " % (seg, i),
                                  " atom %s %s C4\' " % (seg, i), " atom %s %s O4\' " % (seg, i))
    angle3 = universe.select_atoms(" atom %s %s C3\' " % (seg, i), " atom %s %s C4\' " % (seg, i),
                                  " atom %s %s O4\' " % (seg, i), " atom %s %s C1\' " % (seg, i))
    angle4 = universe.select_atoms(" atom %s %s C4\' " % (seg, i), " atom %s %s O4\' " % (seg, i),
                                  " atom %s %s C1\' " % (seg, i), " atom %s %s C2\' " % (seg, i))
    angle5 = universe.select_atoms(" atom %s %s O4\' " % (seg, i), " atom %s %s C1\' " % (seg, i),
                                  " atom %s %s C2\' " % (seg, i), " atom %s %s C3\' " % (seg, i))

    data1 = angle1.dihedral.value()
    data2 = angle2.dihedral.value()
    data3 = angle3.dihedral.value()
    data4 = angle4.dihedral.value()
    data5 = angle5.dihedral.value()

    B = ((data1 * sin(2 * 2 * pi * (1 - 1.) / 5.)) + (data2 * sin(2 * 2 * pi * (2 - 1.) / 5.)) +
         (data3 * sin(2 * 2 * pi * (3 - 1.) / 5.)) + (data4 * sin(2 * 2 * pi * (4 - 1.) / 5.)) +
         (data5 * sin(2 * 2 * pi * (5 - 1.) / 5.))) * -2. / 5.

    A = ((data1 * cos(2 * 2 * pi * (1 - 1.) / 5.)) + (data2 * cos(2 * 2 * pi * (2 - 1.) / 5.)) +
         (data3 * cos(2 * 2 * pi * (3 - 1.) / 5.)) + (data4 * cos(2 * 2 * pi * (4 - 1.) / 5.)) +
         (data5 * cos(2 * 2 * pi * (5 - 1.) / 5.))) * 2. / 5.

    phase_ang = atan2(B, A) * 180. / pi
    if phase_ang < 0:
        phase_ang = phase_ang + 360
    else:
        phase_ang
    return phase_ang


def tors(universe, seg, i):
    """Backbone dihedrals includes alpha, beta, gamma, delta, epsilon, zeta, chi

    The dihedral is computed based on position atoms for resid *i*.

    :Arguments:
      *universe*
         :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory
      *segid*
         segid of resid
      *i*
         resid of the base

    .. NOTE:: If failure occurs be sure to check the segment identification

    """
    a = universe.select_atoms(" atom %s %s O3\' " % (seg, i - 1), " atom %s %s P  " % (seg, i),
                             " atom %s %s O5\' " % (seg, i), " atom %s %s C5\' " % (seg, i))
    b = universe.select_atoms(" atom %s %s P    " % (seg, i), " atom %s %s O5\' " % (seg, i),
                             " atom %s %s C5\' " % (seg, i), " atom %s %s C4\' " % (seg, i))
    g = universe.select_atoms(" atom %s %s O5\' " % (seg, i), " atom %s %s C5\' " % (seg, i),
                             " atom %s %s C4\' " % (seg, i), " atom %s %s C3\' " % (seg, i))
    d = universe.select_atoms(" atom %s %s C5\' " % (seg, i), " atom %s %s C4\' " % (seg, i),
                             " atom %s %s C3\' " % (seg, i), " atom %s %s O3\' " % (seg, i))
    e = universe.select_atoms(" atom %s %s C4\' " % (seg, i), " atom %s %s C3\' " % (seg, i),
                             " atom %s %s O3\' " % (seg, i), " atom %s %s P    " % (seg, i + 1))
    z = universe.select_atoms(" atom %s %s C3\' " % (seg, i), " atom %s %s O3\' " % (seg, i),
                             " atom %s %s P    " % (seg, i + 1), " atom %s %s O5\' " % (seg, i + 1))
    try:
        c = universe.select_atoms(" atom %s %s O4\' " % (seg, i), " atom %s %s C1\' " % (seg, i),
                                 " atom %s %s N1 " % (seg, i), " atom %s %s C2  " % (seg, i))
    except:
        c = universe.select_atoms(" atom %s %s O4\' " % (seg, i), " atom %s %s C1\' " % (seg, i),
                                 " atom %s %s N9 " % (seg, i), " atom %s %s C4  " % (seg, i))

    alpha = a.dihedral.value()
    beta = b.dihedral.value()
    gamma = g.dihedral.value()
    delta = d.dihedral.value()
    epsilon = e.dihedral.value()
    zeta = z.dihedral.value()
    chi = c.dihedral.value()

    if alpha < 0:
        alpha = alpha + 360
    if beta < 0:
        beta = beta + 360
    if gamma < 0:
        gamma = gamma + 360
    if epsilon < 0:
        epsilon = epsilon + 360
    if zeta < 0:
        zeta = zeta + 360
    if chi < 0:
        chi = chi + 360
    return [alpha, beta, gamma, delta, epsilon, zeta, chi]


def tors_alpha(universe, seg, i):
    """alpha backbone dihedral

    The dihedral is computed based on position atoms for resid *i*.

    :Arguments:
      *universe*
         :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory
      *segid*
         segid of resid
      *i*
         resid of the base

    .. versionadded:: 0.7.6
    """
    a = universe.select_atoms(" atom %s %s O3\' " % (seg, i - 1), " atom %s %s P  " % (seg, i),
                             " atom %s %s O5\' " % (seg, i), " atom %s %s C5\' " % (seg, i))
    alpha = a.dihedral.value()
    if alpha < 0:
        alpha = alpha + 360
    return alpha


def tors_beta(universe, seg, i):
    """beta  backbone dihedral

    The dihedral is computed based on position atoms for resid *i*.

    :Arguments:
      *universe*
         :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory
      *segid*
         segid of resid
      *i*
         resid of the base

    .. versionadded:: 0.7.6
    """
    b = universe.select_atoms(" atom %s %s P    " % (seg, i), " atom %s %s O5\' " % (seg, i),
                             " atom %s %s C5\' " % (seg, i), " atom %s %s C4\' " % (seg, i))
    beta = b.dihedral.value()
    if beta < 0:
        beta = beta + 360
    return beta


def tors_gamma(universe, seg, i):
    """ Gamma backbone dihedral

     The dihedral is computed based on position atoms for resid *i*.

    :Arguments:
      *universe*
         :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory
      *segid*
         segid of resid
      *i*
         resid of the base

    .. versionadded:: 0.7.6
    """
    g = universe.select_atoms(" atom %s %s O5\' " % (seg, i), " atom %s %s C5\' " % (seg, i),
                             " atom %s %s C4\' " % (seg, i), " atom %s %s C3\' " % (seg, i))
    gamma = g.dihedral.value()
    if gamma < 0:
        gamma = gamma + 360
    return gamma


def tors_delta(universe, seg, i):
    """delta backbone dihedral

    The dihedral is computed based on position atoms for resid *i*.

    :Arguments:
      *universe*
         :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory
      *segid*
         segid of resid
      *i*
         resid of the base

    .. versionadded:: 0.7.6
    """
    d = universe.select_atoms(" atom %s %s C5\' " % (seg, i), " atom %s %s C4\' " % (seg, i),
                             " atom %s %s C3\' " % (seg, i), " atom %s %s O3\' " % (seg, i))
    delta = d.dihedral.value()
    if delta < 0:
        delta = delta + 360
    return delta


def tors_eps(universe, seg, i):
    """Epsilon backbone dihedral

    The dihedral is computed based on position atoms for resid *i*.

    :Arguments:
      *universe*
         :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory
      *segid*
         segid of resid
      *i*
         resid of the base

    .. versionadded:: 0.7.6
    """
    e = universe.select_atoms(" atom %s %s C4\' " % (seg, i), " atom %s %s C3\' " % (seg, i),
                             " atom %s %s O3\' " % (seg, i), " atom %s %s P    " % (seg, i + 1))
    epsilon = e.dihedral.value()
    if epsilon < 0:
        epsilon = epsilon + 360
    return epsilon


def tors_zeta(universe, seg, i):
    """Zeta backbone dihedral

    The dihedral is computed based on position atoms for resid *i*.

    :Arguments:
      *universe*
         :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory
      *segid*
         segid of resid
      *i*
         resid of the base

    .. versionadded:: 0.7.6
    """
    z = universe.select_atoms(" atom %s %s C3\' " % (seg, i), " atom %s %s O3\' " % (seg, i),
                             " atom %s %s P    " % (seg, i + 1), " atom %s %s O5\' " % (seg, i + 1))
    zeta = z.dihedral.value()
    if zeta < 0:
        zeta = zeta + 360
    return zeta


def tors_chi(universe, seg, i):
    """chi nucleic acid dihedral

     The dihedral is computed based on position atoms for resid *i*.

     :Arguments:
       *universe*
           :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory
       *segid*
           segid of resid
       *i*
           resid of the base

    .. versionadded:: 0.7.6
    """
    try:
        c = universe.select_atoms(" atom %s %s O4\' " % (seg, i), " atom %s %s C1\' " % (seg, i),
                                 " atom %s %s N1 " % (seg, i), " atom %s %s C2  " % (seg, i))
    except:
        c = universe.select_atoms(" atom %s %s O4\' " % (seg, i), " atom %s %s C1\' " % (seg, i),
                                 " atom %s %s N9 " % (seg, i), " atom %s %s C4  " % (seg, i))
    chi = c.dihedral.value()
    if chi < 0:
        chi = chi + 360
    return chi


def hydroxyl(universe, seg, i):
    """2-hydroxyl dihedral. Useful only for RNA calculations.

     .. Note:: This dihedral calculation will only work if using atom names as
               documented by charmm force field parameters.

     :Arguments:
       *universe*
           :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory
       *segid*
           segid of resid
       *i*
           resid of the base

    .. versionadded:: 0.7.6
    """
    h = universe.select_atoms(" atom %s %s C1\' " % (seg, i), " atom %s %s C2\' " % (seg, i),
                             " atom %s %s O2\' " % (seg, i), " atom %s %s H2\'\' " % (seg, i))
    try:
        hydr = h.dihedral.value()
    except ValueError:
        raise ValueError("Resid {0} does not contain atoms C1', C2', O2', H2' but atoms {1}".format(
                i, str(list(h.atoms))))
    if hydr < 0:
        hydr = hydr + 360
    return hydr


def pseudo_dihe_baseflip(universe, bp1, bp2, i, seg1="SYSTEM", seg2="SYSTEM", seg3="SYSTEM"):
    """pseudo dihedral for flipped bases. Useful only for nucleic acid base flipping

     The dihedral is computed based on position atoms for resid *i*

     .. Note:: This dihedral calculation will only work if using atom names as
               documented by charmm force field parameters.

     :Arguments:
       *universe*
           :class:`~MDAnalysis.core.AtomGroup.Universe` containing the trajectory
       *segid1*
           segid of resid base pairing with bp2
       *bp1*
           resid that base pairs with bp2
       *segid2*
           segid same as that of segid of flipping resid
       *bp2*
           resid below the base that flips
       *segid3*
           segid of resid that flips
       *i*
           resid of the base that flips

    .. versionadded:: 0.8.0
    """
    bf1 = universe.select_atoms(
        " ( segid %s and resid %s and nucleicbase ) or ( segid %s and resid %s and nucleicbase ) " % (
        seg1, bp1, seg2, bp2))
    bf4 = universe.select_atoms(" ( segid %s and resid %s and nucleicbase ) " % (seg3, i))
    bf2 = universe.select_atoms(" ( segid %s and resid %s and nucleicsugar ) " % (seg2, bp2))
    bf3 = universe.select_atoms(" ( segid %s and resid %s and nucleicsugar ) " % (seg3, i))
    x = [bf1.center_of_mass(), bf2.center_of_mass(), bf3.center_of_mass(), bf4.center_of_mass()]
    pseudo = mdamath.dihedral(x[0] - x[1], x[1] - x[2], x[2] - x[3])
    pseudo = np.rad2deg(pseudo)
    if pseudo < 0:
        pseudo = pseudo + 360
    return pseudo
