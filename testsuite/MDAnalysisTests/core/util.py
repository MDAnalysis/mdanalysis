# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
from __future__ import absolute_import, division
from six.moves import range

import numpy as np
import MDAnalysis as mda
from MDAnalysis.core import topologyattrs
from MDAnalysis.lib.mdamath import triclinic_box
from MDAnalysis.lib.distances import transform_RtoS


def assert_in_box(positions, box):
    """Asserts that all `positions` are strictly within the primary periodic
    image as defined by `box`
    """
    relpos = transform_RtoS(positions, box)
    assert np.all((relpos >= 0.0) & (relpos < 1.0))


class UnWrapUniverse(object):
    """A Universe containing small molecule clusters with molecules (as well as
    clusters) broken or whole accross periodic boundaries:

    The universe comprises 47 atoms in 12 molecules:
    * Three "molecules" (type A) containing a single atom:
     1. within the central image
     2. in a neighboring image
     3. in an image two boxes away from the central image

    * Four molecules (type B) containing three atoms:
     4. within the central image
     5. in the corner of the box, close to 4., broken accross one boundary
        but whole accross another
     6. outside the box, broken diagonally accross two neighboring images,
        close to 5.
     7. whole accross a boundary opposite to the location of 4., i.e.,
        close to 4. in terms of PBC
    * Two cyclical molecules (type C) containing four atoms:
     8. broken accross the front/back but whole accross the top face
     9. within the central image close to the front and bottom face, close
        to 8. in terms of PBC
    * Three linear chain molecules (type D) containing 8 atoms, spanning
     more than half a box length:
     10. within the central image relatively close to the top boundary
     11. close to 10, broken mid-molecule accross the left/right boundary
     12. close to 11. but in another image, whole accross the same
         boundary as 11. but not mid-molecule

    There are 15 residues in the universe:
    Molecules of type A, B, and C each have a single residue, while each of
    the chain molecules of type D have two residues with 4 atoms per residue.

    Atoms can be selected by their residue's resname. For molecules of type
    A, B, and C, the resnames are A, B, and C, respectively. Molecules of
    type D have the resnames D1 and D2.

    Atoms can also be selected by their moltype attribute, which is identical
    to the corresponding molecule type (A, B, C, D).

    The molecules/residues are contained in 6 segments, whereof the first
    three segments contain all molecules of type A, B, or C, respectively.
    Each of the remaining three segments contain a molecule of type D.

    A projection onto the x-z plane of box (orthorhombic case) looks roughly
    like this::

             :                   :
             :                   :
        6    :        8          :                   :
      - + - -:+-------+----------:- - - - - - - - - -: -
             |5       8          |                   :
             |                   |                   :
             |         (10)      |                   :
             |   o-o-o-o-x-x-x-x |                   :
             +x-x-x-x     o-o-o-o+                   :
             |            (11)   |                   :
             |                   |                   :
             |         1         |       2           : 3
             |                   |                   :
             |          9        |                   :
             |   4      !       7|                   :
             |   !      9       !|                   :
             |   4-4            7+7                  :
             |                   |                   :
            5+5                  |                   :
      - + - -:+------------------:- - - - - - - - - -: -
        6-6  :                   :                   :
             :                   :
             :                   :
             :                   :
             :             (12)  :
             :        x-x-x-x-o-o+o-o
             :                   :

    Note that the cyclic structures of molecules 8 and 9 lie in the x-y plane
    and are therefore not fully visible in the x-z plane projection.

    Parameters
    ----------
    have_bonds : bool, optional
        If ``True``, intramolecular bonds will be present in the topology.
        If ``False``, there will be no bond information.
    have_masses : bool, optional
        If ``True``, each atom will be assigned a mass of 1 u.
        If ``False``, masses won't be present in the topology.
    have_molnums : bool, optional
        If ``True``, the topology will have molecule information (molnums).
        If ``False``, that information will be missing.
    have_charges : bool, optional
        If ``False``, charges won't be present in the topology.
        If ``True``, atoms will carry the following charges::
        * atoms of molecule type A: 2 (total: 6)
        * atoms of molecule type B: -0.5 (total: -6)
        * atoms of molecule type C: -1.5 (total: -12)
        * atoms of molecule type D: 0.5 (total: 12)
    is_triclinic : bool, optional
        If ``False``, the box will be a cube with an edge length of 10 Angstrom.
        If ``True``, the cubic box will be sheared by 1 Angstrom in the x and
        y directions.
    """

    def __new__(cls, have_bonds=True, have_masses=True, have_molnums=True,
                have_charges=True, is_triclinic=False):
        # box:
        a = 10.0 # edge length
        tfac = 0.1  # factor for box vector shift of triclinic boxes (~84°)
        if is_triclinic:
            box = triclinic_box([a, 0.0, 0.0],
                                [a * tfac, a, 0.0],
                                [a * tfac, a * tfac, a])
        else:
            box = np.array([a, a, a, 90.0, 90.0, 90.0], dtype=np.float32)

        # number of atoms, residues, and segments:
        n_atoms = 47
        n_residues = 15
        n_segments = 6

        # resindices:
        residx = np.empty(n_atoms, dtype=np.int64)
        # type A
        rix = 0
        for i in range(0, 3, 1):
            residx[i] = rix
            rix += 1
        # type B
        for i in range(3, 15, 3):
            residx[i:i+3] = rix
            rix += 1
        # type C & D
        for i in range(15, 47, 4):
            residx[i:i+4] = rix
            rix += 1

        # segindices:
        segidx = np.empty(n_residues, dtype=np.int64)
        segidx[0:3] = 0
        segidx[3:7] = 1
        segidx[7:9] = 2
        segidx[9:11] = 3
        segidx[11:13] = 4
        segidx[13:15] = 5

        # universe:
        u = mda.Universe.empty(
            # topology things
            n_atoms=n_atoms,
            n_residues=n_residues,
            n_segments=n_segments,
            atom_resindex=residx,
            residue_segindex=segidx,
            # trajectory things
            trajectory=True,
            velocities=False,
            forces=False,
        )

        # resnames: we always want those for selection purposes
        resnames = ['A'] * 3
        resnames += ['B'] * 4
        resnames += ['C'] * 2
        resnames += ['D1', 'D2'] * 3
        u.add_TopologyAttr(topologyattrs.Resnames(resnames))

        # moltypes: we always want those for selection purposes
        moltypes = ['A'] * 3
        moltypes += ['B'] * 4
        moltypes += ['C'] * 2
        moltypes += ['D'] * 6
        u.add_TopologyAttr(topologyattrs.Moltypes(moltypes))

        # trajectory:
        ts = u.trajectory.ts
        ts.frame = 0
        ts.dimensions = box

        # positions:
        relpos = np.empty((n_atoms, 3), dtype=np.float32)
        # type A
        relpos[0:3, :] = np.array([[0.5, 0.5, 0.5],
                                   [1.4, 0.5, 0.5],
                                   [2.1, 0.5, 0.5]], dtype=np.float32)
        # type B
        relpos[3:15, :] = np.array([[0.1, 0.1, 0.2],
                                    [0.1, 0.1, 0.1],
                                    [0.2, 0.1, 0.1],
                                    [-0.05, 0.2, 0.05],
                                    [0.05, 0.2, 0.05],
                                    [0.05, 0.2, 0.95],
                                    [-0.2, -0.9, 1.05],
                                    [-0.2, 0.1, -0.05],
                                    [-0.1, 0.1, -0.05],
                                    [0.95, 0.2, 0.25],
                                    [0.95, 0.2, 0.15],
                                    [1.05, 0.2, 0.15]], dtype=np.float32)
        # type C
        relpos[15:23, :] = np.array([[0.4, 0.95, 1.05],
                                     [0.4, 0.95, 0.95],
                                     [0.4, 0.05, 0.95],
                                     [0.4, 0.05, 1.05],
                                     [0.6, 0.05, 0.25],
                                     [0.6, 0.05, 0.15],
                                     [0.6, 0.15, 0.15],
                                     [0.6, 0.15, 0.25]], dtype=np.float32)
        # type D
        relpos[23:47, :] = np.array([[0.2, 0.7, 0.8],
                                     [0.3, 0.7, 0.8],
                                     [0.4, 0.7, 0.8],
                                     [0.5, 0.7, 0.8],
                                     [0.6, 0.7, 0.8],
                                     [0.7, 0.7, 0.8],
                                     [0.8, 0.7, 0.8],
                                     [0.9, 0.7, 0.8],
                                     [0.66, 0.75, 0.7],
                                     [0.76, 0.75, 0.7],
                                     [0.86, 0.75, 0.7],
                                     [0.96, 0.75, 0.7],
                                     [0.06, 0.75, 0.7],
                                     [0.16, 0.75, 0.7],
                                     [0.26, 0.75, 0.7],
                                     [0.36, 0.75, 0.7],
                                     [1.14, 0.65, -0.4],
                                     [1.04, 0.65, -0.4],
                                     [0.94, 0.65, -0.4],
                                     [0.84, 0.65, -0.4],
                                     [0.74, 0.65, -0.4],
                                     [0.64, 0.65, -0.4],
                                     [0.54, 0.65, -0.4],
                                     [0.44, 0.65, -0.4]], dtype=np.float32)
        # make a copy, we need the original later
        _relpos = relpos.copy()
        # apply y- and z-dependent shift of x and y coords for triclinic boxes:
        if is_triclinic:
            # x-coord shift depends on y- and z-coords
            _relpos[:, 0] += tfac * _relpos[:, 1:].sum(axis=1)
            # y-coord shift depends on z-coords only
            _relpos[:, 1] += tfac * _relpos[:, 2]
        # scale relative to absolute positions:
        ts.positions = (_relpos * np.array([a, a, a])).astype(np.float32)

        # bonds:
        if have_bonds:
            bonds = []
            # type A has no bonds
            #type B
            for base in range(3, 15, 3):
                for i in range(2):
                    bonds.append((base + i, base + i + 1))
            # type C
            for base in range(15, 23, 4):
                for i in range(3):
                    bonds.append((base + i, base + i + 1))
                bonds.append((0 + base, 3 + base))
            # type D
            for base in range(23, 47, 8):
                for i in range(7):
                    bonds.append((base + i, base + i + 1))
            u.add_TopologyAttr(topologyattrs.Bonds(bonds))

        # masses:
        if have_masses:
            # masses are all set to 1 so that one can cross-check the results of
            # reference='com' with reference='cog' unwrapping
            masses = np.ones(n_atoms)
            u.add_TopologyAttr(topologyattrs.Masses(masses))

        # molnums:
        if have_molnums:
            molnums = [0, 1, 2]
            molnums += [3, 4, 5, 6]
            molnums += [7, 8]
            molnums += [9, 9, 10, 10, 11, 11]
            u.add_TopologyAttr(topologyattrs.Molnums(molnums))

        # charges:
        if have_charges:
            # type A
            charges = [2] * 3
            # type B
            charges += [-0.5] * 12
            # type C
            charges += [-1.5] * 8
            # type C
            charges += [0.5] * 24
            u.add_TopologyAttr(topologyattrs.Charges(charges))

        # shamelessly monkey-patch some custom universe attributes:
        u._is_triclinic = is_triclinic
        u._relpos = relpos
        u._tfac = tfac
        u._box_edge = a
        # bind custom methods to universe:
        u.unwrapped_coords = cls.unwrapped_coords.__get__(u)
        u.wrapped_coords = cls.wrapped_coords.__get__(u)
        u.center = cls.center.__get__(u)

        return u

    def unwrapped_coords(self, compound, reference):
        """Returns coordinates which correspond to the unwrapped system.

        Parameters
        ----------
        compound : {'group', 'segments', 'residues', 'molecules', 'fragments'}
            Which type of component is unwrapped.
        reference : {'com', 'cog', None}
            The reference point of each compound that is shifted into the
            primary unit cell.

        Note
        ----
        This function assumes that all atom masses are equal. Therefore, the
        returned coordinates for ``reference='com'`` and ``reference='cog'`` are
        identical.
        """
        if reference is not None:
            ref = reference.lower()
            if ref not in ['com', 'cog']:
                raise ValueError("Unknown unwrap reference: {}"
                                 "".format(reference))
        comp = compound.lower()
        if comp not in ['group', 'segments', 'residues', 'molecules',
                        'fragments']:
            raise ValueError("Unknown unwrap compound: {}".format(compound))

        # get relative positions:
        relpos = self._relpos.copy()
        # type B
        # molecule 5, atom 2 & molecule 6, atom 1 & 2
        relpos[8, :] = [0.05, 0.2, -0.05]
        relpos[10, :] = [-0.2, -0.9, 0.95]
        relpos[11, :] = [-0.1, -0.9, 0.95]
        # type C
        # molecule 8, atoms 2 & 3
        relpos[17, :] = [0.4, 1.05, 0.95]
        relpos[18, :] = [0.4, 1.05, 1.05]
        # type D
        # molecule 11, residue 1
        relpos[35:39, :] = np.array([[1.06, 0.75, 0.7],
                                     [1.16, 0.75, 0.7],
                                     [1.26, 0.75, 0.7],
                                     [1.36, 0.75, 0.7]], dtype=np.float32)

        # apply image shifts if necessary:
        if reference is None:
            if comp == 'residues':
                #second residue of molecule 11
                relpos[35:39, :] = np.array([[0.06, 0.75, 0.7],
                                             [0.16, 0.75, 0.7],
                                             [0.26, 0.75, 0.7],
                                             [0.36, 0.75, 0.7]],
                                            dtype=np.float32)
        else:
            # molecule 2 & 3
            relpos[1:3, :] = np.array([[0.4, 0.5, 0.5],
                                       [0.1, 0.5, 0.5]], dtype=np.float32)
            # molecule 6
            relpos[9:12, :] = np.array([[0.8, 0.1, 1.05],
                                        [0.8, 0.1, 0.95],
                                        [0.9, 0.1, 0.95]], dtype=np.float32)
            #molecule 8
            relpos[15:19, :] = np.array([[0.4, -0.05, 0.05],
                                         [0.4, -0.05, -0.05],
                                         [0.4, 0.05, -0.05],
                                         [0.4, 0.05, 0.05]], dtype=np.float32)
            if comp == 'residues':
                #molecule 11, residue 1 & molecule 12
                relpos[35:47, :] = np.array([[0.06, 0.75, 0.7],
                                             [0.16, 0.75, 0.7],
                                             [0.26, 0.75, 0.7],
                                             [0.36, 0.75, 0.7],
                                             [1.14, 0.65, 0.6],
                                             [1.04, 0.65, 0.6],
                                             [0.94, 0.65, 0.6],
                                             [0.84, 0.65, 0.6],
                                             [0.74, 0.65,  0.6],
                                             [0.64, 0.65,  0.6],
                                             [0.54, 0.65,  0.6],
                                             [0.44, 0.65,  0.6]],
                                            dtype=np.float32)
            else:
                #molecule 11 & 12
                relpos[31:47, :] = np.array([[-0.34, 0.75, 0.7],
                                             [-0.24, 0.75, 0.7],
                                             [-0.14, 0.75, 0.7],
                                             [-0.04, 0.75, 0.7],
                                             [0.06, 0.75, 0.7],
                                             [0.16, 0.75, 0.7],
                                             [0.26, 0.75, 0.7],
                                             [0.36, 0.75, 0.7],
                                             [1.14, 0.65, 0.6],
                                             [1.04, 0.65, 0.6],
                                             [0.94, 0.65, 0.6],
                                             [0.84, 0.65, 0.6],
                                             [0.74, 0.65, 0.6],
                                             [0.64, 0.65, 0.6],
                                             [0.54, 0.65, 0.6],
                                             [0.44, 0.65, 0.6]],
                                            dtype=np.float32)

        # apply y- and z-dependent shift of x and y coords for triclinic boxes:
        if self._is_triclinic:
            # x-coord shift depends on y- and z-coords
            relpos[:, 0] += self._tfac * relpos[:, 1:].sum(axis=1)
            # y-coord shift depends on z-coords only
            relpos[:, 1] += self._tfac * relpos[:, 2]

        # scale relative to absolute positions:
        a = self._box_edge
        positions = relpos * np.array([a, a, a])

        return positions.astype(np.float32)

    def wrapped_coords(self, compound, center):
        """Returns coordinates which correspond to the wrapped system.

        Parameters
        ----------
        compound : {'atoms', 'group', 'segments', 'residues', 'molecules', \
                    'fragments'}
            Which type of component is unwrapped. Note that for ``'group'``,
            the result will only be correct *if the group is the entire system*.
        center : {'com', 'cog'}
            The reference point of each compound that is shifted into the
            primary unit cell.

        Note
        ----
        This function assumes that all atom masses are equal. Therefore, the
        returned coordinates for ``center='com'`` and ``center='cog'`` are
        identical.
        """
        ctr = center.lower()
        if ctr not in ['com', 'cog']:
            raise ValueError("Unknown unwrap reference: {}".format(center))
        comp = compound.lower()
        if comp not in ['atoms', 'group', 'segments', 'residues',
                            'molecules', 'fragments']:
            raise ValueError("Unknown unwrap compound: {}".format(compound))

        # wrapped relative positions:
        relpos = self._relpos.copy()

        # apply required box shifts:
        if comp == 'atoms':
            # type A
            # type A
            # molecule 2: negative x-shift
            relpos[1, 0] -= 1.0
            # molecule 2: negative double x-shift
            relpos[2, 0] -= 2.0
            # type B
            # molecule 5, atom 0: positive x-shift
            relpos[6, 0] += 1.0
            # molecule 6, atom 0: positive x- and y-shift and negative z-shift
            relpos[9, :] += [1.0, 1.0, -1.0]
            # molecule 6, atom 2 & 3: positive x- and z-shift
            relpos[10:12, :] += [1.0, 0.0, 1.0]
            # molecule 7, atom 2: negative x-shift
            relpos[14, 0] -= 1.0
            # type C
            # molecule 8, atoms 0 & 3: negative z-shift
            relpos[15, 2] -= 1.0
            relpos[18, 2] -= 1.0
            # type D
            # molecule 12, atoms 0 & 1: negative x-shift
            relpos[39:41, 0] -= 1.0
            # molecule 12: positive z-shift
            relpos[39:47, 2] += 1.0
        elif comp == 'group':
            # com or cog of entire system is within box, so no shift
            pass
        elif comp == 'segments':
            # type A
            # molecules 1-3: negative x-shift
            relpos[0:3, 0] -= 1.0
            # type D
            # molecule 12: positive z-shift
            relpos[39:47, 2] += 1.0
        else:  # comp is residues, molecules, or fragments
            # type A
            # molecule 2: negative x-shift
            relpos[1, 0] -= 1.0
            # molecule 2: negative double x-shift
            relpos[2, 0] -= 2.0
            #type B
            # molecule 6: positive x- and y-shift
            relpos[9:12, :2] += 1.0
            #type C
            # molecule 8: negative z-shift
            relpos[15:19, 2] -= 1.0
            #type D
            # molecule 12: positive z-shift
            relpos[39:47, 2] += 1.0

        # apply y- and z-dependent shift of x and y coords for triclinic boxes:
        if self._is_triclinic:
            # x-coord shift depends on y- and z-coords
            relpos[:, 0] += self._tfac * relpos[:, 1:].sum(axis=1)
            # y-coord shift depends on z-coords only
            relpos[:, 1] += self._tfac * relpos[:, 2]

        # scale relative to absolute positions:
        a = self._box_edge
        positions = relpos * np.array([a, a, a])

        return positions.astype(np.float32)

    def center(self, compound):
        """Returns centers which correspond to the unwrapped system.

        Parameters
        ----------
        compound : {'atoms', 'group', 'segments', 'residues', 'molecules', \
                    'fragments'}
            Which type of component is unwrapped. Note that for ``'group'``,
            the result will only be correct *if the group is the entire system*.

        Note
        ----
        This function assumes that all atom masses are equal. Therefore, the
        returned coordinates for ``center='com'`` and ``center='cog'`` are
        identical.
        """

        relpos = self.unwrapped_coords(compound, reference=None)

        comp = compound.lower()
        if comp not in ['group', 'segments', 'residues', 'molecules',
                        'fragments']:
            raise ValueError("Unknown unwrap compound: {}".format(compound))

        pos = 0

        if compound=="residues":
            center_pos = np.zeros((15, 3), dtype=np.float32)
        else:
            center_pos = np.zeros((12, 3), dtype=np.float32)

        for base in range(3):
            loc_center = relpos[base, :]
            center_pos[pos,:] = loc_center
            pos+=1

        for base in range(3, 15, 3):
            loc_center = np.mean(relpos[base:base + 3, :], axis=0)
            center_pos[pos,:] = loc_center
            pos+=1

        if compound=="residues":
            for base in range(15, 47, 4):
                loc_center = np.mean(relpos[base:base + 4, :], axis=0)
                center_pos[pos,:] = loc_center
                pos+=1
        else:
            for base in range(15, 23, 4):
                loc_center = np.mean(relpos[base:base + 4, :], axis=0)
                center_pos[pos,:] = loc_center
                pos+=1
            for base in range(23, 47, 8):
                loc_center = np.mean(relpos[base:base + 8, :], axis=0)
                center_pos[pos,:] = loc_center
                pos+=1

        if compound == "group":
            center_pos = center_pos[11]
        elif compound == "segments":
            center_pos = center_pos[9:]

        return center_pos

