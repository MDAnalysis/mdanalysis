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
from numpy.testing import assert_almost_equal, assert_array_equal
import pytest

import MDAnalysis as mda
from MDAnalysis.core import topologyattrs
from MDAnalysis import NoDataError
from MDAnalysis.lib.mdamath import triclinic_box


def unwrap_test_universe(have_bonds=True, have_masses=True, have_molnums=True,
                         is_triclinic=False):
    """Construct a universe comprising small molecule clusters with molecules
       (as well as clusters) broken or whole accross periodic boundaries:

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
    is_triclinic : bool, optional
        If ``False``, the box will be a cube with an edge length of 10 Angstrom.
        If ``True``, the cubic box will be sheared by 1 Angstrom in the x and
        y directions.

    Returns
    -------
    u : MDAnalysis.Universe
        The test universe as specified above.
    """
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
    relpos[15:23, :] = np.array([[0.4, 0.96, 1.06],
                                 [0.4, 0.96, 0.96],
                                 [0.4, 0.06, 0.96],
                                 [0.4, 0.06, 1.06],
                                 [0.6, 0.06, 0.26],
                                 [0.6, 0.06, 0.16],
                                 [0.6, 0.16, 0.16],
                                 [0.6, 0.16, 0.26]], dtype=np.float32)
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
    # apply y- and z-dependent shift of x and y coords for triclinic boxes:
    if is_triclinic:
        # x-coord shift depends on y- and z-coords
        relpos[:, 0] += tfac * relpos[:, 1:].sum(axis=1)
        # y-coord shift depends on z-coords only
        relpos[:, 1] += tfac * relpos[:, 2]
    # scale relative to absolute positions:
    ts.positions = relpos * np.array([a, a, a], dtype=np.float32)

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

    return u


def unwrapped_coords(compound, reference, is_triclinic):
    """Returns coordinates which correspond to the system returned
    by :func:`unwrap_test_universe` but are unwrapped.

    Parameters
    ----------
    compound : {'group', 'segments', 'residues', 'molecules', 'fragments'}
        Which type of component is unwrapped.
    reference : {'com', 'cog', None}
        The reference point of each compound that is shifted into the primary
        unit cell.
    is_triclinic : bool
        Is the box triclinic?
    
    Note
    ----
    This function assumes that all atom masses are equal. Therefore, the
    returned coordinates for ``reference='com'`` and ``reference='cog'`` are
    identical.
    """
    if reference not in ['com', 'cog', None]:
        raise ValueError("Unknown unwrap reference: {}".format(reference))
    if compound not in ['group', 'segments', 'residues', 'molecules',
                        'fragments']:
        raise ValueError("Unknown unwrap compound: {}".format(compound))

    n_atoms = 47

    # box:
    a = 10.0 # edge length
    tfac = 0.1  # factor for box vector shift of triclinic boxes (~84°)
    if is_triclinic:
        box = triclinic_box([a, 0.0, 0.0],
                            [a * tfac, a, 0.0],
                            [a * tfac, a * tfac, a])
    else:
        box = np.array([a, a, a, 90.0, 90.0, 90.0], dtype=np.float32)

    # unwrapped positions:
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
                                [0.05, 0.2, -0.05],
                                [-0.2, -0.9, 1.05],
                                [-0.2, -0.9, 0.95],
                                [-0.1, -0.9, 0.95],
                                [0.95, 0.2, 0.25],
                                [0.95, 0.2, 0.15],
                                [1.05, 0.2, 0.15]], dtype=np.float32)
    # type C
    relpos[15:23, :] = np.array([[0.4, 0.96, 1.06],
                                 [0.4, 0.96, 0.96],
                                 [0.4, 1.06, 0.96],
                                 [0.4, 1.06, 1.06],
                                 [0.6, 0.06, 0.26],
                                 [0.6, 0.06, 0.16],
                                 [0.6, 0.16, 0.16],
                                 [0.6, 0.16, 0.26]], dtype=np.float32)
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
                                 [1.06, 0.75, 0.7],
                                 [1.16, 0.75, 0.7],
                                 [1.26, 0.75, 0.7],
                                 [1.36, 0.75, 0.7],
                                 [1.14, 0.65, -0.4],
                                 [1.04, 0.65, -0.4],
                                 [0.94, 0.65, -0.4],
                                 [0.84, 0.65, -0.4],
                                 [0.74, 0.65, -0.4],
                                 [0.64, 0.65, -0.4],
                                 [0.54, 0.65, -0.4],
                                 [0.44, 0.65, -0.4]], dtype=np.float32)
    # apply image shifts if necessary:
    if reference is None:
        if compound == 'residues':
            #second residue of molecule 11
            relpos[35:39, :] = np.array([[0.06, 0.75, 0.7],
                                         [0.16, 0.75, 0.7],
                                         [0.26, 0.75, 0.7],
                                         [0.36, 0.75, 0.7]], dtype=np.float32)
    else:
        # molecule 2 & 3
        relpos[1:3, :] = np.array([[0.4, 0.5, 0.5],
                                   [0.1, 0.5, 0.5]], dtype=np.float32)
        # molecule 6
        relpos[9:12, :] = np.array([[0.8, 0.1, 1.05],
                                    [0.8, 0.1, 0.95],
                                    [0.9, 0.1, 0.95]], dtype=np.float32)
        #molecule 8
        relpos[15:19, :] = np.array([[0.4, -0.04, 0.06],
                                     [0.4, -0.04, -0.04],
                                     [0.4, 0.06, -0.04],
                                     [0.4, 0.06, 0.06]], dtype=np.float32)
        if compound == 'residues':
            #molecule 11 & 12
            relpos[31:47, :] = np.array([[0.66, 0.75, 0.7],
                                         [0.76, 0.75, 0.7],
                                         [0.86, 0.75, 0.7],
                                         [0.96, 0.75, 0.7],
                                         [0.06, 0.75, 0.7],
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
                                         [0.44, 0.65,  0.6]], dtype=np.float32)
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
                                         [0.44, 0.65, 0.6]], dtype=np.float32)

    # apply y- and z-dependent shift of x and y coords for triclinic boxes:
    if is_triclinic:
        # x-coord shift depends on y- and z-coords
        relpos[:, 0] += tfac * relpos[:, 1:].sum(axis=1)
        # y-coord shift depends on z-coords only
        relpos[:, 1] += tfac * relpos[:, 2]

    # scale relative to absolute positions:
    positions = relpos * np.array([a, a, a], dtype=np.float32)

    return positions


class TestUnwrap(object):
    """Tests the functionality of *Group.unwrap() using the test universe
    returned by unwrap_test_universe()
    """
    precision = 5

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('fragments', 'molecules', 'residues',
                                          'group', 'segments'))
    @pytest.mark.parametrize('reference', ('com', 'cog', None))
    @pytest.mark.parametrize('is_triclinic', (False, True))
    def test_unwrap_pass(self, level, compound, reference, is_triclinic):
        # get a pristine test universe:
        u = unwrap_test_universe(is_triclinic=is_triclinic)
        # select group appropriate for compound:
        if compound == 'group':
            group = u.atoms[39:47] # molecule 12
        elif compound == 'segments':
            group = u.atoms[23:47] # molecules 10, 11, 12
        else:
            group = u.atoms
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        # store original positions:
        orig_pos = group.atoms.positions
        # get the expected result:
        ref_unwrapped_pos = unwrapped_coords(compound, reference,
                                             is_triclinic=is_triclinic)
        if compound == 'group':
            ref_unwrapped_pos = ref_unwrapped_pos[39:47] # molecule 12
        elif compound == 'segments':
            ref_unwrapped_pos = ref_unwrapped_pos[23:47] # molecules 10, 11, 12
        # first, do the unwrapping out-of-place:
        unwrapped_pos = group.unwrap(compound=compound, reference=reference,
                                     inplace=False)
        # check for correct result:
        assert_almost_equal(unwrapped_pos, ref_unwrapped_pos,
                            decimal=self.precision)
        # make sure atom positions are unchanged:
        assert_array_equal(group.atoms.positions, orig_pos)
        # now, do the unwrapping inplace:
        unwrapped_pos2 = group.unwrap(compound=compound, reference=reference,
                                     inplace=True)
        # check that result is the same as for out-of-place computation:
        assert_array_equal(unwrapped_pos, unwrapped_pos2)
        # check that unwrapped positions are applied:
        assert_array_equal(group.atoms.positions, unwrapped_pos)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('fragments', 'molecules', 'residues',
                                          'group', 'segments'))
    @pytest.mark.parametrize('reference', ('com', 'cog', None))
    @pytest.mark.parametrize('is_triclinic', (False, True))
    def test_wrap_unwrap_cycle(self, level, compound, reference, is_triclinic):
        # get a pristine test universe:
        u = unwrap_test_universe(is_triclinic=is_triclinic)
        # select group appropriate for compound:
        if compound == 'group':
            group = u.atoms[39:47] # molecule 12
        elif compound == 'segments':
            group = u.atoms[23:47] # molecules 10, 11, 12
        else:
            group = u.atoms
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        # wrap:
        group.wrap()
        # store original wrapped positions:
        orig_wrapped_pos = group.atoms.positions
        # unwrap:
        group.unwrap(compound=compound, reference=reference, inplace=True)
        # wrap again:
        group.wrap()
        # make sure wrapped atom positions are as before:
        assert_almost_equal(group.atoms.positions, orig_wrapped_pos,
                            decimal=self.precision)

    @pytest.mark.parametrize('compound', ('fragments', 'molecules', 'residues',
                                          'group', 'segments'))
    @pytest.mark.parametrize('reference', ('com', 'cog', None))
    @pytest.mark.parametrize('is_triclinic', (False, True))
    def test_unwrap_partial_frags(self, compound, reference, is_triclinic):
        # get a pristine test universe:
        u = unwrap_test_universe(is_triclinic=is_triclinic)
        # select group with one atom missing
        group = u.atoms[39:46] # molecule 12 without its last atom
        # store original position of last atom of molecule 12:
        orig_pos = u.atoms[46].position
        # get the expected result:
        ref_unwrapped_pos = unwrapped_coords(compound, reference,
                                             is_triclinic=is_triclinic)[39:46]
        # first, do the unwrapping out-of-place:
        group.unwrap(compound=compound, reference=reference, inplace=True)
        # check for correct result:
        assert_almost_equal(group.positions, ref_unwrapped_pos,
                            decimal=self.precision)
        # make sure the position of molecule 12's last atom is unchanged:
        assert_array_equal(u.atoms[46].position, orig_pos)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('fragments', 'molecules', 'residues',
                                          'group', 'segments'))
    @pytest.mark.parametrize('reference', ('com', 'cog', None))
    @pytest.mark.parametrize('is_triclinic', (False, True))
    def test_unwrap_empty_group(self, level, compound, reference, is_triclinic):
        # get a pristine test universe:
        u = unwrap_test_universe(is_triclinic=is_triclinic)
        if level == 'atoms':
            group = mda.AtomGroup([], u)
        elif level == 'residues':
            group = mda.ResidueGroup([], u)
        elif level == 'segments':
            group = mda.SegmentGroup([], u)
        group.unwrap(compound=compound, reference=reference, inplace=True)
        # check for correct (empty) result:
        assert_array_equal(group.atoms.positions,
                           np.empty((0, 3), dtype=np.float32))

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('fragments', 'molecules', 'residues',
                                          'group', 'segments'))
    @pytest.mark.parametrize('reference', ('com', 'cog', None))
    @pytest.mark.parametrize('is_triclinic', (False, True))
    def test_unwrap_duplicates(self, level, compound, reference, is_triclinic):
        # get a pristine test universe:
        u = unwrap_test_universe(is_triclinic=is_triclinic)
        # select molecule 6:
        group = u.atoms[39:47]
        # select the rest of the universe's atoms:
        rest = u.atoms[:39]
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        # duplicate the group:
        group += group
        # store original positions of the rest:
        orig_rest_pos = rest.positions
        # get the expected result with duplicates:
        ref_unwrapped_pos = unwrapped_coords(compound, reference,
                                             is_triclinic=is_triclinic)[39:47]
        ref_unwrapped_pos = np.vstack((ref_unwrapped_pos, ref_unwrapped_pos))
        # unwrap:
        group.unwrap(compound=compound, reference=reference, inplace=True)
        # check for correct result:
        assert_almost_equal(group.atoms.positions, ref_unwrapped_pos,
                            decimal=self.precision)
        # check that the rest of the atoms are kept unmodified:
        assert_array_equal(rest.positions, orig_rest_pos)

    @pytest.mark.parametrize('compound', ('fragments', 'molecules', 'residues',
                                          'group', 'segments'))
    @pytest.mark.parametrize('is_triclinic', (False, True))
    def test_unwrap_com_cog_difference(self, compound, is_triclinic):
        # get a pristine test universe:
        u = unwrap_test_universe(is_triclinic=is_triclinic)
        # select molecule 5:
        group = u.atoms[6:9]
        # make first atom of molecule 5 much more heavy than the other two.
        # That way, the whole molecule's center of geometry will still lie
        # inside the first unit cell but its center of mass will lie outside
        # the first unit cell in negative x-direction.
        group.masses = [100.0, 1.0, 1.0]
        # unwrap with center of geometry as reference:
        unwrapped_pos_cog = group.unwrap(compound=compound, reference='cog',
                                         inplace=False)
        # get expected result:
        ref_unwrapped_pos = unwrapped_coords(compound, reference='cog',
                                             is_triclinic=is_triclinic)[6:9]
        # check for correctness:
        assert_almost_equal(unwrapped_pos_cog, ref_unwrapped_pos,
                            decimal=self.precision)
        # unwrap with center of mass as reference:
        unwrapped_pos_com = group.unwrap(compound=compound, reference='com',
                                         inplace=False)
        # assert that the com result is shifted with respect to the cog result
        # by one box length in the x-direction:
        shift = np.array([10.0, 0.0, 0.0], dtype=np.float32)
        assert_almost_equal(unwrapped_pos_cog, unwrapped_pos_com - shift,
                            decimal=self.precision)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('fragments', 'molecules', 'residues',
                                          'group', 'segments'))
    def test_unwrap_zero_mass_exception_safety(self, level, compound):
        # get a pristine test universe:
        u = unwrap_test_universe()
        # set masses of molecule 12 to zero:
        u.atoms[39:47].masses = 0.0
        # select group appropriate for compound:
        if compound == 'group':
            group = u.atoms[39:47] # molecule 12
        elif compound == 'segments':
            group = u.atoms[23:47] # molecules 10, 11, 12
        else:
            group = u.atoms
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        # store original positions:
        orig_pos = group.atoms.positions
        # try to unwrap:
        with pytest.raises(ValueError):
            group.unwrap(compound=compound, reference='com',
                         inplace=True)
        # make sure atom positions are unchanged:
        assert_array_equal(group.atoms.positions, orig_pos)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('fragments', 'molecules', 'residues',
                                          'group', 'segments'))
    def test_unwrap_wrong_reference_exception_safety(self, level, compound):
        # get a pristine test universe:
        u = unwrap_test_universe()
        # select group appropriate for compound:
        if compound == 'group':
            group = u.atoms[39:47] # molecule 12
        elif compound == 'segments':
            group = u.atoms[23:47] # molecules 10, 11, 12
        else:
            group = u.atoms
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        # store original positions:
        orig_pos = group.atoms.positions
        # try to unwrap:
        with pytest.raises(ValueError):
            group.unwrap(compound=compound, reference='wrong', inplace=True)
        # make sure atom positions are unchanged:
        assert_array_equal(group.atoms.positions, orig_pos)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('reference', ('com', 'cog', None))
    def test_unwrap_wrong_compound_exception_safety(self, level, reference):
        # get a pristine test universe:
        u = unwrap_test_universe()
        group = u.atoms
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        # store original positions:
        orig_pos = group.atoms.positions
        # try to unwrap:
        with pytest.raises(ValueError):
            group.unwrap(compound='wrong', reference=reference, inplace=True)
        # make sure atom positions are unchanged:
        assert_array_equal(group.atoms.positions, orig_pos)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('fragments', 'molecules', 'residues',
                                          'group', 'segments'))
    def test_unwrap_no_masses_exception_safety(self, level, compound):
        # universe without masses:
        u = unwrap_test_universe(have_masses=False)
        # select group appropriate for compound:
        if compound == 'group':
            group = u.atoms[39:47] # molecule 12
        elif compound == 'segments':
            group = u.atoms[23:47] # molecules 10, 11, 12
        else:
            group = u.atoms
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        # store original positions:
        orig_pos = group.atoms.positions
        # try to unwrap:
        with pytest.raises(NoDataError):
            group.unwrap(compound=compound, reference='com', inplace=True)
        # make sure atom positions are unchanged:
        assert_array_equal(group.atoms.positions, orig_pos)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('fragments', 'molecules', 'residues',
                                          'group', 'segments'))
    @pytest.mark.parametrize('reference', ('com', 'cog', None))
    def test_unwrap_no_bonds_exception_safety(self, level, compound, reference):
        # universe without bonds:
        u = unwrap_test_universe(have_bonds=False)
        # select group appropriate for compound:
        if compound == 'group':
            group = u.atoms[39:47] # molecule 12
        elif compound == 'segments':
            group = u.atoms[23:47] # molecules 10, 11, 12
        else:
            group = u.atoms
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        # store original positions:
        orig_pos = group.atoms.positions
        with pytest.raises(NoDataError):
            group.unwrap(compound=compound, reference=reference, inplace=True)
        # make sure atom positions are unchanged:
        assert_array_equal(group.atoms.positions, orig_pos)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('reference', ('com', 'cog', None))
    def test_unwrap_no_molnums_exception_safety(self, level, reference):
        # universe without molnums:
        u = unwrap_test_universe(have_molnums=False)
        group = u.atoms
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        # store original positions:
        orig_pos = group.atoms.positions
        with pytest.raises(NoDataError):
            group.unwrap(compound='molecules', reference=reference,
                         inplace=True)
        assert_array_equal(group.atoms.positions, orig_pos)
