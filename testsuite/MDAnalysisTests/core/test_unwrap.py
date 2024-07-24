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
import numpy as np
from numpy.testing import (assert_raises, assert_almost_equal,
                           assert_array_equal)
import pytest

import MDAnalysis as mda
from MDAnalysis import NoDataError
from MDAnalysisTests.core.util import UnWrapUniverse
from MDAnalysis.tests.datafiles import CONECT


class TestUnwrap(object):
    """Tests the functionality of *Group.unwrap() using the UnWrapUniverse,
    which is specifically designed for wrapping and unwrapping tests.
    """
    precision = 5

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('fragments', 'molecules', 'residues',
                                          'group', 'segments'))
    @pytest.mark.parametrize('reference', ('com', 'cog', None))
    @pytest.mark.parametrize('is_triclinic', (False, True))
    def test_unwrap_pass(self, level, compound, reference, is_triclinic):
        # get a pristine test universe:
        u = UnWrapUniverse(is_triclinic=is_triclinic)
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
        ref_unwrapped_pos = u.unwrapped_coords(compound, reference)
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
        u = UnWrapUniverse(is_triclinic=is_triclinic)
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
        u = UnWrapUniverse(is_triclinic=is_triclinic)
        # select group with one atom missing
        group = u.atoms[39:46] # molecule 12 without its last atom
        # store original position of last atom of molecule 12:
        orig_pos = u.atoms[46].position
        # get the expected result:
        ref_unwrapped_pos = u.unwrapped_coords(compound, reference)[39:46]
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
        u = UnWrapUniverse(is_triclinic=is_triclinic)
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
        u = UnWrapUniverse(is_triclinic=is_triclinic)
        # select molecule 12:
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
        ref_unwrapped_pos = u.unwrapped_coords(compound, reference)[39:47]
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
        u = UnWrapUniverse(is_triclinic=is_triclinic)
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
        ref_unwrapped_pos = u.unwrapped_coords(compound, 'cog')[6:9]
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
        u = UnWrapUniverse()
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
        u = UnWrapUniverse()
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
        u = UnWrapUniverse()
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
        u = UnWrapUniverse(have_masses=False)
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
        u = UnWrapUniverse(have_bonds=False)
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
        error_message = (
            f"{group.__class__.__name__}.unwrap() not available; this AtomGroup lacks defined bonds. "
            "To resolve this, you can either:\n"
            "1. Guess the bonds at universe creation using `guess_bonds = True`, or\n"
            "2. Create a universe using a topology format where bonds are pre-defined."
        )
        with pytest.raises(NoDataError, match=error_message):
            group.unwrap(compound=compound, reference=reference, inplace=True)
        # make sure atom positions are unchanged:
        assert_array_equal(group.atoms.positions, orig_pos)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('reference', ('com', 'cog', None))
    def test_unwrap_no_molnums_exception_safety(self, level, reference):
        # universe without molnums:
        u = UnWrapUniverse(have_molnums=False)
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


def test_uncontiguous():
    """Real-life case of fragment sparsity that triggers Issue 3352
    """
    precision = 5
    displacement_vec = [14.7, 0., 0.]
    u = mda.Universe(CONECT)
    # This is one of the few residues that has bonds
    ag = u.residues[66].atoms
    ref_pos = ag.positions
    # Let's break it by placing it over the box boundary and re-packing
    u.atoms.positions -= displacement_vec
    u.atoms.pack_into_box()
    # Let's make sure we really broke the fragment
    assert_raises(AssertionError, assert_almost_equal,
                  ref_pos, ag.positions+displacement_vec,
                  decimal=precision)
    # Ok, let's make it whole again and check that we're good
    u.atoms.unwrap()
    assert_almost_equal(ref_pos, ag.positions+displacement_vec,
                        decimal=precision)
