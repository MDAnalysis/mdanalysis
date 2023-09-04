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
from numpy.testing import assert_almost_equal, assert_array_equal
import pytest

import MDAnalysis as mda
from MDAnalysis import NoDataError
from MDAnalysisTests.core.util import UnWrapUniverse, assert_in_box
from MDAnalysisTests.datafiles import TRZ_psf, TRZ

class TestWrap(object):
    """Tests the functionality of *Group.wrap() using the UnWrapUniverse,
    which is specifically designed for wrapping and unwrapping tests.
    """

    precision = 5

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('atoms', 'group', 'segments',
                                          'residues', 'molecules', 'fragments'))
    @pytest.mark.parametrize('center', ('com', 'cog'))
    @pytest.mark.parametrize('is_triclinic', (False, True))
    def test_wrap_pass(self, level, compound, center, is_triclinic):
        # get a pristine test universe:
        u = UnWrapUniverse(is_triclinic=is_triclinic)
        group = u.atoms
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        # store original positions:
        orig_pos = group.atoms.positions
        # get the expected result:
        ref_wrapped_pos = u.wrapped_coords(compound, center)
        # first, do the wrapping out-of-place:
        wrapped_pos = group.wrap(compound=compound, center=center,
                                 inplace=False)
        # check for correct result:
        assert_almost_equal(wrapped_pos, ref_wrapped_pos,
                            decimal=self.precision)
        # make sure atom positions are unchanged:
        assert_array_equal(group.atoms.positions, orig_pos)
        # now, do the wrapping inplace:
        wrapped_pos2 = group.wrap(compound=compound, center=center,
                                  inplace=True)
        # check that result is the same as for out-of-place computation:
        assert_array_equal(wrapped_pos, wrapped_pos2)
        # check that wrapped positions are applied:
        assert_array_equal(group.atoms.positions, wrapped_pos)
        # check that nobody messed with the reference positions,
        # centers of compounds must lie within the primary unit cell:
        if compound == 'atoms':
            assert_in_box(group.atoms.positions, group.dimensions)
        elif center == 'com':
            compos = group.atoms.center_of_mass(wrap=False, compound=compound)
            assert_in_box(compos, group.dimensions)
        else:
            cogpos = group.atoms.center_of_geometry(wrap=False,
                                                    compound=compound)
            assert_in_box(cogpos, group.dimensions)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('atoms', 'group', 'segments',
                                          'residues', 'molecules', 'fragments'))
    @pytest.mark.parametrize('center', ('com', 'cog'))
    @pytest.mark.parametrize('is_triclinic', (False, True))
    def test_unwrap_wrap_cycle(self, level, compound, center, is_triclinic):
        # get a pristine test universe:
        u = UnWrapUniverse(is_triclinic=is_triclinic)
        # set wrapped reference coordinates:
        u.atoms.positions = u.wrapped_coords('atoms', 'com')
        group = u.atoms
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        # unwrap:
        group.unwrap()
        # store original unwrapped positions:
        orig_unwrapped_pos = group.atoms.positions
        # wrap:
        group.wrap(compound=compound, center=center, inplace=True)
        # unwrap again:
        group.unwrap()
        # make sure unwrapped atom positions are as before:
        assert_almost_equal(group.atoms.positions, orig_unwrapped_pos,
                            decimal=self.precision)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('atoms', 'group', 'segments',
                                          'residues', 'molecules', 'fragments'))
    @pytest.mark.parametrize('center', ('com', 'cog'))
    @pytest.mark.parametrize('is_triclinic', (False, True))
    def test_wrap_partial_compound(self, level, compound, center, is_triclinic):
        # get a pristine test universe:
        u = UnWrapUniverse(is_triclinic=is_triclinic)
        group = u.atoms
        ref_wrapped_pos = u.wrapped_coords(compound, center)
        # select topology level with one missing item and get expected result:
        if level == 'atoms':
            # first atom of molecule 12 missing
            missing = group[[-8]]
            group = group[:-8] + group[-7:]
            ref_wrapped_pos = np.concatenate([ref_wrapped_pos[:-8],
                                              ref_wrapped_pos[-7:]])
        elif level == 'residues':
            group = group.residues
            # first residue of molecule 12 missing
            missing = group[-2]
            group = group[:-2] + group[[-1]]
            ref_wrapped_pos = np.concatenate([ref_wrapped_pos[:-8],
                                              ref_wrapped_pos[-4:]])
        elif level == 'segments':
            group = group.segments
            # molecule 12 missing
            missing = group[-1]
            group = group[:-1]
            ref_wrapped_pos = ref_wrapped_pos[:-8]
        # store original positions of the missing item:
        orig_pos = missing.atoms.positions
        # first, do the wrapping out-of-place:
        group.wrap(compound=compound, center=center, inplace=True)
        # check for correct result:
        assert_almost_equal(group.atoms.positions, ref_wrapped_pos,
                            decimal=self.precision)
        # make sure the positions of the missing item are unchanged:
        assert_array_equal(missing.atoms.positions, orig_pos)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('atoms', 'group', 'segments',
                                          'residues', 'molecules', 'fragments'))
    @pytest.mark.parametrize('center', ('com', 'cog'))
    @pytest.mark.parametrize('is_triclinic', (False, True))
    def test_wrap_empty_group(self, level, compound, center, is_triclinic):
        # get a pristine test universe:
        u = UnWrapUniverse(is_triclinic=is_triclinic)
        group = u.atoms[[]]
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        group.wrap(compound=compound, center=center, inplace=True)
        # check for correct (empty) result:
        assert_array_equal(group.atoms.positions,
                           np.empty((0, 3), dtype=np.float32))

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('atoms', 'group', 'segments',
                                          'residues', 'molecules', 'fragments'))
    @pytest.mark.parametrize('center', ('com', 'cog'))
    @pytest.mark.parametrize('is_triclinic', (False, True))
    def test_wrap_duplicates(self, level, compound, center, is_triclinic):
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
        if compound == 'group':
            # reference positions of UnWrapUniverse are known to be incorrect
            # for compound='group' if the group is not the entire system, so we
            # have to correct for that:
            ref_wrapped_pos = u.wrapped_coords('segments', center)[39:47]
        else:
            ref_wrapped_pos = u.wrapped_coords(compound, center)[39:47]
        ref_wrapped_pos = np.vstack((ref_wrapped_pos, ref_wrapped_pos))
        # wrap:
        group.wrap(compound=compound, center=center, inplace=True)
        # check for correct result:
        assert_almost_equal(group.atoms.positions, ref_wrapped_pos,
                            decimal=self.precision)
        # check that the rest of the atoms are kept unmodified:
        assert_array_equal(rest.positions, orig_rest_pos)

    @pytest.mark.parametrize('compound', ('atoms', 'group', 'segments',
                                          'residues', 'molecules', 'fragments'))
    @pytest.mark.parametrize('is_triclinic', (False, True))
    def test_wrap_com_cog_difference(self, compound, is_triclinic):
        # get a pristine test universe:
        u = UnWrapUniverse(is_triclinic=is_triclinic)
        # select molecule 5:
        group = u.atoms[6:9]
        # make first atom of molecule 5 much more heavy than the other two.
        # That way, the whole molecule's center of geometry will still lie
        # inside the first unit cell but its center of mass will lie outside
        # the first unit cell in negative x-direction.
        group.masses = [100.0, 1.0, 1.0]
        # wrap with center='cog':
        wrapped_pos_cog = group.wrap(compound=compound, center='cog',
                                     inplace=False)
        # get expected result:
        ref_wrapped_pos = u.wrapped_coords(compound, 'cog')[6:9]
        # check for correctness:
        assert_almost_equal(wrapped_pos_cog, ref_wrapped_pos,
                            decimal=self.precision)
        # wrap with center='com':
        wrapped_pos_com = group.wrap(compound=compound, center='com',
                                     inplace=False)
        # assert that the com result is shifted with respect to the cog result
        # by one box length in the x-direction:
        shift = np.array([10.0, 0.0, 0.0], dtype=np.float32)
        if compound == 'atoms':
            # center argument must be ignored for compound='atoms':
            shift[0] = 0.0
        assert_almost_equal(wrapped_pos_cog, wrapped_pos_com - shift,
                            decimal=self.precision)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('atoms', 'group', 'segments',
                                          'residues', 'molecules', 'fragments'))
    def test_wrap_zero_mass_exception_safety(self, level, compound):
        # get a pristine test universe:
        u = UnWrapUniverse()
        # set masses of molecule 12 to zero:
        u.atoms[39:47].masses = 0.0
        # select group appropriate for compound:
        if compound == 'group':
            group = u.atoms[39:47] # molecule 12
        else:
            group = u.atoms
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        # store original positions:
        orig_pos = group.atoms.positions
        if compound == 'atoms':
            # wrap() must not care about masses if compound == 'atoms':
            group.wrap(compound=compound, center='com', inplace=True)
            ref_wrapped_pos = u.wrapped_coords(compound, 'com')
            assert_almost_equal(group.atoms.positions, ref_wrapped_pos,
                                decimal=self.precision)
        else:
            # try to wrap:
            with pytest.raises(ValueError):
                group.wrap(compound=compound, center='com', inplace=True)
            # make sure atom positions are unchanged:
            assert_array_equal(group.atoms.positions, orig_pos)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('atoms', 'group', 'segments',
                                          'residues', 'molecules', 'fragments'))
    def test_wrap_wrong_center_exception_safety(self, level, compound):
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
        if compound == 'atoms':
            # wrap() must ignore the center argument if compound == 'atoms':
            group.wrap(compound=compound, center='com', inplace=True)
            ref_wrapped_pos = u.wrapped_coords(compound, 'com')
            assert_almost_equal(group.atoms.positions, ref_wrapped_pos,
                                decimal=self.precision)
        else:
            # try to wrap:
            with pytest.raises(ValueError):
                group.wrap(compound=compound, center='wrong', inplace=True)
            # make sure atom positions are unchanged:
            assert_array_equal(group.atoms.positions, orig_pos)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('center', ('com', 'cog'))
    def test_wrap_wrong_compound_exception_safety(self, level, center):
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
        # try to wrap:
        with pytest.raises(ValueError):
            group.wrap(compound='wrong', center=center, inplace=True)
        # make sure atom positions are unchanged:
        assert_array_equal(group.atoms.positions, orig_pos)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('atoms', 'group', 'segments',
                                          'residues', 'molecules', 'fragments'))
    def test_unwrap_no_masses_exception_safety(self, level, compound):
        # universe without masses:
        u = UnWrapUniverse(have_masses=False)
        group = u.atoms
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        if compound == 'atoms':
            # wrap() must not care about mass presence if compound == 'atoms':
            group.wrap(compound=compound, center='com', inplace=True)
            ref_wrapped_pos = u.wrapped_coords(compound, 'com')
            assert_almost_equal(group.atoms.positions, ref_wrapped_pos,
                                decimal=self.precision)
        else:
            # store original positions:
            orig_pos = group.atoms.positions
            # try to wrap:
            with pytest.raises(NoDataError):
                group.wrap(compound=compound, center='com', inplace=True)
            # make sure atom positions are unchanged:
            assert_array_equal(group.atoms.positions, orig_pos)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('atoms', 'group', 'segments',
                                          'residues', 'molecules', 'fragments'))
    @pytest.mark.parametrize('center', ('com', 'cog'))
    def test_wrap_no_bonds_exception_safety(self, level, compound, center):
        # universe without bonds:
        u = UnWrapUniverse(have_bonds=False)
        group = u.atoms
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        if compound == 'fragments':
            # store original positions:
            orig_pos = group.atoms.positions
            # must raise an exception for fragments
            with pytest.raises(NoDataError):
                group.wrap(compound=compound, center=center, inplace=True)
            # make sure atom positions are unchanged:
            assert_array_equal(group.atoms.positions, orig_pos)
        else:
            # must not care about bonds if compound != 'fragments'
            group.wrap(compound=compound, center=center, inplace=True)
            ref_wrapped_pos = u.wrapped_coords(compound, center)
            assert_almost_equal(group.atoms.positions, ref_wrapped_pos,
                                decimal=self.precision)

    @pytest.mark.parametrize('level', ('atoms', 'residues', 'segments'))
    @pytest.mark.parametrize('compound', ('atoms', 'group', 'segments',
                                          'residues', 'molecules', 'fragments'))
    @pytest.mark.parametrize('center', ('com', 'cog'))
    def test_wrap_no_molnums_exception_safety(self, level, compound, center):
        # universe without bonds:
        u = UnWrapUniverse(have_molnums=False)
        group = u.atoms
        # select topology level:
        if level == 'residues':
            group = group.residues
        elif level == 'segments':
            group = group.segments
        if compound == 'molecules':
            # store original positions:
            orig_pos = group.atoms.positions
            # must raise an exception for molecules
            with pytest.raises(NoDataError):
                group.wrap(compound=compound, center=center, inplace=True)
            # make sure atom positions are unchanged:
            assert_array_equal(group.atoms.positions, orig_pos)
        else:
            # must not care about molnums if compound != 'molecules'
            group.wrap(compound=compound, center=center, inplace=True)
            ref_wrapped_pos = u.wrapped_coords(compound, center)
            assert_almost_equal(group.atoms.positions, ref_wrapped_pos,
                                decimal=self.precision)

class TestWrapTRZ(object):
    """Tests the functionality of AtomGroup.wrap() using a TRZ universe."""

    precision = 5

    @pytest.fixture()
    def u(self):
        return mda.Universe(TRZ_psf, TRZ)

    @pytest.mark.parametrize('box', (np.array([1, 1]), np.zeros(6)))
    def test_wrap_box_fail(self, u, box):
        ag = u.atoms
        with pytest.raises(ValueError):
            ag.wrap(box=box)

    def test_wrap_small_box(self, u):
        ag = u.atoms[:100]
        box = u.dimensions
        box[:3] *= 0.5
        ag.wrap(box=box)
        assert_in_box(ag.positions, box)

    def test_wrap_large_box(self, u):
        ag = u.atoms[:100]
        box = u.dimensions
        box[:3] *= 2.0
        ag.wrap(box=box)
        assert_in_box(ag.positions, box)

    def test_wrap_atoms(self, u):
        ag = u.atoms[100:200]
        rescom = ag.wrap(compound='atoms', center='com', inplace=False)
        assert_in_box(rescom, u.dimensions)
        rescog = ag.wrap(compound='atoms', center='cog', inplace=False)
        assert_almost_equal(rescom, rescog, decimal=self.precision)

    @pytest.mark.parametrize('center', ('com', 'cog'))
    def test_wrap_group(self, u, center):
        ag = u.atoms[:100]
        ag.wrap(compound='group', center=center)
        if center == 'com':
            ctrpos = ag.center_of_mass(wrap=False, compound='group')
        elif center == 'cog':
            ctrpos = ag.center_of_geometry(wrap=False, compound='group')
        assert_in_box(ctrpos, u.dimensions)

    @pytest.mark.parametrize('center', ('com', 'cog'))
    def test_wrap_residues(self, u, center):
        ag = u.atoms[300:400].residues
        ag.wrap(compound='residues', center=center)
        if center == 'com':
            ctrpos = ag.center_of_mass(wrap=False, compound='residues')
        elif center == 'cog':
            ctrpos = ag.center_of_geometry(wrap=False, compound='residues')
        assert_in_box(ctrpos, u.dimensions)

    @pytest.mark.parametrize('center', ('com', 'cog'))
    def test_wrap_segments(self, u, center):
        ag = u.atoms[1000:1200]
        ag.wrap(compound='segments', center=center)
        if center == 'com':
            ctrpos = ag.center_of_mass(wrap=False, compound='segments')
        elif center == 'cog':
            ctrpos = ag.center_of_geometry(wrap=False, compound='segments')
        assert_in_box(ctrpos, u.dimensions)

    @pytest.mark.parametrize('center', ('com', 'cog'))
    def test_wrap_fragments(self, u, center):
        ag = u.atoms[:250]
        ag.wrap(compound='fragments', center=center)
        if center == 'com':
            ctrpos = ag.center_of_mass(wrap=False, compound='fragments')
        elif center == 'cog':
            ctrpos = ag.center_of_geometry(wrap=False, compound='fragments')
        assert_in_box(ctrpos, u.dimensions)
