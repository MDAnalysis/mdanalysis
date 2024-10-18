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
from contextlib import contextmanager

import MDAnalysis as mda
import MDAnalysis.analysis.align as align
from MDAnalysis.analysis.align import HAS_BIOPYTHON
import MDAnalysis.analysis.rms as rms
import os
import numpy as np
import pytest
from MDAnalysis import SelectionError, SelectionWarning
from MDAnalysisTests import executable_not_found
from MDAnalysisTests.datafiles import (PSF, DCD, CRD, FASTA, ALIGN_BOUND,
                                       ALIGN_UNBOUND, PDB_helix)
from numpy.testing import (
    assert_equal,
    assert_array_equal,
    assert_allclose,
)

#Function for Parametrizing conditional raising
@contextmanager
def does_not_raise():
    yield


class TestRotationMatrix:
    a = np.array([[0.1, 0.2, 0.3], [1.1, 1.1, 1.1]])
    b = np.array([[0.1, 0.1, 0.1], [1.1, 1.1, 1.1]])
    w = np.array([1.3, 2.3])

    @pytest.mark.parametrize('a, b, weights, expected', (
            (a, b, None, 0.15785647734415692),
            (a, b, w, 0.13424643502242328),
    ))
    def test_rotation_matrix_input(self, a, b, weights, expected):
        rot, rmsd = align.rotation_matrix(a, b, weights)
        assert_equal(rot, np.eye(3))
        assert rmsd == pytest.approx(expected)

    def test_list_args(self):
        a = [[0.1, 0.2, 0.3], [1.1, 1.1, 1.1]]
        b = [[0.1, 0.1, 0.1], [1.1, 1.1, 1.1]]
        w = [1.3, 2.3]
        rot, rmsd = align.rotation_matrix(a, b, w)
        assert_equal(rot, np.eye(3))
        assert rmsd == pytest.approx(0.13424643502242328)

    def test_exception(self):
        a = [[0.1, 0.2, 0.3],
             [1.1, 1.1, 1.1],
             [2, 2, 2]]
        b = [[0.1, 0.1, 0.1],
             [1.1, 1.1, 1.1]]
        with pytest.raises(ValueError):
            align.rotation_matrix(a, b)


class TestGetMatchingAtoms(object):
    @staticmethod
    @pytest.fixture()
    def universe():
        return mda.Universe(PSF, DCD)

    @staticmethod
    @pytest.fixture()
    def reference():
        return mda.Universe(PSF, DCD)

    @staticmethod
    @pytest.fixture()
    def reference_small(reference):
        return mda.Merge(reference.select_atoms(
            "not name H* and not atom 4AKE 1 CA"))

    @pytest.mark.parametrize("strict", (True, False))
    def test_match(self, universe, reference, strict,
                   selection="protein and backbone"):
        ref = reference.select_atoms(selection)
        mobile = universe.select_atoms(selection)
        groups = align.get_matching_atoms(ref, mobile, strict=strict)
        assert_equal(groups[0].names, groups[1].names)

    @pytest.mark.parametrize("strict", (True, False))
    def test_nomatch_atoms_raise(self, universe, reference,
                                 strict, selection="protein and backbone"):
        # one atom less but same residues; with strict=False should try
        # to get selections (but current code fails, so we also raise SelectionError)
        ref = reference.select_atoms(selection).atoms[1:]
        mobile = universe.select_atoms(selection)
        if strict:
            with pytest.raises(SelectionError):
                groups = align.get_matching_atoms(ref, mobile, strict=strict)
        else:
            with pytest.warns(SelectionWarning):
                with pytest.raises(SelectionError):
                    groups = align.get_matching_atoms(ref, mobile, strict=strict)

    @pytest.mark.parametrize("strict", (True, False))
    def test_nomatch_residues_raise_empty(self, universe, reference_small,
                                          strict, selection="protein and backbone"):
        # one atom less and all residues different: will currently create
        # empty selections with strict=False, see also
        # https://gist.github.com/orbeckst/2686badcd15031e6c946baf9164a683d
        ref = reference_small.select_atoms(selection)
        mobile = universe.select_atoms(selection)
        if strict:
            with pytest.raises(SelectionError):
                groups = align.get_matching_atoms(ref, mobile, strict=strict)
        else:
            with pytest.warns(SelectionWarning):
                with pytest.raises(SelectionError):
                    groups = align.get_matching_atoms(ref, mobile, strict=strict)

    def test_toggle_atom_mismatch_default_error(self, universe, reference):
        selection = ('resname ALA and name CA', 'resname ALA and name O')
        with pytest.raises(SelectionError):
            rmsd = align.alignto(universe, reference, select=selection)

    def test_toggle_atom_mismatch_kwarg_error(self, universe, reference):
        selection = ('resname ALA and name CA', 'resname ALA and name O')
        with pytest.raises(SelectionError):
            rmsd = align.alignto(universe, reference, select=selection, match_atoms=True)

    def test_toggle_atom_nomatch(self, universe, reference):
        selection = ('resname ALA and name CA', 'resname ALA and name O')
        rmsd = align.alignto(universe, reference, select=selection, match_atoms=False)
        assert rmsd[0] > 0.01

    def test_toggle_atom_nomatch_mismatch_atoms(self, universe, reference):
        # mismatching number of atoms, but same number of residues
        u = universe.select_atoms('resname ALA and name CA')
        u += universe.select_atoms('resname ALA and name O')[-1]
        ref = reference.select_atoms('resname ALA and name CA')
        with pytest.raises(SelectionError):
            align.alignto(u, ref, select='all', match_atoms=False)

    @pytest.mark.parametrize('subselection, expectation', [
        ('resname ALA and name CA', does_not_raise()),
        (mda.Universe(PSF, DCD).select_atoms('resname ALA and name CA'), does_not_raise()),
        (1234, pytest.raises(TypeError)),
    ])
    def test_subselection_alignto(self, universe, reference, subselection, expectation):

        with expectation:
            rmsd = align.alignto(universe, reference, subselection=subselection)
            assert_allclose(rmsd[1], 0.0, rtol=0, atol=1.5e-9)

    def test_no_atom_masses(self, universe):
        #if no masses are present
        u = mda.Universe.empty(6, 2, atom_resindex=[0, 0, 0, 1, 1, 1], trajectory=True)
        with pytest.warns(SelectionWarning):
            align.get_matching_atoms(u.atoms, u.atoms)

    def test_one_universe_has_masses(self, universe):
        u = mda.Universe.empty(6, 2, atom_resindex=[0, 0, 0, 1, 1, 1], trajectory=True)
        ref = mda.Universe.empty(6, 2, atom_resindex=[0, 0, 0, 1, 1, 1], trajectory=True)
        ref.add_TopologyAttr('masses')
        with pytest.warns(SelectionWarning):
            align.get_matching_atoms(u.atoms, ref.atoms)

class TestAlign(object):
    @staticmethod
    @pytest.fixture()
    def universe():
        return mda.Universe(PSF, DCD)

    @staticmethod
    @pytest.fixture()
    def reference():
        return mda.Universe(PSF, DCD)

    def test_rmsd(self, universe, reference):
        universe.trajectory[0]  # ensure first frame
        bb = universe.select_atoms('backbone')
        first_frame = bb.positions
        universe.trajectory[-1]
        last_frame = bb.positions
        assert_allclose(rms.rmsd(first_frame, first_frame), 0.0, rtol=0, atol=1.5e-5,
                        err_msg="error: rmsd(X,X) should be 0")
        # rmsd(A,B) = rmsd(B,A) should be exact but spurious failures in the
        # 9th decimal have been observed (see Issue 57 comment #1) so we relax
        # the test to 6 decimals.
        rmsd = rms.rmsd(first_frame, last_frame, superposition=True)
        assert_allclose(
            rms.rmsd(last_frame, first_frame, superposition=True), rmsd, rtol=0, atol=1.5e-6,
            err_msg="error: rmsd() is not symmetric")
        assert_allclose(rmsd, 6.820321761927005, rtol=0, atol=1.5e-5,
                        err_msg="RMSD calculation between 1st and last AdK frame gave wrong answer")     
        # test masses as weights
        last_atoms_weight = universe.atoms.masses
        A = universe.trajectory[0]
        B = reference.trajectory[-1]
        rmsd = align.alignto(universe, reference, weights='mass')
        rmsd_sup_weight = rms.rmsd(A, B, weights=last_atoms_weight, center=True,
                                   superposition=True)
        assert_allclose(rmsd[1], rmsd_sup_weight, rtol=0, atol=1.5e-6)

    def test_rmsd_custom_mass_weights(self, universe, reference):
        last_atoms_weight = universe.atoms.masses
        A = universe.trajectory[0]
        B = reference.trajectory[-1]
        rmsd = align.alignto(universe, reference,
                             weights=reference.atoms.masses)
        rmsd_sup_weight = rms.rmsd(A, B, weights=last_atoms_weight, center=True,
                                   superposition=True)
        assert_allclose(rmsd[1], rmsd_sup_weight, rtol=0, atol=1.5e-6)

    def test_rmsd_custom_weights(self, universe, reference):
        weights = np.zeros(universe.atoms.n_atoms)
        ca = universe.select_atoms('name CA')
        weights[ca.indices] = 1
        rmsd = align.alignto(universe, reference, select='name CA')
        rmsd_weights = align.alignto(universe, reference, weights=weights)
        assert_allclose(rmsd[1], rmsd_weights[1], rtol=0, atol=1.5e-6)

    def test_AlignTraj_outfile_default(self, universe, reference, tmpdir):
        with tmpdir.as_cwd():
            reference.trajectory[-1]
            x = align.AlignTraj(universe, reference)
            try:
                assert os.path.basename(x.filename) == 'rmsfit_adk_dims.dcd'
            finally:
                x._writer.close()

    def test_AlignTraj_outfile_default_exists(
      self, universe, reference, tmpdir, client_AlignTraj
    ):
        reference.trajectory[-1]
        outfile = str(tmpdir.join('align_test.dcd'))
        align.AlignTraj(universe, reference, filename=outfile).run(
          **client_AlignTraj
        )
        fitted = mda.Universe(PSF, outfile)

        # ensure default file exists
        with mda.Writer(str(tmpdir.join("rmsfit_align_test.dcd")),
                        n_atoms=fitted.atoms.n_atoms) as w:
            w.write(fitted.atoms)

        with tmpdir.as_cwd():
            align.AlignTraj(fitted, reference)

            # we are careful now. The default does nothing
            with pytest.raises(IOError):
                align.AlignTraj(fitted, reference, force=False)

    def test_AlignTraj_step_works(
      self, universe, reference, tmpdir, client_AlignTraj
    ):
        reference.trajectory[-1]
        outfile = str(tmpdir.join('align_test.dcd'))
        # this shouldn't throw an exception
        align.AlignTraj(universe, reference, filename=outfile).run(
          step=10, **client_AlignTraj
        )

    def test_AlignTraj_deprecated_attribute(self, universe, reference, tmpdir, client_AlignTraj):
        reference.trajectory[-1]
        outfile = str(tmpdir.join('align_test.dcd'))
        x = align.AlignTraj(universe, reference, filename=outfile).run(stop=2, **client_AlignTraj)

        wmsg = "The `rmsd` attribute was deprecated in MDAnalysis 2.0.0"
        with pytest.warns(DeprecationWarning, match=wmsg):
            assert_equal(x.rmsd, x.results.rmsd)

    def test_AlignTraj(self, universe, reference, tmpdir, client_AlignTraj):
        reference.trajectory[-1]
        outfile = str(tmpdir.join('align_test.dcd'))
        x = align.AlignTraj(universe, reference, filename=outfile).run(**client_AlignTraj)
        fitted = mda.Universe(PSF, outfile)

        assert_allclose(x.results.rmsd[0], 6.9290, rtol=0, atol=1.5e-3)
        assert_allclose(x.results.rmsd[-1], 5.2797e-07, rtol=0, atol=1.5e-3)

        # RMSD against the reference frame
        # calculated on Mac OS X x86 with MDA 0.7.2 r689
        # VMD: 6.9378711
        self._assert_rmsd(reference, fitted, 0, 6.929083044751061)
        self._assert_rmsd(reference, fitted, -1, 0.0)

    def test_AlignTraj_weighted(self, universe, reference, tmpdir, client_AlignTraj):
        outfile = str(tmpdir.join('align_test.dcd'))
        x = align.AlignTraj(universe, reference,
                            filename=outfile, weights='mass').run(**client_AlignTraj)
        fitted = mda.Universe(PSF, outfile)
        assert_allclose(x.results.rmsd[0], 0, rtol=0, atol=1.5e-3)
        assert_allclose(x.results.rmsd[-1], 6.9033, rtol=0, atol=1.5e-3)

        self._assert_rmsd(reference, fitted, 0, 0.0,
                          weights=universe.atoms.masses)
        self._assert_rmsd(reference, fitted, -1, 6.929083032629219,
                          weights=universe.atoms.masses)

    def test_AlignTraj_custom_weights(self, universe, reference, tmpdir, client_AlignTraj):
        weights = np.zeros(universe.atoms.n_atoms)
        ca = universe.select_atoms('name CA')
        weights[ca.indices] = 1

        outfile = str(tmpdir.join('align_test.dcd'))

        x = align.AlignTraj(universe, reference,
                            filename=outfile, select='name CA').run(**client_AlignTraj)
        x_weights = align.AlignTraj(universe, reference,
                                    filename=outfile, weights=weights).run(**client_AlignTraj)

        assert_allclose(x.results.rmsd, x_weights.results.rmsd, rtol=0, atol=1.5e-7)

    def test_AlignTraj_custom_mass_weights(self, universe, reference, tmpdir, client_AlignTraj):
        outfile = str(tmpdir.join('align_test.dcd'))
        x = align.AlignTraj(universe, reference,
                            filename=outfile,
                            weights=reference.atoms.masses).run(**client_AlignTraj)
        fitted = mda.Universe(PSF, outfile)
        assert_allclose(x.results.rmsd[0], 0, rtol=0, atol=1.5e-3)
        assert_allclose(x.results.rmsd[-1], 6.9033, rtol=0, atol=1.5e-3)

        self._assert_rmsd(reference, fitted, 0, 0.0,
                          weights=universe.atoms.masses)
        self._assert_rmsd(reference, fitted, -1, 6.929083032629219,
                          weights=universe.atoms.masses)

    def test_AlignTraj_partial_fit(self, universe, reference, tmpdir, client_AlignTraj):
        outfile = str(tmpdir.join('align_test.dcd'))
        # fitting on a partial selection should still write the whole topology
        align.AlignTraj(universe, reference, select='resid 1-20',
                        filename=outfile, weights='mass').run(**client_AlignTraj)
        mda.Universe(PSF, outfile)

    def test_AlignTraj_in_memory(self, universe, reference, tmpdir, client_AlignTraj):
        outfile = str(tmpdir.join('align_test.dcd'))
        reference.trajectory[-1]
        x = align.AlignTraj(universe, reference, filename=outfile,
                            in_memory=True).run(**client_AlignTraj)
        assert x.filename is None
        assert_allclose(x.results.rmsd[0], 6.9290, rtol=0, atol=1.5e-3)
        assert_allclose(x.results.rmsd[-1], 5.2797e-07, rtol=0, atol=1.5e-3)

        # check in memory trajectory
        self._assert_rmsd(reference, universe, 0, 6.929083044751061)
        self._assert_rmsd(reference, universe, -1, 0.0)

    def test_AlignTraj_writer_kwargs(self, universe, reference, tmpdir, client_AlignTraj):
        # Issue 4564
        writer_kwargs = dict(precision=2)
        with tmpdir.as_cwd():
            aligner = align.AlignTraj(universe, reference,
                                      select='protein and name CA',
                                      filename='aligned_traj.xtc',
                                      writer_kwargs=writer_kwargs,
                                      in_memory=False).run(**client_AlignTraj)
            assert_equal(aligner._writer.precision, 2)

    def _assert_rmsd(self, reference, fitted, frame, desired, weights=None):
        fitted.trajectory[frame]
        rmsd = rms.rmsd(reference.atoms.positions, fitted.atoms.positions,
                        superposition=True)
        assert_allclose(rmsd, desired, rtol = 0, atol=1.5e-5,
                            err_msg="frame {0:d} of fit does not have "
                                    "expected RMSD".format(frame))

    def test_alignto_checks_selections(self, universe, reference):
        """Testing that alignto() fails if selections do not
        match (Issue 143)"""
        u = universe

        def different_size():
            a = u.atoms[10:100]
            b = u.atoms[10:101]
            return align.alignto(a, b)

        with pytest.raises(SelectionError):
            different_size()

        def different_atoms():
            a = u.atoms[10:20]
            b = u.atoms[10:17] + u.atoms[18:21]
            return align.alignto(a, b)

        with pytest.raises(SelectionError):
            different_atoms()

    def test_alignto_partial_universe(self, universe, reference):
        u_bound = mda.Universe(ALIGN_BOUND)
        u_free = mda.Universe(ALIGN_UNBOUND)
        selection = 'segid B'

        segB_bound = u_bound.select_atoms(selection)
        segB_free = u_free.select_atoms(selection)
        segB_free.translate(segB_bound.centroid() - segB_free.centroid())

        align.alignto(u_free, u_bound, select=selection)
        assert_allclose(segB_bound.positions, segB_free.positions,
                                  rtol=0, atol=1.5e-3)


def _get_aligned_average_positions(ref_files, ref, select="all", **kwargs, client_AlignTraj):
    u = mda.Universe(*ref_files, in_memory=True)
    prealigner = align.AlignTraj(u, ref, select=select, **kwargs).run(**client_AlignTraj)
    ag = u.select_atoms(select)
    reference_coordinates = u.trajectory.timeseries(asel=ag).mean(axis=1)
    rmsd = sum(prealigner.results.rmsd/len(u.trajectory))
    return reference_coordinates, rmsd

class TestAverageStructure(object):

    ref_files = (PSF, DCD)

    @pytest.fixture
    def universe(self):
        return mda.Universe(*self.ref_files)

    @pytest.fixture
    def reference(self):
        return mda.Universe(PSF, CRD)

    def test_average_structure_deprecated_attrs(self, universe, reference, client_AverageStructure):
        # Issue #3278 - remove in MDAnalysis 3.0.0
        avg = align.AverageStructure(universe, reference).run(stop=2, **client_AverageStructure)

        wmsg = "The `universe` attribute was deprecated in MDAnalysis 2.0.0"
        with pytest.warns(DeprecationWarning, match=wmsg):
            assert_equal(avg.universe.atoms.positions,
                         avg.results.universe.atoms.positions)

        wmsg = "The `positions` attribute was deprecated in MDAnalysis 2.0.0"
        with pytest.warns(DeprecationWarning, match=wmsg):
            assert_equal(avg.positions, avg.results.positions)

        wmsg = "The `rmsd` attribute was deprecated in MDAnalysis 2.0.0"
        with pytest.warns(DeprecationWarning, match=wmsg):
            assert avg.rmsd == avg.results.rmsd

    def test_average_structure(self, universe, reference, client_AverageStructure):
        ref, rmsd = _get_aligned_average_positions(self.ref_files, reference)
        avg = align.AverageStructure(universe, reference).run(**client_AverageStructure)
        assert_allclose(avg.results.universe.atoms.positions, ref, rtol=0, atol=1.5e-4)
        assert_allclose(avg.results.rmsd, rmsd, rtol=0, atol=1.5e-7)

    def test_average_structure_mass_weighted(self, universe, reference, client_AverageStructure):
        ref, rmsd = _get_aligned_average_positions(self.ref_files, reference, weights='mass')
        avg = align.AverageStructure(universe, reference, weights='mass').run(**client_AverageStructure)
        assert_allclose(avg.results.universe.atoms.positions, ref,
                            rtol=0, atol=1.5e-4)
        assert_allclose(avg.results.rmsd, rmsd, rtol=0, atol=1.5e-7)

    def test_average_structure_select(self, universe, reference, client_AverageStructure):
        select = 'protein and name CA and resid 3-5'
        ref, rmsd = _get_aligned_average_positions(self.ref_files, reference, select=select)
        avg = align.AverageStructure(universe, reference, select=select).run(**client_AverageStructure)
        assert_allclose(avg.results.universe.atoms.positions, ref,
                            rtol=0, atol=1.5e-4)
        assert_allclose(avg.results.rmsd, rmsd, rtol=0, atol=1.5e-7)

    def test_average_structure_no_ref(self, universe, client_AverageStructure):
        ref, rmsd = _get_aligned_average_positions(self.ref_files, universe)
        avg = align.AverageStructure(universe).run(**client_AverageStructure)
        assert_allclose(avg.results.universe.atoms.positions, ref,
                        rtol=0, atol=1.5e-4)
        assert_allclose(avg.results.rmsd, rmsd, rtol=0, atol=1.5e-7)

    def test_average_structure_no_msf(self, universe, client_AverageStructure):
        avg = align.AverageStructure(universe).run(**client_AverageStructure)
        assert not hasattr(avg, 'msf')

    def test_mismatch_atoms(self, universe):
        u = mda.Merge(universe.atoms[:10])
        with pytest.raises(SelectionError):
            align.AverageStructure(universe, u)

    def test_average_structure_ref_frame(self, universe, client_AverageStructure):
        ref_frame = 3
        u = mda.Merge(universe.atoms)

        # change to ref_frame
        universe.trajectory[ref_frame]
        u.load_new(universe.atoms.positions)

        # back to start
        universe.trajectory[0]
        ref, rmsd = _get_aligned_average_positions(self.ref_files, u)
        avg = align.AverageStructure(universe, ref_frame=ref_frame).run(**client_AverageStructure)
        assert_allclose(avg.results.universe.atoms.positions, ref,
                        rtol=0, atol=1.5e-4)
        assert_allclose(avg.results.rmsd, rmsd, rtol=0, atol=1.5e-7)

    def test_average_structure_in_memory(self, universe, client_AverageStructure):
        avg = align.AverageStructure(universe, in_memory=True).run(**client_AverageStructure)
        reference_coordinates = universe.trajectory.timeseries().mean(axis=1)
        assert_allclose(avg.results.universe.atoms.positions,
                        reference_coordinates, rtol=0, atol=1.5e-4)
        assert avg.filename is None


class TestAlignmentProcessing:
    seq = FASTA
    error_msg = "selection string has unexpected length"

    @pytest.mark.skipif(HAS_BIOPYTHON, reason='biopython is installed')
    def test_importerror_biopython(self):
        errmsg = "The `fasta2select` method requires an installation"
        with pytest.raises(ImportError, match=errmsg):
            _ = align.fasta2select(self.seq, is_aligned=True)

    @pytest.mark.skipif(not HAS_BIOPYTHON, reason='requires biopython')
    def test_fasta2select_aligned(self):
        """test align.fasta2select() on aligned FASTA (Issue 112)"""
        sel = align.fasta2select(self.seq, is_aligned=True)
        # length of the output strings, not residues or anything real...
        assert len(sel['reference']) == 30623, self.error_msg
        assert len(sel['mobile']) == 30623, self.error_msg

    @pytest.mark.skipif(
        executable_not_found("clustalw2") or not HAS_BIOPYTHON,
        reason="Test skipped because clustalw2 executable not found")
    def test_fasta2select_file(self, tmpdir):
        """test align.fasta2select() on a non-aligned FASTA with default
        filenames"""
        with tmpdir.as_cwd():
            sel = align.fasta2select(self.seq, is_aligned=False,
                                     alnfilename=None, treefilename=None)
            assert len(sel['reference']) == 23080, self.error_msg
            assert len(sel['mobile']) == 23090, self.error_msg

    @pytest.mark.skipif(
        executable_not_found("clustalw2") or not HAS_BIOPYTHON,
        reason="Test skipped because clustalw2 executable not found")
    def test_fasta2select_ClustalW(self, tmpdir):
        """MDAnalysis.analysis.align: test fasta2select() with ClustalW
        (Issue 113)"""
        alnfile = str(tmpdir.join('alignmentprocessing.aln'))
        treefile = str(tmpdir.join('alignmentprocessing.dnd'))
        sel = align.fasta2select(self.seq, is_aligned=False,
                                 alnfilename=alnfile, treefilename=treefile)
        # numbers computed from alignment with clustalw 2.1 on Mac OS X
        # [orbeckst] length of the output strings, not residues or anything
        # real...
        assert len(sel['reference']) == 23080, self.error_msg
        assert len(sel['mobile']) == 23090, self.error_msg

    @pytest.mark.skipif(not HAS_BIOPYTHON, reason='requires biopython')
    def test_fasta2select_resids(self, tmpdir):
        """test align.fasta2select() when resids provided (Issue #3124)"""
        resids = [x for x in range(705)]
        sel = align.fasta2select(self.seq, is_aligned=True,
                                 ref_resids=resids, target_resids=resids)
        # length of the output strings, not residues or anything real...
        assert len(sel['reference']) == 30621, self.error_msg
        assert len(sel['mobile']) == 30621, self.error_msg


class TestSequenceAlignmentFunction:
    # remove 3.0

    @staticmethod
    @pytest.fixture
    def atomgroups():
        universe = mda.Universe(PSF)
        reference = universe.atoms
        mobile = universe.select_atoms("resid 122-159")
        return reference, mobile

    @pytest.mark.skipif(HAS_BIOPYTHON, reason='biopython installed')
    def test_biopython_import_error(self, atomgroups):
        ref, mob = atomgroups
        errmsg = "The `sequence_alignment` method requires an installation of"
        with pytest.raises(ImportError, match=errmsg):
            align.sequence_alignment(mob, ref)

    @pytest.mark.skipif(not HAS_BIOPYTHON, reason='requires biopython')
    @pytest.mark.filterwarnings("ignore:`sequence_alignment` is deprecated!")
    def test_sequence_alignment(self, atomgroups):
        reference, mobile = atomgroups
        aln = align.sequence_alignment(mobile, reference)

        assert len(aln) == 5, "return value has wrong tuple size"

        seqA, seqB, score, begin, end = aln
        assert_equal(seqA, reference.residues.sequence(format="string"),
                     err_msg="reference sequence mismatch")
        assert mobile.residues.sequence(
            format="string") in seqB, "mobile sequence mismatch"
        assert score  == pytest.approx(54.6)
        assert_array_equal([begin, end], [0, reference.n_residues])

    @pytest.mark.skipif(not HAS_BIOPYTHON, reason='requires biopython')
    def test_sequence_alignment_deprecation(self, atomgroups):
        reference, mobile = atomgroups
        wmsg = ("`sequence_alignment` is deprecated!\n"
                "`sequence_alignment` will be removed in release 3.0.")
        with pytest.warns(DeprecationWarning, match=wmsg):
            align.sequence_alignment(mobile, reference)


class TestIterativeAverage(object):
    @pytest.fixture()
    def mobile(self):
        u = mda.Universe(PSF, DCD)
        return u

    @pytest.fixture()
    def reference(self):
        u = mda.Universe(PSF, DCD)
        return u

    def test_iterative_average_default(self, mobile):
        res = align.iterative_average(mobile, select="bynum 1:10")
        assert_allclose(
            res.results.positions,
            [
                [11.93075595, 8.6729893, -10.49887605],
                [12.60587898, 7.91673117, -10.73327464],
                [12.45662411, 9.51900517, -10.35551193],
                [11.27452274, 8.83003843, -11.2619057],
                [11.25808119, 8.26794477, -9.23340715],
                [12.02767222, 7.95332228, -8.57249317],
                [10.54679871, 9.49505306, -8.61215292],
                [9.99500556, 9.16624224, -7.75231192],
                [9.83897407, 9.93134598, -9.29541129],
                [11.45760169, 10.5857071, -8.13037669]
            ],
            atol=1e-5,
        )

    def test_iterative_average_eps_high(self, mobile):
        res = align.iterative_average(mobile, select="bynum 1:10",
                                      eps=1e-5)
        assert_allclose(
            res.results.positions,
            [
                [11.93075595, 8.6729893, -10.49887605],
                [12.60587898, 7.91673117, -10.73327464],
                [12.45662411, 9.51900517, -10.35551193],
                [11.27452274, 8.83003843, -11.2619057],
                [11.25808119, 8.26794477, -9.23340715],
                [12.02767222, 7.95332228, -8.57249317],
                [10.54679871, 9.49505306, -8.61215292],
                [9.99500556, 9.16624224, -7.75231192],
                [9.83897407, 9.93134598, -9.29541129],
                [11.45760169, 10.5857071, -8.13037669]
            ],
            atol=1e-5,
        )

    def test_iterative_average_weights_mass(self, mobile, reference):
        res = align.iterative_average(mobile, reference,
                                      select="bynum 1:10",
                                      niter=10, weights="mass")
        assert_allclose(
            res.results.positions,
            [
                [11.96438784, 8.85426235, -10.24735737],
                [12.75920431, 8.27294545, -10.54295766],
                [12.3285704, 9.72083717, -9.9419435],
                [11.33941507, 9.03249423, -11.01306158],
                [11.30988499, 8.14958885, -9.1205501],
                [12.09108655, 7.85155906, -8.46681943],
                [10.37499697, 9.13535837, -8.3732586],
                [9.83883314, 8.57939098, -7.6195549],
                [9.64405257, 9.55924307, -9.04315991],
                [11.0678934, 10.27798773, -7.64881842]
            ],
            atol=1e-5,
        )

    def test_iterative_average_convergence_failure(self, mobile, reference):
        with pytest.raises(RuntimeError):
            _ = align.iterative_average(mobile, reference,
                                        niter=1, eps=0)

    def test_iterative_average_convergence_verbose(self, mobile, reference):
        _ = align.iterative_average(mobile, select="bynum 1:10",
                                    verbose=True)


def test_alignto_reorder_atomgroups():
    # Issue 2977
    u = mda.Universe(PDB_helix)
    mobile = u.atoms[:4]
    ref = u.atoms[[3, 2, 1, 0]]
    rmsd = align.alignto(mobile, ref, select='bynum 1-4')
    assert_allclose(rmsd, (0.0, 0.0))
