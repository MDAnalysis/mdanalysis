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
from __future__ import absolute_import

import pytest
import parmed as pmd
import MDAnalysis as mda

from numpy.testing import (assert_equal,
                           assert_array_almost_equal,
                           assert_almost_equal)

from MDAnalysisTests.coordinates.base import _SingleFrameReader
from MDAnalysisTests.coordinates.reference import RefAdKSmall

from MDAnalysisTests.datafiles import (
    GRO
)

class TestParmEdReaderGRO:
    __test__ = True
    ref_filename = GRO

    universe = mda.Universe(pmd.load_file(GRO))
    ref = mda.Universe(GRO)
    prec = 3

    def test_dimensions(self):
        assert_array_almost_equal(
            self.universe.trajectory.ts.dimensions, self.ref.trajectory.ts.dimensions,
            self.prec,
            "ParmEdReader failed to get unitcell dimensions from ParmEd")
    
    def test_coordinates(self):
        up = self.universe.atoms.positions
        rp = self.ref.atoms.positions
        assert_array_almost_equal(up, rp, decimal=3)
    


class BaseTestParmEdReader(_SingleFrameReader):

    def setUp(self):
        self.universe = mda.Universe(pmd.load_file(self.ref_filename))
        self.ref = mda.Universe(self.ref_filename)
        self.prec = 3

    def test_dimensions(self):
        assert_array_almost_equal(
            self.universe.trajectory.ts.dimensions, self.ref.trajectory.ts.dimensions,
            self.prec,
            "ParmEdReader failed to get unitcell dimensions from ParmEd")
    
    def test_coordinates(self):
        up = self.universe.atoms.positions
        rp = self.ref.atoms.positions
        assert_array_almost_equal(up, rp, decimal=3)
    

class TestParmEdReaderPDB(BaseTestParmEdReader):
    __test__ = True

    ref_filename = RefAdKSmall.filename
    
    def test_uses_ParmEdReader(self):
        from MDAnalysis.coordinates.ParmEd import ParmEdReader

        assert isinstance(self.universe.trajectory, ParmEdReader), "failed to choose ParmEdReader"




class BaseTestParmEdWriter:

    prec = 0.01

    def __init__(self):
        self.ref = pmd.load_file(self.ref_filename)
        self.universe = mda.Universe(self.ref)
        self.output = self.universe.atoms.write(file_format='PARMED')
    
    def test_equivalent_connectivity_counts(self):
        for attr in ('atoms', 'bonds', 'angles', 'dihedrals', 'impropers',
                     'cmaps', 'urey_bradleys'):
            r = getattr(self.ref, attr)
            o = getattr(self.output, attr)
            assert len(r) == len(o)
    
    def test_equivalent_connectivity_atom_counts(self):
        for i in range(len(self.ref.atoms)):
            assert len(self.ref.atoms[i].bonds) == len(self.output.atoms[i].bonds)

    def test_equivalent_coordinates(self):
        assert_array_almost_equal(self.ref.coordinates, self.output.coordinates)
    
    def test_equivalent_atoms(self):
        for r, o in zip(self.ref.atoms, self.output.atoms):
            for attr in ('name',
                         'number', 'altloc', 
                         'atomic_number'):
                ra = getattr(r, attr)
                oa = getattr(o, attr)
                assert ra == oa
            
            for attr in ('mass', 'charge', 'occupancy', 'rmin', 'epsilon'):
                ra = getattr(r, attr)
                oa = getattr(o, attr)
                assert abs(ra-oa) < self.prec

            
        

class TestParmEdReaderGRO(BaseTestParmEdWriter):
    ref_filename = GRO
