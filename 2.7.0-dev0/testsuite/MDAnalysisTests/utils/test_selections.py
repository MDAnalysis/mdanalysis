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

# Test the selection exporters in MDAnalysis.selections
# use StringIO and NamedStream to write to memory instead to temp files
import pytest
from io import StringIO
from unittest.mock import patch

import re

import numpy as np
from numpy.testing import assert_equal, assert_array_equal

from MDAnalysis.tests.datafiles import PSF, DCD

import MDAnalysis
import MDAnalysis.selections as selections
from MDAnalysis.lib.util import NamedStream


class _SelectionWriter(object):

    filename = None
    max_number = 357  # to keep fixtures smallish, only select CAs up to number 357

    @staticmethod
    @pytest.fixture()
    def universe():
        return MDAnalysis.Universe(PSF, DCD)

    @pytest.fixture()
    def namedfile(self):
        return NamedStream(StringIO(), self.filename)

    def _selection(self, universe):
        return universe.select_atoms("protein and name CA and bynum 1-{0}".format(self.max_number))

    def _write(self, universe, namedfile, **kwargs):
        g = self._selection(universe)
        g.write(namedfile, **kwargs)
        return g

    def _write_selection(self, universe, namedfile, **kwargs):
        g = self._selection(universe)
        g.write(namedfile, **kwargs)
        return g

    def _write_with(self, universe, namedfile, **kwargs):
        g = self._selection(universe)
        with self.writer(namedfile, **kwargs) as outfile:
            outfile.write(g)
        return g

    def test_write_bad_mode(self, universe, namedfile,):
        with pytest.raises(ValueError):
            self._write(universe, namedfile, name=self.ref_name, mode='a+')

    def test_write(self, universe, namedfile):
        self._write(universe, namedfile, name=self.ref_name)
        self._assert_selectionstring(namedfile)

    def test_writeselection(self, universe, namedfile):
        self._write_selection(universe, namedfile, name=self.ref_name)
        self._assert_selectionstring(namedfile)

    def test_write_with(self, universe, namedfile):
        self._write_with(universe, namedfile, name=self.ref_name)
        self._assert_selectionstring(namedfile)


def ndx2array(lines):
    """Convert Gromacs NDX text file lines to integer array"""
    return np.array(" ".join(lines).replace("\n", "").split(), dtype=int)


def lines2one(lines):
    """Join lines and squash all whitespace"""
    return " ".join(" ".join(lines).split())


class TestSelectionWriter_Gromacs(_SelectionWriter):

    writer = MDAnalysis.selections.gromacs.SelectionWriter
    filename = "CA.ndx"
    ref_name = "CA_selection"
    ref_indices = ndx2array(
        [ '5 22 46 65 84 103 122 129 141 153 160 170 \n',
          '177 199 206 220 237 247 264 284 303 320 335 357 \n',
          ]
        )

    def _assert_selectionstring(self, namedfile):
        header = namedfile.readline().strip()
        assert_equal(header, "[ {0} ]".format(self.ref_name),
                     err_msg="NDX file has wrong selection name")
        indices = ndx2array(namedfile.readlines())
        assert_array_equal(indices, self.ref_indices,
                           err_msg="indices were not written correctly")


class TestSelectionWriter_Charmm(_SelectionWriter):

    writer = MDAnalysis.selections.charmm.SelectionWriter
    filename = "CA.str"
    ref_name = "CA_selection"
    ref_selectionstring = lines2one([
        """! MDAnalysis CHARMM selection
           DEFINE CA_selection SELECT -
           BYNUM 5 .or. BYNUM 22 .or. BYNUM 46 .or. BYNUM 65 .or. -
           BYNUM 84 .or. BYNUM 103 .or. BYNUM 122 .or. BYNUM 129 .or. -
           BYNUM 141 .or. BYNUM 153 .or. BYNUM 160 .or. BYNUM 170 .or. -
           BYNUM 177 .or. BYNUM 199 .or. BYNUM 206 .or. BYNUM 220 .or. -
           BYNUM 237 .or. BYNUM 247 .or. BYNUM 264 .or. BYNUM 284 .or. -
           BYNUM 303 .or. BYNUM 320 .or. BYNUM 335 .or. BYNUM 357 END
        """])

    def _assert_selectionstring(self, namedfile):
        selectionstring = lines2one(namedfile.readlines())
        assert_equal(selectionstring, self.ref_selectionstring,
                     err_msg="Charmm selection was not written correctly")


class TestSelectionWriter_PyMOL(_SelectionWriter):

    writer = MDAnalysis.selections.pymol.SelectionWriter
    filename = "CA.pml"
    ref_name = "CA_selection"
    ref_selectionstring = lines2one([
        """# MDAnalysis PyMol selection\n select CA_selection, \\
           index 5 | index 22 | index 46 | index 65 | index 84 | index 103 | \\
           index 122 | index 129 | index 141 | index 153 | index 160 | index 170 | \\
           index 177 | index 199 | index 206 | index 220 | index 237 | index 247 | \\
           index 264 | index 284 | index 303 | index 320 | index 335 | index 357
        """])

    def _assert_selectionstring(self, namedfile):
        selectionstring = lines2one(namedfile.readlines())
        assert_equal(selectionstring, self.ref_selectionstring,
                     err_msg="PyMOL selection was not written correctly")


class TestSelectionWriter_VMD(_SelectionWriter):

    writer = MDAnalysis.selections.vmd.SelectionWriter
    filename = "CA.vmd"
    ref_name = "CA_selection"
    ref_selectionstring = lines2one([
        """# MDAnalysis VMD selection atomselect macro CA_selection {index 4 21 45 64 83 102 121 128 \\
           140 152 159 169 176 198 205 219 \\
           236 246 263 283 302 319 334 356 }
        """])

    def _assert_selectionstring(self, namedfile):
        selectionstring = lines2one(namedfile.readlines())
        assert_equal(selectionstring, self.ref_selectionstring,
                     err_msg="PyMOL selection was not written correctly")


def spt2array(line):
    """Get name of and convert Jmol SPT definition to integer array"""
    match = re.search(r'\@~(\w+) \(\{([\d\s]*)\}\)', line)
    return match.group(1), np.array(match.group(2).split(), dtype=int)


class TestSelectionWriter_Jmol(_SelectionWriter):

    writer = MDAnalysis.selections.jmol.SelectionWriter
    filename = "CA.spt"
    ref_name, ref_indices = spt2array(
        ( '@~ca ({4 21 45 64 83 102 121 128 140 152 159 169 176 198 205 219 236'
          ' 246 263 283 302 319 334 356});')
        )

    def _assert_selectionstring(self, namedfile):
        header, indices = spt2array(namedfile.readline())
        assert_equal(header, self.ref_name,
                     err_msg="SPT file has wrong selection name")
        assert_array_equal(indices, self.ref_indices,
                           err_msg="SPT indices were not written correctly")


class TestSelections:
    @patch.object(selections, "_SELECTION_WRITERS", {"FOO": "BAR"})
    def test_get_writer_valid(self):
        writer = selections.get_writer(None, "FOO")
        assert writer == "BAR"
