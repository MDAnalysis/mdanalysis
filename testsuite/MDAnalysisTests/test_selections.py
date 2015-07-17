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

# Test the selection exporters in MDAnalysis.selections

# use cStringIO and NamedStream to write to memory instead to temp files
import cStringIO

import numpy as np
from numpy.testing import TestCase, assert_equal, assert_array_equal
from nose.plugins.attrib import attr

from MDAnalysisTests.plugins.knownfailure import knownfailure
from MDAnalysis.lib.util import NamedStream

import MDAnalysis
from MDAnalysis.tests.datafiles import PSF, DCD, PDB

from re import search

# add more tests, see Issue #353

class _SelectionWriter(TestCase):
    filename = None

    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        stream = cStringIO.StringIO()
        self.namedfile = NamedStream(stream, self.filename)

    def tearDown(self):
        del self.universe
        del self.namedfile

def ndx2array(lines):
    """Convert Gromacs NDX text file lines to integer array"""
    return np.array(" ".join(lines).replace("\n", "").split(), dtype=int)

class TestSelectionWriter_Gromacs(_SelectionWriter):
    filename = "CA.ndx"

    ref_name = "CA_selection"
    ref_indices = ndx2array(
        [ '5 22 46 65 84 103 122 129 141 153 160 170 \n',
          '177 199 206 220 237 247 264 284 303 320 335 357 \n',
          '378 385 406 418 435 454 465 479 486 498 515 534 \n',
          '558 568 578 594 616 627 634 645 660 679 686 708 \n',
          '725 735 757 769 788 805 817 827 834 856 875 891 \n',
          '905 917 932 951 967 986 996 1015 1031 1053 1068 1092 \n',
          '1111 1121 1138 1153 1165 1176 1200 1214 1221 1241 1260 1279 \n',
          '1291 1298 1320 1332 1356 1370 1391 1403 1420 1430 1442 1452 \n',
          '1469 1491 1506 1516 1523 1542 1556 1572 1584 1605 1621 1640 \n',
          '1655 1675 1687 1705 1717 1729 1744 1763 1782 1798 1810 1834 \n',
          '1853 1869 1876 1900 1924 1940 1957 1969 1981 1992 1999 2023 \n',
          '2039 2060 2077 2093 2115 2135 2151 2165 2177 2199 2215 2230 \n',
          '2237 2259 2271 2283 2299 2313 2320 2335 2350 2369 2383 2397 \n',
          '2421 2443 2455 2467 2484 2499 2514 2528 2544 2568 2590 2614 \n',
          '2633 2649 2664 2685 2702 2719 2736 2750 2762 2774 2793 2812 \n',
          '2819 2840 2861 2872 2894 2909 2919 2934 2944 2951 2965 2979 \n',
          '3001 3022 3032 3054 3070 3082 3089 3103 3127 3139 3155 3165 \n',
          '3180 3196 3220 3230 3242 3261 3276 3298 3317 3336 \n']
        )

    def test_write_ndx(self):
        CA = self.universe.selectAtoms("protein and name CA")
        CA.write(self.namedfile, name=self.ref_name)

        header = self.namedfile.readline().strip()
        assert_equal(header, "[ {0} ]".format(self.ref_name),
                     err_msg="NDX file has wrong selection name")
        indices = ndx2array(self.namedfile.readlines())
        assert_array_equal(indices, self.ref_indices,
                           err_msg="indices were not written correctly")

    def test_writeselection_ndx(self):
        CA = self.universe.selectAtoms("protein and name CA")
        CA.write_selection(self.namedfile, name=self.ref_name)

        header = self.namedfile.readline().strip()
        assert_equal(header, "[ {0} ]".format(self.ref_name),
                     err_msg="NDX file has wrong selection name")
        indices = ndx2array(self.namedfile.readlines())
        assert_array_equal(indices, self.ref_indices,
                           err_msg="indices were not written correctly")

def spt2array(line):
    """Get name of and convert Jmol SPT definition to integer array"""
    match = search(r'\@~(\w+) \(\{([\d\s]*)\}\)', line)
    return match.group(1), np.array(match.group(2).split(), dtype=int)

class TestSelectionWriter_Jmol(_SelectionWriter):
    filename = "CA.spt"

    def setUp(self):
        self.universe = MDAnalysis.Universe(PDB)
        stream = cStringIO.StringIO()
        self.namedfile = NamedStream(stream, self.filename)

    ref_name, ref_indices = spt2array(
        ( '@~ca ({4 21 45 64 83 102 121 128 137 152 159 169 176 198 205 219 236'
          ' 246 263 283 302 319 334 356 377 384 402 417 434 453 464 478 485 497'
          ' 514 533 557 567 577 593 615 626 633 644 659 678 685 707 724 734 756'
          ' 768 787 804 816 826 833 855 874 890 904 916 931 950 966 985 995 101'
          '4 1030 1052 1067 1091 1110 1120 1137 1152 1164 1175 1199 1213 1220 1'
          '240 1259 1278 1290 1297 1316 1331 1355 1369 1387 1402 1419 1429 1441'
          ' 1451 1468 1490 1505 1515 1522 1541 1555 1571 1583 1604 1620 1639 16'
          '54 1674 1686 1701 1716 1728 1743 1762 1781 1797 1809 1833 1852 1868 '
          '1875 1899 1923 1939 1956 1965 1980 1991 1998 2022 2038 2059 2076 209'
          '2 2114 2134 2147 2161 2176 2198 2214 2229 2236 2258 2270 2282 2298 2'
          '312 2319 2334 2349 2368 2382 2396 2420 2442 2454 2466 2483 2498 2513'
          ' 2527 2543 2567 2589 2613 2632 2648 2663 2684 2701 2718 2735 2749 27'
          '58 2773 2792 2811 2818 2839 2860 2871 2893 2908 2918 2933 2943 2950 '
          '2964 2978 3000 3021 3031 3053 3069 3081 3088 3102 3123 3138 3154 316'
          '4 3179 3195 3219 3229 3241 3260 3275 3297 3316 3335});')
        )

    def test_write_spt(self):
        CA = self.universe.selectAtoms("protein and name CA")
        CA.write(self.namedfile, name=self.ref_name)

        header, indices = spt2array(self.namedfile.readline())

        assert_equal(header, self.ref_name,
                     err_msg="SPT file has wrong selection name")
        assert_array_equal(indices, self.ref_indices,
                           err_msg="indices were not written correctly")

    def test_writeselection_spt(self):
        CA = self.universe.selectAtoms("protein and name CA")
        CA.write_selection(self.namedfile, name=self.ref_name)

        header, indices = spt2array(self.namedfile.readline())

        assert_equal(header, self.ref_name,
                     err_msg="SPT file has wrong selection name")
