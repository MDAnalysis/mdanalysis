# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2020 The MDAnalysis Development Team and contributors
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

"""HOLE Analysis --- :mod:`MDAnalysis.analysis.hole2.templates`
=====================================================================================

:Author: Lily Wang
:Year: 2020
:Copyright: GNU Public License v3

.. versionadded:: 1.0

Templates used in :mod:`MDAnalysis.analysis.hole2.hole`
"""

exe_err = ('HOLE binary {name} not found. {name} must be on the '
           'PATH, or the path must provided with the keyword '
           'argument: {kw}')

IGNORE_RESIDUES = ["SOL", "WAT", "TIP", "HOH", "K  ", "NA ", "CL "]


#: Built-in HOLE radii (based on ``simple.rad`` from the HOLE_ distribution):
#: van der Waals radii are AMBER united atom from Weiner et al. (1984), JACS, vol 106 pp765-768.
#: *Simple* - Only use one value for each element C O H etc.
#: Added radii for K+, NA+, CL- (Pauling hydration radius from Hille 2002).
#: The data file can be written with the convenience function :func:`write_simplerad2`.
SIMPLE2_RAD = r"""
remark: Time-stamp: <2005-11-21 13:57:55 oliver> [OB]
remark: van der Waals radii: AMBER united atom
remark: from Weiner et al. (1984), JACS, vol 106 pp765-768
remark: Simple - Only use one value for each element C O H etc.
remark: van der Waals radii
remark: general last
VDWR C??? ??? 1.85
VDWR O??? ??? 1.65
VDWR S??? ??? 2.00
VDWR N??? ??? 1.75
VDWR H??? ??? 1.00
VDWR H?   ??? 1.00
VDWR P??? ??? 2.10
remark: ASN, GLN polar H (odd names for these atoms in xplor)
VDWR E2?  GLN 1.00
VDWR D2?  ASN 1.00
remark: amber lone pairs on sulphurs
VDWR LP?? ??? 0.00
remark: for some funny reason it wants radius for K even though
remark: it is on the exclude list
remark: Use Pauling hydration radius (Hille 2001) [OB]
VDWR K?   ??? 1.33
VDWR NA?  ??? 0.95
VDWR CL?  ??? 1.81
remark: funny hydrogens in gA structure [OB]
VDWR 1H?  ??? 1.00
VDWR 2H?  ??? 1.00
VDWR 3H?  ??? 1.00
remark: need bond rad for molqpt option
BOND C??? 0.85
BOND N??? 0.75
BOND O??? 0.7
BOND S??? 1.1
BOND H??? 0.5
BOND P??? 1.0
BOND ???? 0.85
"""

hole_input = """
! Input file for Oliver Smart's HOLE program
! written by MDAnalysis.analysis.hole2.hole
! filename = {filename}
COORD  {coordinates}
RADIUS {radius}
SPHPDB {sphpdb}
SAMPLE {sample:f}
ENDRAD {end_radius:f}
IGNORE {ignore}
SHORTO {output_level:d}
"""

hole_lines = {
    'cpoint': 'CPOINT {:.10f} {:.10f} {:.10f}\n',
    'cvect': 'CVECT {:.10f} {:.10f} {:.10f}\n',
    'random_seed': 'RASEED {}\n',
    'dcd': 'CHARMD {dcd}\nCHARMS {iniskip:d} {step:d}\n',
}

vmd_script_array = """\
set no_water_color {no_water_color}
set one_water_color {one_water_color}
set double_water_color {double_water_color}
array set triangles {{}}
"""

vmd_script_function = r"""
global vmd_frame;
trace add variable vmd_frame([molinfo top]) write drawFrame

proc drawFrame { name element op } {
            global vmd_frame triangles no_water_color one_water_color double_water_color;
            set frame $vmd_frame([molinfo top])
            draw delete all;

            draw color $no_water_color;
            foreach shape [lindex $triangles($frame) 0] {
                draw trinorm {*}$shape
            }
            draw color $one_water_color;
            foreach shape [lindex $triangles($frame) 1] {
                draw trinorm {*}$shape
            }
            draw color $double_water_color;
            foreach shape [lindex $triangles($frame) 2] {
                draw trinorm {*}$shape
            }
        }

drawFrame 0 0 0
"""

