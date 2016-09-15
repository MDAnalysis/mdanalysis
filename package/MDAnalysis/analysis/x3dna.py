# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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

# stub for legacy.x3dna
# See https://github.com/MDAnalysis/mdanalysis/issues/906

# Remove this file in version 1.0

__all__ = ["mean_std_from_x3dnaPickle", "X3DNA", "X3DNAtraj"]

import warnings

with warnings.catch_warnings():
    warnings.simplefilter('always', DeprecationWarning)
    warnings.warn(("The MDAnalysis.analysis.x3dna module is not fully "
                   "maintained anymore and was moved to the legacy "
                   "module: To import it use "
                   "'from MDAnalysis.analysis.legacy import x3dna'."
                   "The MDAnalysis.analysis.x3dna module will be "
                   "removed in the 1.0 release."),
                  DeprecationWarning)


from MDAnalysis.analysis.legacy.x3dna import (
    mean_std_from_x3dnaPickle, X3DNA, X3DNAtraj, BaseX3DNA)
