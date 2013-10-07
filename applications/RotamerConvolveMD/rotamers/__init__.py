# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# Convolve MTSS rotamers with MD trajectory.
# Copyright (c) 2011-2013 Philip Fowler, Oliver Beckstein
# Published under the GNU Public Licence, version 2 (or higher)
#
# Includes a rotamer library for MTSS at 298 K by Gunnar Jeschke,
# which is published under the same licence by permission.
"""\
Convolve Rotamers
=================

:Author:  Philip Fowler, Oliver Beckstein
:Year:    2011-2013
:Licence: GNU Public Licence, version 2 (or higher)

Calculate a distribution of spin label distances from an MD trajectory
or an arbitrary ensemble of conformations by fitting Gunnar Jeschke's
rotamer library of MTSS (at 298 K).

"""

VERSION = 1,0

