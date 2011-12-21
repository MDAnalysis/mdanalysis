# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

from Scientific.Geometry import Transformation

class Transformation:
    def __init__(self):
        pass
    def transform(self):
        raise NotImplemented()

class Recenter(Transformation):
    def __init__(self, system, asel):
        self.system = system
        self.asel = asel
    def transform(self):
        com = self.asel.centerOfMass()
        self.system.coord -= self.asel.centerOfMass()

class RMSOrient(Transformation):
    def __init__(self, system, asel):
        self.system = system
        self.asel = asel
    #def transform(self):
    #    # XXX Not complete yet
    #    com = self.asel.centerOfMass()

