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

import MDAnalysis
import numpy
from numpy.testing import *
from MDAnalysis.tests.datafiles import GRO_velocity

class TestGROVelocities(TestCase):
    def setUp(self):
        #reference velocities for the full 6-atom test case:
        self.reference_velocities = numpy.array(
            [[-101.227    ,  -0.57999998,   0.43400002],
             [  8.08500004,   3.19099998,  -7.79099989],
             [ -9.04500008, -26.46899986,  13.17999935],
             [  2.51899981,   3.1400001 ,  -1.73399997],
             [-10.64100075, -11.34899998,   0.257     ],
             [ 19.42700005,  -8.21600056 ,  -0.24399999 ]], dtype=numpy.float32)
        self.prec = 3

    def testParse_velocities(self):
        #read the velocities from the GRO_velocity file and compare the AtomGroup and individual Atom velocities parsed with the reference values:
        u = MDAnalysis.Universe(GRO_velocity)
        all_atoms = u.selectAtoms('all')
        #check for read-in and unit conversion for .gro file velocities for the entire AtomGroup:
        assert_almost_equal(all_atoms.velocities(), self.reference_velocities, self.prec,
                            err_msg="problem reading .gro file velocities")
        #likewise for each individual atom (to be robust--in case someone alters the individual atom property code):
        assert_almost_equal(all_atoms[0].velocity,self.reference_velocities[0], self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[1].velocity,self.reference_velocities[1], self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[2].velocity,self.reference_velocities[2], self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[3].velocity,self.reference_velocities[3], self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[4].velocity,self.reference_velocities[4], self.prec,
                            err_msg="problem reading .gro file velocities")
        assert_almost_equal(all_atoms[5].velocity,self.reference_velocities[5], self.prec,
                            err_msg="problem reading .gro file velocities")
