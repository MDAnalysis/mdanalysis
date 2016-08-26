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

import numpy as np
from numpy.testing import TestCase
import random

seed = 3856349111 # from: np.uint32(hash("MDAnalysis"))

class RandomState(TestCase):
    """
    Class that sets up random number generators of defined state for each test.

    Use as base from test classes that require defined state random number
    generators.
    
    If an inheriting class overrides :meth:`setUp`, be sure to also call
    :meth:`RandomState.setUp` via :func:`super`.

    At any point in a test each of the RNGs can be reseeded by using the
    respective `seed` method.
    """
    def setUp(self):
        self.np_random = np.random.RandomState(seed)
        self.random = random.Random(seed)

