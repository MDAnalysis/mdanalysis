# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""
Useful functions for running tests

"""

import __builtin__
from functools import wraps
import mock


def block_import(package):
    """Block import of a given package

    eg:

    @blocker('numpy')
    def try_and_do_something():
        import numpy as np  # this will fail!

    Shadows the builtin import method, sniffs import requests
    and blocks the designated package.
    """
    def blocker_wrapper(func):
        @wraps(func)
        def func_wrapper(*args, **kwargs):
            with mock.patch('__builtin__.__import__',
                            wraps=__builtin__.__import__) as mbi:
                def blocker(*args, **kwargs):
                    if package in args:
                        raise ImportError("Blocked by block_import")
                    else:
                        # returning DEFAULT allows the real function to continue
                        return mock.DEFAULT
                mbi.side_effect = blocker
                func(*args, **kwargs)
        return func_wrapper
    return blocker_wrapper


