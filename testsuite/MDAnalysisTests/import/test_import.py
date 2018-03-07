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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import, print_function
import sys
import os
import warnings
import importlib
import subprocess
import pytest
import mock

"""Test if importing MDAnalysis has unwanted side effects (PR #1794)."""

class TestMDAImport(object):
    # Tests concerning importing MDAnalysis.
    def test_os_dot_fork_restored(self):
        import MDAnalysis
        reload(MDAnalysis)
        try:
            assert os.fork is not None
        except AssertionError as err:
            reload(os)
            raise err

    def test_os_dot_fork_not_called(self):
        """Tests whether os.fork() is called as a side effect when importing
        MDAnalysis. See PR #1794 for details.
        """
        # Modules triggering calls to os.fork() when imported for which
        # there are corresponding workarounds in MDAnalysis:
        modules_with_workarounds = ['uuid']
        # Get a list of all modules which are imported by MDAnalysis:
        mda_modules_command = \
        "from __future__ import absolute_import, print_function\n"
        "import sys\n"
        "modules = set(sys.modules.keys())\n"
        "import MDAnalysis\n"
        "modules = sorted(list(set(sys.modules.keys()) - modules))\n"
        "for module in modules:\n"
        "    print(module)"
        echo = subprocess.Popen(["echo", mda_modules_command],
                                stdout=subprocess.PIPE)
        mda_modules = subprocess.check_output(["python"], stdin=echo.stdout)
        mda_modules = filter(bool, mda_modules.split("\n"))
        # Remove MDAnalysis itself and the modules known to trigger a call to
        # os.fork() from the list of modules imported by MDAnalysis:
        mda_modules = [module for module in mda_modules \
                       if not module.startswith('MDAnalysis') \
                       and module not in modules_with_workarounds]
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ImportWarning)
            # Try reloading each module from the mda_modules list and assert
            # that this doesn't trigger calls to os.fork():
            for module in mda_modules:
                try:
                    with mock.patch('os.fork') as os_dot_fork:
                        try:
                            m = importlib.import_module(module)
                            reload(m)
                        except (ImportError, TypeError, RuntimeError):
                            pass
                        assert not os_dot_fork.called
                except AssertionError as err:
                    print(module)
                    raise err
            # "unload" modules known to trigger a call to os.fork() on import:
            unloaded = {}
            for module in modules_with_workarounds:
                unloaded[module] = sys.modules[module]
                del sys.modules[module]
            # assert that os.fork() isn't called when importing MDAnalysis:
            with mock.patch('os.fork') as os_dot_fork:
                import MDAnalysis
                reload(MDAnalysis)
                assert not os_dot_fork.called
            # restore unloaded modules with originals (prevents breaking other
            # tests):
            for module in modules_with_workarounds:
                sys.modules[module] = unloaded[module]
            del unloaded
