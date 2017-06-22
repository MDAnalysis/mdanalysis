# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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
from __future__ import absolute_import

import os
import MDAnalysisTests

class TestRelativeImports(object):
    """Relative imports are banned in unit testing modules (Issue #189), so run tests to enforce this policy."""

    path_to_testing_modules = MDAnalysisTests.__path__[0]
    # Exclusion path relative to MDAnalysisTests
    exclusions = ['/plugins']

    @classmethod
    def is_excluded(cls, path):
        leaf = path[len(cls.path_to_testing_modules):]
        return leaf in cls.exclusions

    @staticmethod
    def _run_test_relative_import(testing_module):
        with open(testing_module, 'r') as test_module_file_object:
            for lineno, line in enumerate(test_module_file_object, start=1):
                if 'from .' in line and 'import' in line \
                        and not 'test_imports' in testing_module:
                    raise AssertionError(
                        "A relative import statement was found in "
                        "module {testing_module} at linenumber {lineno}.".format(**vars()))

    def test_relative_imports(self):
        for dirpath, dirnames, files in os.walk(self.path_to_testing_modules):
            if self.is_excluded(dirpath):
                continue
            for f in filter(lambda x: x.endswith('.py'), files):
                fpath = os.path.join(dirpath, f)
                if self.is_excluded(fpath):
                    continue
                yield self._run_test_relative_import, fpath
