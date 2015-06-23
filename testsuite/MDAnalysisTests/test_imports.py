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

import os
import glob

class TestRelativeImports(object):
    '''Relative imports are banned in unit testing modules (Issue #189), so run tests to enforce this policy.'''

    path_to_testing_modules = os.curdir

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
        list_testing_modules = glob.glob(os.path.join(self.path_to_testing_modules, '*.py'))
        for testing_module in list_testing_modules:
            yield self._run_test_relative_import, testing_module
