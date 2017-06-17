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

from __future__ import print_function, absolute_import

import collections
import os

import psutil

from nose.plugins.base import Plugin
from nose.pyversion import exc_to_unicode, force_unicode
from nose.util import ln


class ReportOpenFiles(Plugin):
    """Follow open files during the test suite.

    Open files are listed when a test fails or errors. They also are listed at
    the end of the test suite.

    The plugin is enable by default. It can be disabled with the
    --no-open-files option or with the NOSE_NO_MDA_FILES environment variable.
    """
    enabled = True
    env_opt = 'NOSE_NO_MDA_FILES'
    name = 'open-files'

    def options(self, parser, env):
        """Registers the commandline option, defaulting to enabled.
        """
        parser.add_option(
            "--no-{}".format(self.name), action="store_false",
            default=not env.get(self.env_opt), dest="do_mda_open_files",
            help="Don't display open file handlers. [{}]".format(self.env_opt))

    def configure(self, options, conf):
        """Configure the plugin.
        """
        super(ReportOpenFiles, self).configure(options, conf)
        self.config = conf # This will let other tests know about config settings.
        try:
            self.enabled = options.do_mda_open_files
        except AttributeError:
            self.enabled = False

    def formatError(self, test, err):
        """List the open files when a test errors.
        """
        open_files = _gather_open_files()

        if not open_files:
            return err

        ec, ev, tb = err
        handle_couter = collections.Counter(open_files)
        new_ev = u'\n'.join(
            [ev, ln(u'>>> There {} open file handlers for {} files: <<<'
                  .format(len(open_files), len(handle_couter)))]
            + [ln(u'* ({}) {}'.format(count, path))
               for path, count in handle_couter.items()]
            + [ln(u'>>> End of the open file list. <<<')]
        )
        return (ec, new_ev, tb)

    def formatFailure(self, test, err):
        """List the open files when a test fails.
        """
        return self.formatError(test, err)

    def report(self, stream):
        """Display the number of open files at the end of the test suite.
        """
        open_files = _gather_open_files()
        handle_couter = collections.Counter(open_files)
        print('\n', file=stream)
        print('By the end of the tests, there are {} open handles for {} files:'
              .format(len(open_files), len(handle_couter)), file=stream)
        for path, count in handle_couter.items():
            print('* ({}) {}'.format(count, path), file=stream)


def _gather_open_files():
    """Return a list of the path to open files handled by the process.

    If several file handle are open for a given file, then the path will
    appear as many time as they are open file handle for that file.
    """
    process = psutil.Process(os.getpid())
    handlers = process.open_files()
    path = [f.path for f in handlers]
    return path


plugin_class = ReportOpenFiles
