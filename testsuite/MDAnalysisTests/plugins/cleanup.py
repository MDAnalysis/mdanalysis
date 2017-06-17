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

"""Plugin that cleans up all .npz created during testing.
"""
# We may want to make the glob pattern an option, for extensibility.
from __future__ import absolute_import

import nose
from nose.plugins.base import Plugin 
from os import walk, unlink, path
from pkg_resources import resource_filename
import warnings

class Cleanup(Plugin):
    """Removes XDR offset files after all testing is done."""
    enabled = True
    name = "mda_cleanup"
    env_opt = 'NOSE_NO_MDA_CLEANUP'
    score = 2000

    def options(self, parser, env):
        """Registers the commandline option, defaulting to enabled.
        """
        parser.add_option("--no-%s" % self.name,
                          action="store_false",
                          dest="do_mda_cleanup",
                          default=not env.get(self.env_opt),
                          help="Disables the cleanup of MDAnalysis offset "
                          "(*.npz) files. [{}]".format(self.env_opt))

    def configure(self, options, conf):
        super(Cleanup, self).configure(options, conf)
        self.config = conf # This will let other tests know about config settings.
        try:
            self.enabled = options.do_mda_cleanup
        except AttributeError:
            self.enabled = False
        self.verbosity = options.verbosity

    def report(self, stream):
        from .. import __name__ as mdatestsname
        dirname = resource_filename(mdatestsname, 'data')
        if self.verbosity > 0:
            stream.write("Cleanup: deleting offset files in "
                "{}\n".format(dirname))
        for root, dirs, fnames in walk(dirname):
            for fname in fnames:
                if fname.endswith('_offsets.npz'):
                    fullname = path.join(root, fname)
                    if self.verbosity > 1:
                        stream.write("Cleanup: deleting offset file "
                            "{}\n".format(fullname))
                    try:
                        unlink(fullname)
                    except OSError:
                        warnings.warn("Cleanup couldn't delete offset file "
                                "{}".format(fullname))


plugin_class = Cleanup

