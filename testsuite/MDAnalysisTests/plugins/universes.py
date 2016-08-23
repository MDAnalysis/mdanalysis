# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning,
# Oliver Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

from __future__ import print_function

import collections
import os
import random
import sys
import weakref

import gc
import objgraph

from nose.plugins.base import Plugin
from nose.pyversion import exc_to_unicode, force_unicode
from nose.util import ln


class ReportUniverse(Plugin):
    """Track leaking Universes.
    """
    enabled = False
    env_opt = 'NOSE_MDA_UNIVERSE'
    name = 'universe'

    def options(self, parser, env):
        """Registers the commandline option, defaulting to disabled.
        """
        parser.add_option(
            "--with-{}".format(self.name), action="store_true",
            default=env.get(self.env_opt, False), dest="do_mda_universe",
            help="Track leaking Universes. [{}]".format(self.env_opt))

    def configure(self, options, conf):
        """Configure the plugin.
        """
        super(ReportUniverse, self).configure(options, conf)
        self.config = conf # This will let other tests know about config settings.
        try:
            self.enabled = options.do_mda_universe
        except AttributeError:
            self.enabled = False

    def report(self, stream):
        # We need to import MDAnalysis to have the list of the living Universes.
        # Yet, we should avoid importing the library at the beginning of the
        # test as it could mask some errors and mess up with coverage tracking.
        # Therefore, MDAnalysis is imported here--in the function--only were
        # it is needed, and only if the plugin is enabled.
        import MDAnalysis as mda

        print('\n', file=stream)
        print('There are {} living Universes at the end of the tests.'
              .format(len(mda._anchor_universes)), file=stream)
        if mda._anchor_universes:
            print('These univerese are still referenced:', file=stream)
            for universe in mda._anchor_universes:
                reader = universe.trajectory
                referrers = gc.get_referrers(universe)
                ref_counter = collections.Counter(r.__class__ for r in referrers)
                weak_referrers = weakref.getweakrefs(universe)
                # Exclude the reference we created in the loop.
                n_references = sys.getrefcount(universe) - 1
                n_weakrefs = weakref.getweakrefcount(universe)
                print('* {}'.format(universe), file=stream)
                print('  reader: {}'.format(reader), file=stream)
                print('  {} references to the universe'.format(n_references),
                      file=stream)
                print('  referenced by ' + str(ref_counter), file=stream)
                print('  {} weak references to the universe'.format(n_weakrefs),
                      file=stream)
                print('  weak references by ' + str(weak_referrers), file=stream)
        
        living_atoms = [o for o in gc.get_objects() if isinstance(o, mda.core.AtomGroup.Atom)]
        print('There are {} living atoms.'.format(len(living_atoms)), file=stream)

        objgraph.show_refs(living_atoms[0], filename='sample-graph.png')
        objgraph.show_backrefs(living_atoms[0], filename='backrefs.png', max_depth=10)


plugin_class = ReportUniverse
