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
from __future__ import print_function

from six.moves import range

import MDAnalysis
import MDAnalysis.analysis.hole
from MDAnalysis.analysis.hole import HOLEtraj

from numpy.testing import TestCase, dec
import numpy as np
import nose
from nose.plugins.attrib import attr

import errno

from MDAnalysisTests.datafiles import PDB_HOLE, XTC_HOLE
from MDAnalysisTests import executable_not_found, tempdir

def rlimits_missing():
    # return True if resources module not accesible (ie setting of rlimits)
    try:
        # on Unix we can manipulate our limits: http://docs.python.org/2/library/resource.html
        import resource

        soft_max_open_files, hard_max_open_files = resource.getrlimit(resource.RLIMIT_NOFILE)
    except ImportError:
        return True
    return False

class TestHoleModule(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(PDB_HOLE, XTC_HOLE)
        try:
            # on Unix we can manipulate our limits: http://docs.python.org/2/library/resource.html
            import resource
            self.soft_max_open_files, self.hard_max_open_files = resource.getrlimit(resource.RLIMIT_NOFILE)
        except ImportError:
            pass

    @attr('slow')
    @attr('issue')
    @dec.skipif(rlimits_missing, msg="Test skipped because platform does not allow setting rlimits")
    @dec.skipif(executable_not_found("hole"), msg="Test skipped because HOLE not found")
    def test_hole_module_fd_closure(self):
        """MDAnalysis.analysis.hole: Issue 129: ensure low level file descriptors to PDB files used by Hole program are properly closed"""
        # If Issue 129 isn't resolved, this function will produce an OSError on
        # the system, and cause many other tests to fail as well.
        #
        # Successful test takes ~10 s, failure ~2 s.
        try:
            # Hasten failure by setting "ulimit -n 64" (can't go too low because of open modules etc...)
            import resource
            resource.setrlimit(resource.RLIMIT_NOFILE, (64, self.hard_max_open_files))
        except ImportError:
            raise NotImplementedError("Test cannot be run without the resource module.")

        with tempdir.in_tempdir():
            try:
                # will need to have the 'hole' command available in the path
                H = HOLEtraj(self.universe, cvect=[0, 1, 0], sample=20.0)
            except OSError as err:
                if err.errno == errno.ENOENT:
                    raise OSError(errno.ENOENT, "HOLE binary not found")
                raise
            finally:
                self._restore_rlimits()

            # pretty unlikely that the code will get through 2 rounds if the MDA
            # issue 129 isn't fixed, although this depends on the file descriptor
            # open limit for the machine in question
            try:
                for i in range(2):
                    # will typically get an OSError for too many files being open after
                    # about 2 seconds if issue 129 isn't resolved
                    H.run()
            except OSError as err:
                if err.errno == errno.EMFILE:
                    raise AssertionError("HOLEtraj does not close file descriptors (Issue 129)")
                elif err.errno == errno.ENOENT:
                    raise OSError(errno.ENOENT, "HOLE binary not found")
                raise
            finally:
                # make sure to restore open file limit !!
                self._restore_rlimits()

    def _restore_rlimits(self):
        try:
            import resource

            resource.setrlimit(resource.RLIMIT_NOFILE, (self.soft_max_open_files, self.hard_max_open_files))
        except ImportError:
            pass

    def tearDown(self):
        self._restore_rlimits()
        del self.universe
