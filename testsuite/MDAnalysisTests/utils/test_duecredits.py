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
from __future__ import absolute_import
import logging
import pytest
import os,sys
os.environ['DUECREDIT_ENABLE']='yes'

logger = logging.getLogger("MDAnalysisTests.utils")
#Testing the import#
try:
    import MDAnalysis.due
except ImportError:
    logger.info('Could not find due.py,all tests will be skipped')
    pass


class TestDuecredits():
        
    def test_duecredit_active(self):
        import MDAnalysis as mda
        assert mda.due.active == True

    
    def test_duecredit_collector_citations(self):
        import MDAnalysis as mda
        assert mda.due.citations[('MDAnalysis/',
                                  'gowers2016')].cites_module == True        
        assert mda.due.citations[('MDAnalysis/',
                                  '10.1002/jcc.21787')].cites_module == True
        
    


