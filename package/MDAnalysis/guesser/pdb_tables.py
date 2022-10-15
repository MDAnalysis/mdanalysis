# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

from .tables import kv2dict

amino_acids = ["ACE", "ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "FOR", "GLN", "GLU", "GLX", "GLY", "HIS",
               "HYP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "PCA", "SER", "THR", "TRP", "TYR", "UNK", "VAL"]

nucleic_acids = [" A", " C", " G", " T", " U", " +U", " YG", "1MA",
                 "1MG", "2MG", "5MC", "5MU", "7MG", "H2U", "M2G", "OMC", "OMG", "PSU"]

known_groups = ["101", "12A", "1AR", "1GL", "2AS", "2GL", "3AA", "3AT", "3DR", "3PO", "6HA", "6HC", "6HG", "6HT", "A26", "AA6", "ABD", "AC1", "ACO", "AIR", "AMU", "AMX", "AP5", "AMG", "APU", "B9A", "BCA", "BNA", "CAA",
                "CBS", "CGS", "CMC", "CND", "CO8", "COA", "COF", "COS", "DCA", "DGD", "FAB", "FAD", "FAG", "FAM", "FDA", "GPC", "IB2", "NAD", "NAH", "NAI", "NAL", "NAP", "NBD", "NDP", "PAD", "SAD", "SAE", "T5A", "tRE", "UP5","ZID"]

standard_groups = amino_acids + nucleic_acids + known_groups
