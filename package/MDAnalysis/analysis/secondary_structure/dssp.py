# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2020 The MDAnalysis Development Team and contributors
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

"""Secondary structure analysis --- :mod:`MDAnalysis.analysis.secondary_structure.dssp`
=======================================================================================

:Authors: Lily Wang
:Year: 2020
:Copyright: GNU Public License v3

.. versionadded:: 1.0.0

This module contains a class to compute secondary structure using the DSSP
algorithm [Kabsch1983]_, as implemented by MDTraj_ [McGibbon2015]_ .


Classes
-------

.. autoclass:: DSSP
    :members:

"""

from __future__ import absolute_import, division

import numpy as np

from ...due import due, Doi
from .base import SecondaryStructureBase

due.cite(Doi("10.1016/j.bpj.2015.08.015"),
         description="MDTraj",
         path="MDAnalysis.analysis.secondary_structure.dssp",
         cite_module=True)

due.cite(Doi("10.1002/bip.360221211"),
         description="DSSP",
         path="MDAnalysis.analysis.secondary_structure.dssp",
         cite_module=True)

due.cite(Doi("10.1093/nar/gku1028"),
         description="DSSP (new)",
         path="MDAnalysis.analysis.secondary_structure.dssp",
         cite_module=True)


class DSSP(SecondaryStructureBase):
    """Class for computing secondary structure using the DSSP algorithm,
    as implemented in `MDTraj <http://mdtraj.org/1.9.3/>`. Please cite
    [Kabsch1983]_ , [Touw2015]_ and [McGibbon2015]_ if you use this in 
    published work.


    Parameters
    ----------
    universe: Universe or AtomGroup
        The Universe or AtomGroup to apply the analysis to. As secondary 
        structure is a residue property, the analysis is applied to every 
        residue in your chosen atoms.
    select: string, optional
        The selection string for selecting atoms from ``universe``. The 
        analysis is applied to the residues in this subset of atoms.
    c_name: string, optional
        Atom name for amide carbon in protein backbone. If there are
        multiple possibilities, separate the names with a space, e.g. 
        "C C2 Cz". This gets passed into :meth:`Universe.select_atoms`.
        Names earlier in the string will get preferenced if multiple atoms
        fitting the selection are found in one residue. If multiple atoms
        have the same name, the earlier one (i.e. with lower index) is chosen.
    o_name: string, optional
        Atom name for amide oxygen in protein backbone. Analogous to 
        ``c_name``.
    n_name: string, optional
        Atom name for amide nitrogen in protein backbone. Analogous to 
        ``c_name``.
    ca_name: string, optional
        Atom name for alpha-carbon in protein residues. Analogous to 
        ``c_name``.
    pro_name: string, optional
        Residue name for prolines. Analogous to ``c_name``.
    add_topology_attr: bool, optional
        Whether to add the most common secondary structure as a topology 
        attribute ``secondary_structure`` to your residues. 
    verbose: bool, optional
        Turn on more logging and debugging.


    Attributes
    ----------
    residues: :class:`~MDAnalysis.core.groups.ResidueGroup`
        The residues to which the analysis is applied.
    ss_codes: :class:`numpy.ndarray` of shape (n_frames, n_residues)
        Single-letter secondary structure codes
    ss_names: :class:`numpy.ndarray` of shape (n_frames, n_residues)
        Secondary structure names
    ss_simple: :class:`numpy.ndarray` of shape (n_frames, n_residues)
        Simplified secondary structure names
    ss_counts: dict of {code: :class:`numpy.ndarray` of shape (n_frames,)}
        Dictionary that counts number of residues with each secondary
        structure code, for each frame
    simple_counts: dict of {name: :class:`numpy.ndarray` of shape (n_frames,)}
        Dictionary that counts number of residues with each simplified 
        secondary structure, for each frame
    ss_mode: :class:`numpy.ndarray` of shape (n_residues,)
        The most common secondary structure for each residue
    simple_mode: :class:`numpy.ndarray` of shape (n_residues,)
        The most common simple secondary structure for each residue
    ss_codes_to_names: dict of {code: name}
        Dictionary converting each single-letter code to the full name of
        the secondary structure
    ss_codes_to_simple: dict of {code: name}
        Dictionary converting each single-letter code to simplified 
        secondary structures
    atomgroup: :class:`~MDAnalysis.core.groups.AtomGroup`
        An ordered AtomGroup of backbone atoms
    nco_indices: :class:`numpy.ndarray` of shape (n_residues, 3)
        An array of the atom indices of N, C, O atoms, *relative* to the ``atomgroup``.
        An index is -1 if it is not found in the residue.
    ca_indices: :class:`numpy.ndarray` of shape (n_residues,)
        An array of the atom indices of CA atoms, *relative* to the ``atomgroup``.
        An index is -1 if it is not found in the residue.
    is_proline: :class:`numpy.ndarray` of shape (n_residues,)
        An integer array indicating whether a residue is a proline. 1 is True,
        0 is False.
    not_protein: :class:`numpy.ndarray` of shape (n_residues,)
        A boolean array indicating whether a residue has all of the N, C, O, CA atoms.
    chain_ids: :class:`numpy.ndarray` of shape (n_residues,)
        An integer array representing the index of the chainID or segment each residue 
        belongs to.


    Example
    -------

    ::

        import MDAnalysis as mda
        from MDAnalysis.tests.datafiles import PSF, DCD
        from MDAnalysis.analysis import secondary_structure as ss

        u = mda.Universe(PSF, DCD)
        dssp = ss.DSSP(u, select='backbone', o_name='O O1',
                       add_topology_attr=True)
        dssp.run()
        ss0 = u.residues[0].secondary_structure
        print('The first residue has secondary structure {}'.format(ss0))

    The full secondary structure codes can be accessed at :attr`DSSP.ss_codes`.


    """

    def __init__(self, universe, select='backbone', verbose=False,
                 add_topology_attr=False, c_name='C', o_name='O O1',
                 n_name='N', ca_name='CA', pro_name='PRO'):
        # TODO: implement this on its own w/o mdtraj?
        try:
            from mdtraj.geometry._geometry import _dssp
        except ImportError:
            raise ValueError('DSSP requires mdtraj to be installed')
        else:
            self._mdtraj_dssp = _dssp

        super(DSSP, self).__init__(universe, select=select, verbose=verbose,
                                   add_topology_attr=add_topology_attr)
        self.names = {'C': c_name,
                      'O': o_name,
                      'N': n_name,
                      'CA': ca_name,
                      'PRO': pro_name}

        self._setup_atomgroups()

    def _setup_atomgroups(self):
        """
        Set up atomgroup and indices for MDTraj.
        """
        ag = self.residues.atoms
        ridx = self.residues.resindices

        # grab n, c, o
        ncoca = ag[[]]
        ncoca_indices = np.array([[-1]*self.n_residues]*4, dtype=np.int32)
        for i, atom in enumerate(['N', 'C', 'O', 'CA']):
            name = self.names[atom]
            selections = ['name {}'.format(n) for n in name.split()]
            atoms = ag.select_atoms(*selections)

            # index of atoms, within atomgroup
            residx, idx, counts = np.unique(atoms.resindices,
                                            return_index=True,
                                            return_counts=True)

            if np.count_nonzero(counts > 1):
                # some residues have atoms with the same name
                # take the first one? That's what MDTraj does
                atoms = atoms[idx]

            # -1 means missing
            found_idx = np.isin(ridx, residx, assume_unique=True)
            ncoca_indices[i][found_idx] = np.arange(len(atoms)) + len(ncoca)
            ncoca += atoms

        self.atomgroup = ncoca
        self.nco_indices = np.ascontiguousarray(ncoca_indices[:3].T)
        self.ca_indices = np.ascontiguousarray(ncoca_indices[-1])
        is_protein = (ncoca_indices != -1).sum(axis=0) == 4
        self.not_protein = ~is_protein

        pro = ag.select_atoms('resname {}'.format(self.names['PRO']))
        pro_idx = pro.residues.resindices
        is_proline = np.isin(ridx, pro_idx, assume_unique=True)
        self.is_proline_indices = is_proline.astype(np.int32)

        # chains need to be ints for mdtraj
        if hasattr(self._universe._topology, 'chainIDs'):
            chains = [r.atoms.chainIDs[0] for r in self.residues]
        else:
            chains = self.residues.segindices
        cids, cidx = np.unique(chains, return_inverse=True)
        self.chain_ids = np.arange(len(cids), dtype=np.int32)[cidx]

    def _compute_dssp(self):
        # mdtraj uses angstrom
        xyz = self.atomgroup.positions.reshape((1, -1, 3))/10
        xyz = np.ascontiguousarray(xyz.astype(np.float32))
        value = self._mdtraj_dssp(xyz, self.nco_indices, self.ca_indices,
                                  self.is_proline_indices,
                                  self.chain_ids)
        arr = np.fromiter(value, dtype='<U1')
        arr[arr == ' '] = 'C'  # replace random coil
        arr[self.not_protein] = ''  # replace non-protein
        self.ss_codes[self._frame_index] = arr
