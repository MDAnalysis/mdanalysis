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

"""Secondary structure analysis --- :mod:`MDAnalysis.analysis.secondary_structure.wrapper`
==========================================================================================

:Authors: Lily Wang
:Year: 2020
:Copyright: GNU Public License v3

.. versionadded:: 1.0.0

This module contains classes for computing secondary structure using external programs. 
MDAnalysis provides a wrapper for the  DSSP_ and STRIDE_ programs to be run on an MD trajectory.

Both DSSP_ and STRIDE_ need to be installed separately.

.. _DSSP: https://swift.cmbi.umcn.nl/gv/dssp/
.. _STRIDE: http://webclu.bio.wzw.tum.de/stride/


Classes and functions
---------------------

.. autoclass:: DSSPWrapper
.. autoclass:: STRIDEWrapper

"""

import os
import errno
import tempfile
import subprocess
import warnings
import numpy as np

from ...due import due, Doi
from ...lib import util
from .base import SecondaryStructureBase

due.cite(Doi("10.1002/bip.360221211"),
         description="DSSP",
         path="MDAnalysis.analysis.secondary_structure.wrapper",
         cite_module=True)

due.cite(Doi("10.1093/nar/gku1028"),
         description="DSSP (new)",
         path="MDAnalysis.analysis.secondary_structure.wrapper",
         cite_module=True)

due.cite(Doi("10.1002/prot.340230412"),
         description="STRIDE",
         path="MDAnalysis.analysis.secondary_structure.wrapper",
         cite_module=True)


class SecondaryStructureWrapper(SecondaryStructureBase):
    """Base class for secondary structure analysis that wraps an 
    external program.
    """
    cmd = ''

    @property
    def exe_name(self):
        raise NotImplementedError

    def __init__(self, executable, universe, select='protein', verbose=False,
                 add_topology_attr=False):

        exe = util.which(executable)
        if exe is None:
            msg = ('{name} executable not found at {executable}. '
                   '{exe_name} must be on the PATH, or the path must '
                   'be provided with the keyword argument executable')
            raise OSError(errno.ENOENT, msg.format(name=type(self).__name__,
                                                   executable=executable,
                                                   exe_name=self.exe_name))
        self.exe = exe
        super(SecondaryStructureWrapper, self).__init__(universe,
                                                        select=select,
                                                        verbose=verbose,
                                                        add_topology_attr=add_topology_attr)

    def _execute(self):
        """Run the wrapped program."""
        # ignore PDB warnings
        warnings.filterwarnings("ignore",
                                message="Found no information for attr")

        fd, pdbfile = tempfile.mkstemp(suffix='.pdb')
        os.close(fd)
        try:
            self.atomgroup.write(pdbfile)
            cmd_args = [self.exe] + self.cmd.format(pdb=pdbfile).split()
            proc = subprocess.Popen(cmd_args, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
        finally:
            try:
                os.unlink(pdbfile)
            except OSError:
                pass
        return stdout.decode('utf-8')

    def _prepare(self):
        super(SecondaryStructureWrapper, self)._prepare()
        nf = self.n_frames
        self.phi = np.zeros((nf, self.n_residues), dtype=float)
        self.psi = np.zeros((nf, self.n_residues), dtype=float)
        self.sasa = np.zeros((nf, self.n_residues), dtype=float)

    def _compute_dssp(self):
        output = self._execute()
        self._process_output(output)


class DSSPWrapper(SecondaryStructureWrapper):
    """Runs :program:`mkdssp` on a trajectory.

    :program:`mkdssp` implements the DSSP algorithm to determine secondary
    structure. Please cite [Kabsch1983]_ and [Touw2015]_ if you use this in 
    published work.

    This class creates temporary PDB files for each frame and runs ``mkdssp`` on them.

    Parameters
    ----------
    universe: Universe or AtomGroup
        The Universe or AtomGroup to apply the analysis to. As secondary 
        structure is a residue property, the analysis is applied to every 
        residue in your chosen atoms.
    select: string, optional
        The selection string for selecting atoms from ``universe``. The 
        analysis is applied to the residues in this subset of atoms.
    executable: str, optional
        Path to the ``mkdssp`` executable.
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
    phi: :class:`numpy.ndarray` of shape (n_frames, n_residues)
        Phi angles
    psi: :class:`numpy.ndarray` of shape (n_frames, n_residues)
        Psi angles
    sasa: :class:`numpy.ndarray` of shape (n_frames, n_residues)
        Solvent-accessible surface area
    ss_counts: dict of {code: :class:`numpy.ndarray` of shape (n_frames,)}
        Dictionary that counts number of residues with each secondary
        structure code, for each frame
    simple_counts: dict of {name: :class:`numpy.ndarray` of shape (n_frames,)}
        Dictionary that counts number of residues with each simplified 
        secondary structure, for each frame
    ss_mode: :class:`numpy.ndarray` of shape (n_residues,)
        The most common secondary structure for each residue
    ss_codes_to_names: dict of {code: name}
        Dictionary converting each single-letter code to the full name of
        the secondary structure
    ss_codes_to_simple: dict of {code: name}
        Dictionary converting each single-letter code to simplified 
        secondary structures
    """

    exe_name = 'mkdssp'
    cmd = '-i {pdb}'

    columns = {
        'ss': 16,
        'sasa': (35, 39),
        'phi': (104, 109),
        'psi': (111, 116),
    }

    def __init__(self, universe, select='protein', executable='mkdssp',
                 verbose=False, add_topology_attr=False):
        super(DSSPWrapper, self).__init__(executable, universe,
                                          select=select, verbose=verbose,
                                          add_topology_attr=add_topology_attr)

    def _process_output(self, output):
        frame = self._frame_index
        per_res = output.split('  #  RESIDUE AA STRUCTURE BP1 BP2  ACC')[1]
        i = 0
        for line in per_res.split('\n')[1:]:
            try:
                if line[13] == '!':
                    continue
            except IndexError:
                continue
            ss = line[self.columns['ss']]

            if ss == ' ':
                ss = 'C'
            self.ss_codes[frame][i] = ss
            for kw in ('phi', 'psi', 'sasa'):
                _i, _j = self.columns[kw]
                getattr(self, kw)[frame][i] = line[_i:_j]
            i += 1


class STRIDEWrapper(SecondaryStructureWrapper):
    """Runs :program:`stride` on a trajectory.

    :program:`stride` implements the STRIDE algorithm to determine secondary
    structure. Please cite [Heinig2004]_ if you use this in published work.

    This class creates temporary PDB files for each frame and runs ``stride`` on them.

    Parameters
    ----------
    universe: Universe or AtomGroup
        The Universe or AtomGroup to apply the analysis to. As secondary 
        structure is a residue property, the analysis is applied to every 
        residue in your chosen atoms.
    select: string, optional
        The selection string for selecting atoms from ``universe``. The 
        analysis is applied to the residues in this subset of atoms.
    executable: str, optional
        Path to the ``stride`` executable.
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
    phi: :class:`numpy.ndarray` of shape (n_frames, n_residues)
        Phi angles
    psi: :class:`numpy.ndarray` of shape (n_frames, n_residues)
        Psi angles
    sasa: :class:`numpy.ndarray` of shape (n_frames, n_residues)
        Solvent-accessible surface area
    ss_counts: dict of {code: :class:`numpy.ndarray` of shape (n_frames,)}
        Dictionary that counts number of residues with each secondary
        structure code, for each frame
    simple_counts: dict of {name: :class:`numpy.ndarray` of shape (n_frames,)}
        Dictionary that counts number of residues with each simplified 
        secondary structure, for each frame
    ss_mode: :class:`numpy.ndarray` of shape (n_residues,)
        The most common secondary structure for each residue
    ss_codes_to_names: dict of {code: name}
        Dictionary converting each single-letter code to the full name of
        the secondary structure
    ss_codes_to_simple: dict of {code: name}
        Dictionary converting each single-letter code to simplified 
        secondary structures
    """

    exe_name = 'stride'
    cmd = '{pdb}'

    def __init__(self, universe, select='protein', executable='stride',
                 verbose=False, add_topology_attr=False):
        super(STRIDEWrapper, self).__init__(executable, universe,
                                            select=select, verbose=verbose,
                                            add_topology_attr=add_topology_attr)

    def _process_output(self, output):
        lines = output.split('\nASG ')[1:]
        lines[-1] = lines[-1].split('\n')[0]
        frame = self._frame_index
        for i, line in enumerate(lines):
            # resname chain resid resnum ss ssname phi psi area ~~~~
            fields = line.split()
            self.ss_codes[frame][i] = fields[4]
            self.phi[frame][i] = fields[6]
            self.psi[frame][i] = fields[7]
            self.sasa[frame][i] = fields[8]
