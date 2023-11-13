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


"""GROMOS11 trajectory reader --- :mod:`MDAnalysis.coordinates.TRC`
====================================================================

Reads coordinates, timesteps and box-sizes from GROMOS11 TRC trajectories.
The trajectory format is documented in the GROMOS Manual Vol. 4, chapter 2 and 4.
The manual can be downloaded here: 
https://gromos.net/gromos11_pdf_manuals/vol4.pdf
The code has been tested with GROMOS11 version 1.6 (Dec. 2023)

To load the trajectory into :class:`~MDAnalysis.core.universe.Universe`,
you need to provide topology information using a pdb::

    import MDAnalysis as mda
    u = mda.Universe("topology.pdb", ["md_1.trc.gz","md_2.trc.gz"], 
    continuous=True)

.. Note::
   The GROMOS trajectory format organizes its data in blocks. A block starts 
   with a blockname in capital letters and ends with a line containing only ''END''.
   Only the TITLE-block at the beginning of each file is mandatory, 
   others blocks can be chosen depending on the task. 

   This reader is designed to read the blocks "TIMESTEP", "POSITIONRED" and 
   "GENBOX" from the trajectory which covers most standard trajectories.
   
   MDAnalysis requires the blocks of each frame to be in the same order 
   and ignores non-supported blocks.
   
Classes
-------

.. autoclass:: TRCReader
   :members:

"""

import os
import errno
import warnings
import numpy as np

from . import base
from .timestep import Timestep
from ..lib import util
from ..lib.util import cached, store_init_arguments
from ..exceptions import NoDataError
from ..version import __version__

import logging
logger = logging.getLogger("MDAnalysis.coordinates.GROMOS11")


class TRCReader(base.ReaderBase):
    """Coordinate reader for the GROMOS11 format
    """

    format = 'TRC'
    units = {'time': 'ps', 'length': 'nm'}
    _Timestep = Timestep
    
    @store_init_arguments
    def __init__(self, filename, **kwargs):
        super(TRCReader, self).__init__(filename, **kwargs)
        
        self.SUPPORTED_BLOCKS = ['TITLE', 'TIMESTEP', 'POSITIONRED', 
        'GENBOX']
        self.NOT_SUPPORTED_BLOCKNAMES = ['POSITION', 'REFPOSITION', 
        'VELOCITY', 'VELOCITYRED', 'FREEFORCE', 'FREEFORCERED', 
        'CONSFORCE', 'CONSFORCERED']
                
        # GROMOS11 trajectories can be either *.trc or *.trc.gz.
        root, ext = os.path.splitext(self.filename)
        self.trcfile = util.anyopen(self.filename)
        self.compression = ext[1:] if ext[1:] != "trc" else None

        # Read and calculate some information about the trajectory
        self.traj_properties = self._read_traj_properties()
                
        self._cache = {}
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        
        self._reopen()
        self.ts.dt = self.traj_properties["dt"]
        
        self._read_frame(0)
        
    @property
    @cached('n_atoms')
    def n_atoms(self):
        """The number of atoms in one frame."""
        try:
            return self._read_atom_count()
        except IOError:
            return 0

    def _read_atom_count(self):
        n_atoms = self.traj_properties["n_atoms"] 
        return n_atoms
        
    @property
    @cached('n_frames')
    def n_frames(self):
        """The number of frames in the trajectory."""
        try:
            return self._read_frame_count()
        except IOError:
            return 0

    def _read_frame_count(self):
        n_frames = self.traj_properties["n_frames"] 
        return n_frames
                
    def _frame_to_ts(self, frameDat, ts):
        """Convert a frame to a :class: TimeStep"""
        ts.frame = self._frame
        ts.time = frameDat["time"]
        
        ts.data['time'] = frameDat["time"]
        ts.data['step'] = frameDat["step"]
        
        ts.dimensions = frameDat["dimensions"]
        ts.positions = frameDat["positions"]

        #
        # Convert the units
        #        
        if self.convert_units:
            if ts.has_positions:
                self.convert_pos_from_native(ts.positions)
            if ts.dimensions is not None:
                self.convert_pos_from_native(ts.dimensions[:3])
            if ts.has_velocities:
                self.convert_velocities_from_native(ts.velocities)
      
        return ts
    
    def _read_traj_properties(self):
        """
        * Reads the number of atoms per frame (n_atoms)
        * Reads the number of frames (n_frames) 
        * Reads the startposition of the timestep block 
          for each frame (l_blockstart_offset) 
        """
        
        traj_properties = {}
        
        #
        # Check which of the supported blocks comes first win the trajectory
        #
        first_block = None
        with util.anyopen(self.filename) as f:   
            for line in iter(f.readline,''):
                for blockname in self.SUPPORTED_BLOCKS:
                    if (blockname in line) and (blockname != 'TITLE'):
                        #Save the name of the first non-Title block 
                        #in the trajectory file
                        first_block = blockname
                        
                if (first_block is not None): break #First block found
        
        #
        # Calculate meta-data of the trajectory
        #                               
        in_positionred_block = False
        lastline_was_timestep = False
        
        atom_counter = 0
        n_atoms = 0
        frame_counter = 0  
 
        l_blockstart_offset = []
        l_timestep_timevalues = []
        
        with util.anyopen(self.filename) as f:
            for line in iter(f.readline,''):

                #
                # First block of frame
                #
                if (first_block in line):
                    l_blockstart_offset.append(f.tell()-len(line))    
                    frame_counter += 1
                
                # 
                # Timestep-Block
                # 
                if "TIMESTEP" in line:
                    lastline_was_timestep = True
                     
                elif (lastline_was_timestep == True):
                    l_timestep_timevalues.append(float(line.split()[1]))
                    lastline_was_timestep = False       

                #
                # Coordinates-Block
                #
                if "POSITIONRED" in line:
                    in_positionred_block = True
 
                if (in_positionred_block == True) and (n_atoms == 0):
                    if (len(line.split()) == 3):
                        atom_counter += 1  
                
                if ("END" in line) and (in_positionred_block == True):
                    n_atoms = atom_counter
                    in_positionred_block = False
                     

        traj_properties["n_atoms"] = n_atoms
        traj_properties["n_frames"] = frame_counter
        traj_properties["l_blockstart_offset"] = l_blockstart_offset
        
        if (len(l_timestep_timevalues) >= 2):
            traj_properties["dt"] = (l_timestep_timevalues[1] 
                                     - l_timestep_timevalues[0])
        else:
            traj_properties["dt"] = 0
            warnings.warn("The trajectory does not contain TIMESTEP \
            information!", UserWarning)
            
        return traj_properties

    def _read_GROMOS11_trajectory(self):
    
        frameDat = {}
        f = self.trcfile
        if (f.closed):
            raise Exception("The trajectory has been closed before reading.")

        # Read trajectory    
        for line in iter(f.readline,''):

            if ("TIMESTEP" in line): 
                tmp_step, tmp_time = f.readline().split()
                frameDat["step"] = int(tmp_step)
                frameDat["time"] = float(tmp_time)
                
            elif ("POSITIONRED" in line): 
                tmp_buf = []
                while (True): 
                    coords_str = f.readline()
                    if '#' in coords_str:
                        continue
                    elif "END" in coords_str:
                        break
                    else:
                        tmp_buf.append(coords_str.split())
                    
                if (np.array(tmp_buf).shape[0] == self.n_atoms):
                    frameDat["positions"] = np.asarray(tmp_buf, 
                                                       dtype=np.float64)
                else:
                    raise ValueError("The trajectory contains \
                                     the wrong number of atoms!")                 

            elif ("GENBOX" in line): 
                ntb_setting = int(f.readline())
                if (ntb_setting == 0):
                    frameDat["dimensions"] = None
                    self.periodic = False
                                    
                elif (ntb_setting in [-1, 1]):  
                    tmp_a, tmp_b, tmp_c = f.readline().split()
                    tmp_alpha, tmp_beta, tmp_gamma = f.readline().split()
                    frameDat["dimensions"] = [float(tmp_a), 
                                              float(tmp_b), 
                                              float(tmp_c), 
                                              float(tmp_alpha),
                                              float(tmp_beta), 
                                              float(tmp_gamma)]
                    self.periodic = True
                                
                    line3 = f.readline().split()
                    line4 = f.readline().split()          
                    for v in (line3 + line4):
                        if (float(v) != 0.0):
                            raise NotImplementedError("This reader \
                                                      supports neither \
                                                      triclinic and/or \
                                                      (yawed,pitched, \
                                                      rolled) boxes!")
                else:
                    raise NotImplementedError("This reader does only support\
                                               vacuum and rectangular boxes!")                             
                break
                
            elif any(non_supp_bn in line for non_supp_bn in self.NOT_SUPPORTED_BLOCKNAMES):
                for non_supp_bn in self.NOT_SUPPORTED_BLOCKNAMES:
                    if (non_supp_bn in line):
                        warnings.warn("Block "+non_supp_bn+" is not supported!", UserWarning)

                        
        return frameDat

    
    def _read_frame(self, i):
        """read frame i"""
        self._frame = i - 1

        # Move position in file just (-2 byte) before the start of the block 
        self.trcfile.seek(self.traj_properties["l_blockstart_offset"][i]-2, 0)

        return self._read_next_timestep()

    def _read_next_timestep(self, ts=None):
        if ts is None:
            ts = self.ts
        
        self._frame += 1
        if (self._frame >= self.n_frames):
            raise EOFError('Trying to go over trajectory limit')

        raw_framedata = self._read_GROMOS11_trajectory()        
        self._frame_to_ts(raw_framedata, ts)
        self.ts = ts
        
        return ts
        
    def _reopen(self):
        """Close and reopen the trajectory"""
        self.close()
        self.open_trajectory()

    def open_trajectory(self):
        if self.trcfile is not None:
            raise IOError(
                errno.EALREADY, 'TRC file already opened', self.filename)

        # Reload trajectory file
        self.trcfile = util.anyopen(self.filename)

        # Reset ts
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        
        # Set frame to -1, so next timestep is zero
        self._frame = -1 

        return self.trcfile

    def close(self):
        """Close the trc trajectory file if it was open."""
        if self.trcfile is None:
            return
        self.trcfile.close()
        self.trcfile = None
        
