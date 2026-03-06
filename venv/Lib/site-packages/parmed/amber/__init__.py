""" Package that provides an API to Amber-specific files """

__author__ = "Jason Swails <jason.swails@gmail.com>"

__all__ = ['AmberFormat', 'FortranFormat', 'AmberAsciiRestart', 'AmberMdcrd',
           'AmberMask', 'NetCDFTraj', 'NetCDFRestart', 'AmberOFFLibrary',
           'AmberParameterSet', 'AmberParm', 'ChamberParm', 'AmoebaParm',
           'Rst7', 'BeemanRestart', 'ConvertFromPSF', 'LoadParm', 'AMBERHOME',
           'titratable_residues']

from .amberformat import AmberFormat, FortranFormat
from .asciicrd import AmberAsciiRestart, AmberMdcrd
from .mask import AmberMask
from .netcdffiles import NetCDFTraj, NetCDFRestart
from .offlib import AmberOFFLibrary
from .parameters import AmberParameterSet
from .readparm import AmberParm, ChamberParm, AmoebaParm, Rst7, BeemanRestart, ConvertFromPSF, LoadParm
from . import titratable_residues

# See if there is an AMBERHOME defined, which we will use by default. Otherwise,
# set it to the empty string
import os as _os
AMBERHOME = _os.getenv('AMBERHOME') or ''
del _os
