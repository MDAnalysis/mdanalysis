"""
This module contains all of the rism namelist variables for the
amber programs and automatically loads those dictionaries with the   
default values found in that program (sander only).                  

                           GPL LICENSE INFO                             

Copyright (C) 2023 Valdés-Tresanco MS, Valdés-Tresanco ME and Jason Swails

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.
"""
from typing import Dict

from .typing import InputDataType

class rism:

    def __init__(self):
        self.sander: Dict[str, InputDataType] = {
            'closure': ["kh"], 'noasympcorr': 1, 'buffer': 14.0, 'solvcut': -1.0, 'grdspc': [0.5, 0.5, 0.5],
            'ng3': [-1, -1, -1], 'solvbox': [-1, -1, -1], 'tolerance': [1e-05], 'ljtolerance': -1.0,
            'asympkspacetolerance': -1.0, 'treedcf': 1, 'treetcf': 1, 'treecoulomb': 0, 'treedcfmac': 0.1,
            'treetcfmac': 0.1, 'treecoulombmac': 0.1, 'treedcforder': 2, 'treetcforder': 2, 'treecoulomborder': 2,
            'treedcfn0': 500, 'treetcfn0': 500, 'treecoulombn0': 500, 'mdiis_del': 0.7, 'mdiis_nvec': 5,
            'mdiis_restart': 10.0, 'maxstep': 10000, 'npropagate': 5, 'polardecomp': 0, 'entropicdecomp': 0,
            'verbose': 0, 'gfcorrection': 0, 'pcpluscorrection': 0
        }
