"""
This module contains all of the gbnsr6 namelist variables for the
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
class gbnsr6:

    def __init__(self):
        self.gbnsr6: Dict[str, InputDataType] = {
            'b': 0.028, 'alpb': 1, 'epsin': 1.0, 'epsout': 78.5, 'istrng': 0.0, 'rs': 0.52, 'dprob': 1.4,
            'space': 0.5, 'arcres': 0.2, 'rbornstat': 0, 'dgij': 0, 'radiopt': 0, 'chagb': 0, 'roh': 0.586,
            'tau': 1.47, 'cavity_surften': 0.005
        }
