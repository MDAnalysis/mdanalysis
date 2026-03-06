"""
This module contains parsers for sdf file format 
extension described at https://www.cas.org/content/chemical-suppliers/example-sdf
"""
import linecache

from .registry import FileFormatType
from .. import rdkit

class SDFFile(metaclass=FileFormatType):
    """ Class to read SDF file """

    @staticmethod
    def id_format(filename):
        """ Identify the file as a SDF file format or not

        Parameters
        ----------
        filename : str
            Name of the file to test whether or not it is a sdf file

        Returns
        -------
        is_fmt : bool
            True if it is a sdf file, False otherwise
        """
        words = linecache.getline(filename, lineno=4).split()
        return len(words) >= 3 and words[-1] in ('V2000', 'V3000')

    @staticmethod
    def parse(filename, structure=False):
        return rdkit.from_sdf(filename, structure=structure)
