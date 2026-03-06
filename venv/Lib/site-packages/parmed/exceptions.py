"""
This module contains all of the exceptions that are used in the chemistry
package
"""
import warnings

# Hard errors

class ParmedError(Exception):
    """ Base class for all non-trivial exceptions raised by ParmEd """

class ParsingError(ParmedError):
    """ If there was a problem parsing any kind of file """

class PDBError(ParsingError):
    """ If there was a problem parsing a PDB file """

class Mol2Error(ParsingError):
    """ If there was a problem parsing a Mol2 file """

class MaskError(ParmedError):
    """ Error when a Mask is poorly formed """

class OpenMMError(ParmedError):
    """ If there's a problem making an OpenMM system """

class AmberError(ParsingError):
    """ This is a generic AmberParmError """

class TinkerError(ParsingError):
    """ Raised when one of the TINKER parsing routines hits a bad file """

class CharmmError(ParsingError):
    """ If there is a problem parsing CHARMM PSF files """

class ResidueError(ParmedError):
    """ For when there are problems defining a residue """

class IncompatiblePatchError(ParmedError):
    """ For when applying a PatchTemplate to a ResidueTemplate fails """

class ParameterError(ParmedError):
    """ If a parameter is missing from a database """

class GromacsError(ParmedError):
    """ If there is a problem parsing GROMACS topology files """

class DlpolyError(ParmedError):
    """ If there is a problem parsing DLPOLY topology files """

class FormatNotFound(ParmedError):
    """ If the file format does not have a registered parser with it """

class RosettaError(ParmedError):
    """ If there is a problem loading a Rosetta pose object """

class PreProcessorError(ParmedError):
    """ If there is a problem running the C-like preprocessor """

class MoleculeError(ParmedError):
    """ If there is a problem defining a molecule via the bond graph """

class PdbxError(ParmedError):
    """ Class for catching general errors with PDBx/mmCIF parsing """

class PdbxSyntaxError(PdbxError):
    """ Class for catching errors in mmCIF/PDBx syntax """
    def __init__(self, lineNumber='-1', text=''):
        Exception.__init__(self)
        self.lineNumber = lineNumber
        self.text = text

    def __str__(self):
        return "%%ERROR - [at line: %d] %s" % (self.lineNumber, self.text)

class InputError(ParmedError):
    """ When there is an error with input """

# Warnings

class ParmedWarning(Warning):
    """ Base class for all warnings raised by ParmEd """

# Make sure that all warnings are always printed
warnings.filterwarnings('always', category=ParmedWarning)

class PDBWarning(ParmedWarning):
    """ A non-fatal error to indicate a problematic PDB file """

class AmberWarning(ParmedWarning):
    """ If there is something that is non-fatal """

class SplitResidueWarning(ParmedWarning):
    """ For if a residue with the same number but different names is split """

class CharmmWarning(ParmedWarning):
    """ For non-fatal PSF parsing issues """

class GromacsWarning(ParmedWarning):
    " If we are uncertain about something regarding the GROMACS topology file "

class DlpolyWarning(ParmedWarning):
    " If we are uncertain about something regarding the DLPOLY topology file "

class ParameterWarning(ParmedWarning):
    """ If a type of parameter is missing, but you don't want it to be fatal """

class TinkerWarning(ParmedWarning):
    pass

class PreProcessorWarning(ParmedWarning):
    """ If there is something we should warn about in preprocessing """

class OpenMMWarning(ParmedWarning):
    """ If there is something we should warn when processing OpenMM objects """

# Control flow exceptions

class CharmmPsfEOF(ParmedError):
    """ If we hit an end-of-file in parsing CHARMM files """
