""" Exceptions used in parmed script """
__all__ = [
    'ParmedError', 'ParmedWarning', 'InputError', 'ParmError', 'ParmWarning', 'SeriousParmWarning',
    'ChangeRadiiError', 'WriteOFFError', 'ParmedUtilsError', 'ParmedChangeError',
    'ParmedAddLJTypeError', 'ChangeLJPairError', 'LJ_TypeError', 'ParmedMoleculeError',
    'ChangeStateError', 'TiMergeError', 'WarningList', 'DeleteDihedralError', 'ArgumentError',
    'NoArgument', 'InterpreterError', 'AmberIncompatibleWarning', 'BadParmWarning', 
    'FixableParmWarning', 'NonfatalWarning', 'NonUniversalWarning', 'MissingDisulfide',
    'LongBondWarning', 'NonexistentParm', 'ParmFileNotFound', 'DuplicateParm', 'AmbiguousParmError',
    'IncompatibleParmsError', 'AddPDBWarning', 'AddPDBError', 'LJ12_6_4Error',
    'DuplicateParamWarning', 'HMassRepartitionError', 'SimulationError', 'SimulationWarning',
    'UnhandledArgumentWarning', 'ParmIndexError', 'FileExists', 'FileDoesNotExist', 'ChamberError',
    'TiMergeError', 'WarningList'
]

from sys import stderr
from ..exceptions import ParmedError, ParmedWarning, InputError
import warnings

class ParmError(ParmedError):
    """ Base parmed error """
    def __init__(self, msg='parmed error'):
        self.msg = msg
    def __str__(self):
        return self.msg

class ParmWarning(ParmedWarning):
    """ Base parmed warning """
    def __init__(self, msg='parmed warning'):
        self.msg = msg
    def __str__(self):
        return self.msg

class SeriousParmWarning(ParmWarning):
    """ These warnings are more serious, and are fatal in strict operation """

# By default, make SeriousParmWarning fatal
warnings.filterwarnings('error', category=SeriousParmWarning)

class ChangeRadiiError(ParmError):
    pass

class WriteOFFError(ParmError):
    pass

class ParmedUtilsError(ParmError):
    pass

class ParmedChangeError(ParmError):
    pass

class ParmedAddLJTypeError(ParmError):
    pass

class ChangeLJPairError(ParmError):
    pass

class LJ_TypeError(ParmError):
    pass

class ParmedMoleculeError(ParmError):
    pass

class ChangeStateError(ParmError):
    pass

class SetParamError(ParmError):
    pass

class DeleteDihedralError(ParmError):
    pass

class ArgumentError(ParmError):
    pass

class NoArgument(ParmError):
    pass

class InterpreterError(ParmError):
    pass

class AmberIncompatibleWarning(ParmWarning):
    pass

class BadParmWarning(ParmWarning):
    pass

class FixableParmWarning(ParmWarning):
    pass

class NonfatalWarning(ParmWarning):
    pass

class NonUniversalWarning(ParmWarning):
    pass

class MissingDisulfide(ParmWarning):
    pass

class LongBondWarning(ParmWarning):
    pass

class NonexistentParm(ParmError):
    pass

class ParmFileNotFound(ParmError, FileNotFoundError):
    pass

class DuplicateParm(ParmError):
    pass

class AmbiguousParmError(ParmError):
    pass

class IncompatibleParmsError(ParmError):
    pass

class AddPDBWarning(ParmWarning):
    pass

class AddPDBError(ParmError):
    pass

class LJ12_6_4Error(ParmError):
    pass

class DuplicateParamWarning(SeriousParmWarning):
    pass

class HMassRepartitionError(ParmError):
    pass

class SimulationError(ParmError):
    pass

class SimulationWarning(ParmWarning):
    pass

class UnhandledArgumentWarning(SeriousParmWarning):
    pass

class ParmIndexError(ParmError, IndexError):
    pass

class FileExists(ParmError):
    pass

class FileDoesNotExist(ParmError, IOError):
    pass

class ChamberError(ParmError):
    pass

class TiMergeError(ParmError):
    pass

class WarningList(list):
    """ List of warnings """
   
    def __init__(self, empty_msg='No warnings found'):
        self._empty_msg = empty_msg
        super().__init__()

    def warn(self, msg, exc_type=ParmWarning):
        """ Adds a warning to the list """
        self.append((exc_type, msg))

    def dump(self, dest=stderr, ncols=80):
        """ Dump a list of all warnings to the destination """
        if len(self) == 0:
            dest.write(self._empty_msg + '\n')
            return

        dest.write(f'{len(self)} total warnings\n\n')

        for w in self:
            words = f'{w[0].__name__}: {w[1]}'.split()
            prstr = words[0] + ' '
            indent_chars = len(words[0]) + 1
            i = 1
            while i < len(words):
                if prstr and len(prstr) + len(words[i]) > ncols:
                    dest.write(prstr + '\n')
                    prstr = ' ' * indent_chars
                prstr += words[i] + ' '
                i += 1
            if prstr:
                dest.write(prstr + '\n')
            dest.write('\n')
