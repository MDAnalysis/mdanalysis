""" A list of helpful decorators for use throughout ParmEd """
from functools import wraps
import warnings

try:
    import openmm as mm
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False

__all__ = ['needs_openmm']

def needs_openmm(fcn):
    @wraps(fcn)
    def new_fcn(*args, **kwargs):
        if not HAS_OPENMM:
            raise ImportError('Could not find or import OpenMM version 7.0+')
        return fcn(*args, **kwargs)

    return new_fcn

def deprecated(fcn):
    @wraps(fcn)
    def new_fcn(*args, **kwargs):
        warnings.warn(
            f'{fcn.__name__} is deprecated and will be removed in the future', DeprecationWarning
        )
        return fcn(*args, **kwargs)
    return new_fcn
