import warnings

from .groups import (Atom, AtomGroup, Residue, ResidueGroup,
                     Segment, SegmentGroup)
from . import universe


def deprecate_class(class_new, message):
    """utility to deprecate a class"""
    class new_class(class_new):
        def __init__(self, *args, **kwargs):
            super(new_class, self).__init__(*args, **kwargs)
            warnings.warn(message, DeprecationWarning)
    return new_class


Universe = deprecate_class(universe.Universe,
                           "MDAnalysis.core.AtomGroup.Universe has been removed."
                           "Please use MDAnalaysis.core.universe.Universe."
                           "This stub will be removed in 1.0")

_group_message = ("MDAnalysis.core.AtomGroup.{0} has been removed."
                  "Please use MDAnalaysis.groups.{0}"
                  "This stub will be removed in 1.0")

Atom = deprecate_class(Atom, message=_group_message.format('Atom'))
AtomGroup = deprecate_class(AtomGroup,
                            message=_group_message.format('AtomGroup'))

Residue = deprecate_class(Residue, message=_group_message.format('Residue'))
ResidueGroup = deprecate_class(ResidueGroup,
                               message=_group_message.format('ResidueGroup'))

Segment = deprecate_class(Segment, message=_group_message.format('Segment'))
SegmentGroup = deprecate_class(SegmentGroup,
                               message=_group_message.format('SegmentGroup'))

__all__ = ['Universe', 'Atom', 'AtomGroup', 'Residue', 'ResidueGroup',
           'Segment', 'SegmentGroup']
