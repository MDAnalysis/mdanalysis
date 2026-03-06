"""
This module is simply a namespace for all Amber topology-like classes. There are
different variants; for instance the standard Amber topology, the chamber-made
topology, and the tinker_to_amber-made topology. These classes are all defined
in their own private modules, but imported here to simplify the API.
"""

from ..constants import PrmtopPointers
from .amberformat import AmberFormat
from ._amberparm import AmberParm, Rst7
from ._chamberparm import ChamberParm, ConvertFromPSF
from ._tinkerparm import AmoebaParm, BeemanRestart

__all__ = ['AmberFormat', 'AmberParm', 'ChamberParm', 'LoadParm', 'Rst7']

# Supply a function to load a topology file in the 'correct' format
def LoadParm(parmname, xyz=None, box=None):
    """
    Loads a topology file using the correct class.

    Parameters
    ----------
    parmname : ``str``
        The name of the topology file to load
    xyz : str or array, optional
        If provided, the coordinates, velocities and unit cell dimensions from the provided
        Amber inpcrd/restart file will be loaded into the molecule, or the
        coordinates will be loaded from the coordinate array
    box : array, optional
        If provided, the unit cell information will be set from the provided
        unit cell dimensions (a, b, c, alpha, beta, and gamma, respectively)

    Returns
    -------
    parm : :class:`AmberParm` (or subclass)
        This function parses the topology file, determines if it is an
        Amber-style (i.e., *traditional* Amber force field), Chamber-style
        (i.e., CHARMM force field), or Amoeba-style (i.e., Amoeba force field),
        and then returns an instance of the appropriate type.
    """
    from .. import load_file
    parm = AmberFormat(parmname)
    if 'CTITLE' in parm.flag_list:
        parm = parm.view_as(ChamberParm)
    elif 'AMOEBA_FORCEFIELD' in parm.flag_list:
        parm = parm.view_as(AmoebaParm)
    else:
        parm = parm.view_as(AmberParm)

    if isinstance(xyz, str):
        f = load_file(xyz)
        if not hasattr(f, 'coordinates') or f.coordinates is None:
            raise TypeError(f'{xyz} does not have coordinates')
        parm.coordinates = f.coordinates
        if hasattr(f, 'velocities') and f.velocities is not None:
            parm.velocities = f.velocities
        if hasattr(f, 'box') and f.box is not None and box is None:
            parm.box = f.box

    else:
        parm.coordinates = xyz
    if box is not None:
        parm.box = box

    # If all else fails, set the box from the prmtop file
    if parm.parm_data['POINTERS'][PrmtopPointers.IFBOX] > 0 and parm.box is None:
        box = parm.parm_data['BOX_DIMENSIONS']
        parm.box = list(box[1:]) + [box[0], box[0], box[0]]

    parm.hasbox = parm.box is not None

    return parm
