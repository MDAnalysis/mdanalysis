""" A collection of utilities for use with Pandas objects """
import numpy as np
try:
    import pandas as pd
except ImportError:
    pd = None
from ..exceptions import ParameterWarning
import warnings

# Utility function for generating a DataFrame based on a collection of items.
# The passed object must have an `atoms` attribute

def create_dataframe(obj):
    """ Creates a pandas.DataFrame object from the current instance

    Returns
    -------
    dataframe : :class:`pandas.DataFrame`

    Notes
    -----
    The DataFrame will be over all atoms. The columns will be the attributes
    of the atom (as well as its containing residue). Some columns will
    *always* exist. Others will only exist if those attributes have been set
    on the Atom instances (see the :class:`Atom` docs for possible
    attributes and their meaning). The columns that will always be present
    are:

        - number : int
        - name : str
        - type : str
        - atomic_number : int
        - charge : float
        - mass : float
        - nb_idx : int
        - solvent_radius : float
        - screen : float
        - occupancy : float
        - bfactor : float
        - altloc : str
        - tree : str
        - join : int
        - irotat : int
        - rmin : float
        - epsilon : float
        - rmin_14 : float
        - epsilon_14 : float
        - resname : str (name of the containing residue)
        - resid : int (Sequential index of the containing residue)
        - resnum : int (original residue number in the input structure)
        - chain : str (chain ID that the containing residue belongs to)

    The following attributes are optionally present if they were present in
    the original file defining the structure:

        - xx : float (x-coordinate position)
        - xy : float (y-coordinate position)
        - xz : float (z-coordinate position)
        - vx : float (x-coordinate velocity)
        - vy : float (y-coordinate velocity)
        - vz : float (z-coordinate velocity)
        - type_idx : int (integer type index for AMOEBA)
        - class_idx : int (integer class type index for AMOEBA)
        - multipole_111 : float (Monopole)
        - multipole_211 : float (1,1 Dipole component)
        - multipole_212 : float (1,2 Dipole component)
        - multipole_222 : float (2,2 Dipole component)
        - multipole_411 : float (1,1 Quadrupole component)
        - multipole_412 : float (1,2 Quadrupole component)
        - multipole_422 : float (2,2 Quadrupole component)
        - multipole_413 : float (1,3 Quadrupole component)
        - multipole_423 : float (2,3 Quadrupole component)
        - multipole_433 : float (3,3 Quadrupole component)
        - polarizability : float (dipole polarizability)
        - vdw_parent : int (index of the vdW parent atom of this atom)
        - U11 : float (U[1][1] of anisotropic b-factor tensor)
        - U22 : float (U[2][2] of anisotropic b-factor tensor)
        - U33 : float (U[3][3] of anisotropic b-factor tensor)
        - U12 : float (U[1][2] of anisotropic b-factor tensor)
        - U13 : float (U[1][3] of anisotropic b-factor tensor)
        - U23 : float (U[2][3] of anisotropic b-factor tensor)
    """
    if pd is None:
        raise ImportError('pandas is not available; cannot create a DataFrame')
    ret = pd.DataFrame()

    atoms = obj.atoms

    ret['number'] = [atom.number for atom in atoms]
    ret['name'] = [atom.name for atom in atoms]
    ret['type'] = [atom.type for atom in atoms]
    ret['atomic_number'] = [atom.atomic_number for atom in atoms]
    ret['charge'] = [atom.charge for atom in atoms]
    ret['mass'] = [atom.mass for atom in atoms]
    ret['nb_idx'] = [atom.nb_idx for atom in atoms]
    ret['solvent_radius'] = [atom.solvent_radius for atom in atoms]
    ret['screen'] = [atom.screen for atom in atoms]
    ret['occupancy'] = [atom.occupancy for atom in atoms]
    ret['bfactor'] = [atom.bfactor for atom in atoms]
    ret['altloc'] = [atom.altloc for atom in atoms]
    ret['tree'] = [atom.tree for atom in atoms]
    ret['join'] = [atom.join for atom in atoms]
    ret['irotat'] = [atom.irotat for atom in atoms]
    ret['rmin'] = [atom.rmin for atom in atoms]
    ret['epsilon'] = [atom.epsilon for atom in atoms]
    ret['rmin_14'] = [atom.rmin_14 for atom in atoms]
    ret['epsilon_14'] = [atom.epsilon_14 for atom in atoms]
    ret['resname'] = [atom.residue.name for atom in atoms]
    ret['resid'] = [atom.residue.idx for atom in atoms]
    ret['resnum'] = [atom.residue.number for atom in atoms]
    ret['chain'] = [atom.residue.chain for atom in atoms]
    ret['segid'] = [atom.residue.segid for atom in atoms]

    # Now for optional attributes
    # Coordinates
    try:
        coords = pd.DataFrame(
                [[atom.xx, atom.xy, atom.xz] for atom in atoms],
                columns=['xx', 'xy', 'xz']
        )
    except AttributeError:
        pass
    else:
        ret = ret.join(coords)
    # Velocities
    try:
        vels = pd.DataFrame(
                [[atom.vx, atom.vy, atom.vz] for atom in atoms],
                columns=['vx', 'vy', 'vz']
        )
    except AttributeError:
        pass
    else:
        ret = ret.join(vels)
    # AMOEBA LJ type
    try:
        ret['type_idx'] = [atom.type_idx for atom in atoms]
    except AttributeError:
        pass
    # AMOEBA class type
    try:
        ret['class_idx'] = [atom.class_idx for atom in atoms]
    except AttributeError:
        pass
    # AMOEBA multipoles
    try:
        multipoles = pd.DataFrame(
            [atom.multipoles for atom in atoms],
            columns=['multipole_111', 'multipole_211', 'multipole_212',
                     'multipole_222', 'multipole_411', 'multipole_412',
                     'multipole_422', 'multipole_413', 'multipole_423',
                     'multipole_433']
        )
    except AttributeError:
        pass
    else:
        ret = ret.join(multipoles)
    # AMOEBA polarizabilities
    try:
        ret['polarizability'] = [atom.polarizability for atom in atoms]
    except AttributeError:
        pass
    # AMOEBA vdw parent atom
    try:
        ret['vdw_parent'] = [atom.vdw_parent.idx for atom in atoms]
    except AttributeError:
        pass
    # anisotropic b-factors
    none6 = [None] * 6
    anisos = [atom.anisou for atom in atoms]
    for i, aniso in enumerate(anisos):
        if hasattr(aniso, 'tolist'):
            anisos[i] = aniso.tolist()
    all_nones = True
    for i, aniso in enumerate(anisos):
        if aniso is None:
            anisos[i] = none6
        elif all_nones:
            all_nones = False
    if not all_nones:
        ret = ret.join(pd.DataFrame(anisos, columns=['U11', 'U22', 'U33', 'U12', 'U13', 'U23']))
    return ret

def load_dataframe(obj, dataframe):
    """
    Loads a DataFrame into the current object, setting atomic properties based
    on the entries of the DataFrame.  Supported atomic properties are:

        - number : int
        - name : str
        - type : str
        - atomic_number : int
        - charge : float
        - mass : float
        - nb_idx : int
        - solvent_radius : float
        - screen : float
        - occupancy : float
        - bfactor : float
        - altloc : str
        - tree : str
        - join : int
        - irotat : int
        - rmin : float
        - epsilon : float
        - rmin_14 : float
        - epsilon_14 : float
        - xx : float (x-coordinate position)
        - xy : float (y-coordinate position)
        - xz : float (z-coordinate position)
        - vx : float (x-coordinate velocity)
        - vy : float (y-coordinate velocity)
        - vz : float (z-coordinate velocity)
        - type_idx : int (integer type index for AMOEBA)
        - class_idx : int (integer class type index for AMOEBA)
        - multipole_111 : float (Monopole)
        - multipole_211 : float (1,1 Dipole component)
        - multipole_212 : float (1,2 Dipole component)
        - multipole_222 : float (2,2 Dipole component)
        - multipole_411 : float (1,1 Quadrupole component)
        - multipole_412 : float (1,2 Quadrupole component)
        - multipole_422 : float (2,2 Quadrupole component)
        - multipole_413 : float (1,3 Quadrupole component)
        - multipole_423 : float (2,3 Quadrupole component)
        - multipole_433 : float (3,3 Quadrupole component)
        - polarizability : float (dipole polarizability)
        - vdw_parent : int (index of the vdW parent atom of this atom)
        - segid : segment ID (similar to chain, but for CHARMM)
        - U11 : float (U[1][1] of anisotropic b-factor tensor)
        - U22 : float (U[2][2] of anisotropic b-factor tensor)
        - U33 : float (U[3][3] of anisotropic b-factor tensor)
        - U12 : float (U[1][2] of anisotropic b-factor tensor)
        - U13 : float (U[1][3] of anisotropic b-factor tensor)
        - U23 : float (U[2][3] of anisotropic b-factor tensor)

    The resname, resid, chain, and resnum attributes are ignored. Other
    attributes emit a ParameterWarning
    """
    atoms = obj.atoms

    def set_attribute(attr, data):
        """ Set the attribute list from the data """
        if len(data) != len(atoms):
            raise ValueError('Data does not match length of atoms list')
        for atom, x in zip(atoms, data):
            setattr(atom, attr, x)
    def set_residue_attr(attr, data):
        if len(data) != len(atoms):
            raise ValueError('Data does not match length of atoms list')
        for atom, x in zip(atoms, data):
            setattr(atom.residue, attr, x)

    multipoles = [None for i in range(10)]
    anisous = [None for i in range(6)]
    for key, data in dataframe.items():
        if key == 'number':
            set_attribute('number', data)
        elif key == 'name':
            set_attribute('name', data)
        elif key == 'type':
            set_attribute('type', data)
        elif key == 'atomic_number':
            set_attribute('atomic_number', data)
        elif key == 'charge':
            set_attribute('charge', data)
        elif key == 'mass':
            set_attribute('mass', data)
        elif key == 'nb_idx':
            set_attribute('nb_idx', data)
        elif key == 'solvent_radius':
            set_attribute('solvent_radius', data)
        elif key == 'screen':
            set_attribute('screen', data)
        elif key == 'occupancy':
            set_attribute('occupancy', data)
        elif key == 'bfactor':
            set_attribute('bfactor', data)
        elif key == 'altloc':
            set_attribute('altloc', data)
        elif key == 'tree':
            set_attribute('tree', data)
        elif key == 'join':
            set_attribute('join', data)
        elif key == 'irotat':
            set_attribute('irotat', data)
        elif key == 'rmin':
            set_attribute('rmin', data)
        elif key == 'epsilon':
            set_attribute('epsilon', data)
        elif key == 'rmin_14':
            set_attribute('rmin_14', data)
        elif key == 'epsilon_14':
            set_attribute('epsilon_14', data)
        elif key == 'xx':
            set_attribute('xx', data)
        elif key == 'xy':
            set_attribute('xy', data)
        elif key == 'xz':
            set_attribute('xz', data)
        elif key == 'vx':
            set_attribute('vx', data)
        elif key == 'vy':
            set_attribute('vy', data)
        elif key == 'vz':
            set_attribute('vz', data)
        elif key == 'type_idx':
            set_attribute('type_idx', data)
        elif key == 'class_idx':
            set_attribute('class_idx', data)
        elif key == 'vdw_parent':
            if len(data) != len(atoms):
                raise ValueError('vdw_parent length not equal to natom')
            for atom, parent in zip(atoms, data):
                atom.vdw_parent = atoms[parent]
        elif key == 'polarizability':
            set_attribute('polarizability', data)
        elif key == 'multipole_111':
            multipoles[0] = data
        elif key == 'multipole_211':
            multipoles[1] = data
        elif key == 'multipole_212':
            multipoles[2] = data
        elif key == 'multipole_213':
            multipoles[3] = data
        elif key == 'multipole_411':
            multipoles[4] = data
        elif key == 'multipole_412':
            multipoles[5] = data
        elif key == 'multipole_422':
            multipoles[6] = data
        elif key == 'multipole_413':
            multipoles[7] = data
        elif key == 'multipole_423':
            multipoles[8] = data
        elif key == 'multipole_433':
            multipoles[9] = data
        elif key == 'U11':
            anisous[0] = data
        elif key == 'U22':
            anisous[1] = data
        elif key == 'U33':
            anisous[2] = data
        elif key == 'U12':
            anisous[3] = data
        elif key == 'U13':
            anisous[4] = data
        elif key == 'U23':
            anisous[5] = data
        elif key in ('resname', 'resid', 'resnum', 'chain', 'segid'):
            set_residue_attr(key, data)
            continue
        else:
            warnings.warn(f'Atomic property {key} not recognized', ParameterWarning)

    # Now combine the multipoles and anisous if they are all specified
    if not None in anisous:
        set_attribute('anisou', np.vstack(anisous).T)
    if not None in multipoles:
        set_attribute('multipoles', np.vstack(multipoles).T)
