"""
This module contains classes for reading various TINKER-style files
"""
from ..exceptions import TinkerError
from ..formats.registry import FileFormatType, load_file
from ..periodic_table import element_by_name, AtomicNum
from ..structure import Structure
from ..topologyobjects import Atom, Bond, Residue
from ..utils.io import genopen

class KeywordControlFile:
    """ Reads and processes a keyword control file for TINKER simulations """
    _datatypes = {'PARAMETERS' : str, 'A-AXIS' : float, 'B-AXIS' : float,
                  'C-AXIS' : float, 'ALPHA' : float, 'BETA' : float,
                  'GAMMA' : float}

    def __init__(self, fname):
        # The currently recognized/parsed keywords (these are the only ones
        # necessary for running Amoeba in Amber)
        self.keywords = {'PARAMETERS' : None, 'A-AXIS' : None, 'B-AXIS' : None,
                         'C-AXIS' : None, 'ALPHA' : None, 'BETA' : None, 'GAMMA' : None}
        # Parse the file
        for line in open(fname, 'r'):
            # Skip over any blank lines
            if not line.strip(): continue
            words = line.split()
            key = words[0].upper()
            # Get rid of the keyword and all whitespace
            result = line.replace(words[0], '').strip()
            try:
                self.keywords[key] = self._datatypes[key](result)
            except KeyError:
                pass # All non-keyword lines are ignored comments
            except ValueError as err:
                raise TinkerError(
                    'Malformed keyword control file! Could not convert the value of '
                    f'{key} into type {self._datatypes[key]}'
                ) from err

class XyzFile(Structure, metaclass=FileFormatType):
    """ Reads and processes a Tinker XYZ file

    Parameters
    ----------
    fname : str or file-like
        Name of the file, or the file object containing the XYZ file contents
    seq : str, optional
        Name of the file containing the residue (and chain) sequence. Default is
        None (so every atom will be part of the same residue)
    """

    @staticmethod
    def _check_atom_record(words):
        """ Checks that the tokenized line is consistent with an atom record

        Parameters
        ----------
        words : list of str
            Result of line.split() on a line of the file

        Returns
        -------
        is_atom : bool
            If it is consistent with an atom record, return True. Otherwise,
            False
        """
        if len(words) < 6:
            return False
        # First token is the atom index
        try:
            int(words[0])
        except ValueError:
            return False
        # Second token is the atom name (NOT numeric)
        try:
            float(words[1])
        except ValueError:
            try:
                int(words[1])
            except ValueError:
                pass
            else:
                return False
        else:
            return False
        # Third, fourth, and fifth tokens are floats, but NOT integers
        try:
            [float(w) for w in words[2:5]]
        except ValueError:
            return False
        else:
            try:
                [int(w) for w in words[2:5]]
            except ValueError:
                pass
            else:
                return False
        # Remaining numbers are integers
        try:
            [int(w) for w in words[5:]]
        except ValueError:
            return False
        return True

    @staticmethod
    def id_format(filename):
        """ Identify the file as a Tinker XYZ file

        Parameters
        ----------
        filename : str
            Name of the file to test whether or not it is a mol2 file

        Returns
        -------
        is_fmt : bool
            True if it is a xyz file, False otherwise
        """
        f = genopen(filename, 'r')
        words = f.readline().split() # natom and title
        if not words:
            return False
        try:
            natom = int(words[0])
        except (ValueError, IndexError):
            return False
        else:
            if natom <= 0:
                return False
        words = f.readline().split()
        # Either a box line or a line with the first atom
        if len(words) == 6:
            try:
                [float(w) for w in words]
            except ValueError:
                if XyzFile._check_atom_record(words):
                    return True
                return False
            else:
                # Next line needs to be an atom record
                words = f.readline().split()
        return XyzFile._check_atom_record(words)

    def __init__(self, fname, seq=None):
        super(XyzFile, self).__init__()
        if isinstance(fname, str):
            fxyz = genopen(fname, 'r')
            own_handle_xyz = True
        else:
            fxyz = fname
            own_handle_xyz = False
        if seq is not None:
            seqstruct = load_file(seq)
        # Now parse the file
        try:
            natom = int(fxyz.readline().split()[0])
        except (ValueError, IndexError):
            raise TinkerError('Bad XYZ file format; first line')
        if seq is not None and natom != len(seqstruct.atoms):
            raise ValueError(f'Sequence file {seq} # of atoms does not match the # of atoms in the XYZ file')
        words = fxyz.readline().split()
        if len(words) == 6 and not XyzFile._check_atom_record(words):
            self.box = [float(w) for w in words]
            words = fxyz.readline().split()
        residue = Residue('SYS')
        residue.number = 1
        residue._idx = 0
        if seq is not None:
            residue = seqstruct.residues[0]
            atomic_number = _guess_atomic_number(words[1], residue)
        else:
            atomic_number = AtomicNum[element_by_name(words[1])]
        atom = Atom(atomic_number=atomic_number, name=words[1], type=words[5])
        atom.xx, atom.xy, atom.xz = [float(w) for w in words[2:5]]
        self.add_atom(
            atom, residue.name, residue.number, residue.chain, residue.insertion_code, residue.segid
        )
        bond_ids = [[int(w) for w in words[6:]]]
        for i, line in enumerate(fxyz):
            words = line.split()
            if seq is not None:
                residue = seqstruct.atoms[i+1].residue
                atomic_number = _guess_atomic_number(words[1], residue)
            else:
                atomic_number = AtomicNum[element_by_name(words[1])]
            atom = Atom(atomic_number=atomic_number, name=words[1], type=words[5])
            atom.xx, atom.xy, atom.xz = [float(w) for w in words[2:5]]
            self.add_atom(
                atom, residue.name, residue.number, residue.chain, residue.insertion_code, residue.segid
            )
            bond_ids.append([int(w) for w in words[6:]])
        # All of the bonds are stored now -- go ahead and make them now
        for atom, bonds in zip(self.atoms, bond_ids):
            i = atom.idx + 1
            for idx in bonds:
                if idx > i:
                    self.bonds.append(Bond(atom, self.atoms[idx-1]))
        if seq is None:
            # Try to improve atomic number prediction for monoatomic species
            # (like ions) if no sequence as loaded
            for atom in self.atoms:
                if len(atom.bonds) == 0: # not bonded to anybody else
                    atom.atomic_number = _guess_atomic_number(atom.name)
        if own_handle_xyz:
            fxyz.close()

class DynFile:
    """ Reads and processes a Tinker DYN file """
    def __init__(self, fname=None):
        if fname is not None:
            self.read(fname)

    def read(self, fname):
        """ Parses the .dyn file """
        f = open(fname, 'r')
        try:
            if f.readline().strip() != 'Number of Atoms and Title :':
                raise TinkerError(f'{fname} is not a recognized TINKER .dyn file')
            line = f.readline()
            self.natom, self.title = int(line[0:8].strip()), line[8:].strip()
            if f.readline().strip() != 'Periodic Box Dimensions :':
                raise TinkerError('No periodic box dimension line')
            self.box = [0.0 for i in range(6)]
            words = f.readline().upper().split()
            self.box[0] = float(words[0].replace('D', 'E'))
            self.box[1] = float(words[1].replace('D', 'E'))
            self.box[2] = float(words[2].replace('D', 'E'))
            words = f.readline().upper().split()
            self.box[3] = float(words[0].replace('D', 'E'))
            self.box[4] = float(words[1].replace('D', 'E'))
            self.box[5] = float(words[2].replace('D', 'E'))
            if f.readline().strip() != 'Current Atomic Positions :':
                raise TinkerError('No atomic positions in .dyn file')
            self.positions = [[0.0, 0.0, 0.0] for i in range(self.natom)]
            DynFile._read_section(f, self.positions, self.natom)
            line = f.readline().strip()
            if line == 'Current Translational Velocities :':
                self.rigidbody = True
                self.translational_velocities = [[0.0, 0.0, 0.0] for i in range(self.natom)]
                DynFile._read_section(f, self.translational_velocities, self.natom)
                if f.readline().strip() != 'Current Angular Velocities :':
                    raise TinkerError('Could not find Angular velocity section in .dyn file')
                self.angular_velocities = [[0.0, 0.0, 0.0] for i in range(self.natom)]
                DynFile._read_section(f, self.angular_velocities, self.natom)
                if f.readline().strip() != 'Current Angular Momenta :':
                    raise TinkerError('Could not find angular momenta section in .dyn file')
                self.angular_momenta = [[0.0, 0.0, 0.0] for i in range(self.natom)]
                DynFile._read_section(f, self.angular_momenta, self.natom)
            elif line == 'Current Atomic Velocities :':
                self.rigidbody = False
                self.velocities = [[0.0, 0.0, 0.0] for i in range(self.natom)]
                DynFile._read_section(f, self.velocities, self.natom)
                if f.readline().strip() != 'Current Atomic Accelerations :':
                    raise TinkerError(f'Could not find accelerations in {fname}')
                self.accelerations = [[0.0, 0.0, 0.0] for i in range(self.natom)]
                DynFile._read_section(f, self.accelerations, self.natom)
                if f.readline().strip() != 'Alternate Atomic Accelerations :':
                    raise TinkerError(f'Could not find old accelerations in {fname}')
                self.old_accelerations = [[0.0, 0.0, 0.0] for i in range(self.natom)]
                DynFile._read_section(f, self.old_accelerations, self.natom)
            else:
                raise TinkerError(f'No velocities in {fname}')
        finally:
            f.close()

    @staticmethod
    def _read_section(f, container, natom):
        """
        Reads a section of an open file into the passed container. The
        container should be a 2-dimensional list of length natom x 3
        """
        for i in range(natom):
            words = f.readline().upper().split()
            try:
                container[i][0] = float(words[0].replace('D', 'E'))
                container[i][1] = float(words[1].replace('D', 'E'))
                container[i][2] = float(words[2].replace('D', 'E'))
            except (IndexError, ValueError):
                raise TinkerError('Could not parse values from dyn file')

def is_float(thing):
    try:
        float(thing)
        return True
    except ValueError:
        return False

def _guess_atomic_number(name, residue=None):
    """ Guesses the atomic number """
    # Special-case single-atom residues, which are almost always ions
    name = ''.join(c for c in name if c.isalpha())
    if residue is None or len(residue.atoms) == 1:
        if len(name) > 1:
            try:
                return AtomicNum[name[0].upper() + name[1].lower()]
            except KeyError:
                return AtomicNum[element_by_name(name)]
    return AtomicNum[element_by_name(name)]
