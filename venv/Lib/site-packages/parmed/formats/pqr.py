"""
This module contains classes for reading and writing PQR files
"""
from contextlib import closing
import numpy as np
from .registry import FileFormatType
from .pdb import _standardize_resname, PDBFile, _is_hetatm
from ..exceptions import PDBError, PDBWarning
from ..periodic_table import AtomicNum, Mass, Element, element_by_name
from ..structure import Structure
from ..topologyobjects import Atom, ExtraPoint
from ..utils.io import genopen
import warnings

class PQRFile(metaclass=FileFormatType):
    """ Standard PDB file format parser and writer """
    #===================================================

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as a PDB file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is a PQR file
        """
        with closing(genopen(filename, 'r')) as f:
            for line in f:
                words = line.split()
                if not words:
                    continue
                elif words[0] in ('CRYST1', 'END', 'END', 'HEADER', 'NUMMDL',
                        'MASTER', 'AUTHOR', 'CAVEAT', 'COMPND', 'EXPDTA',
                        'MDLTYP', 'KEYWDS', 'OBSLTE', 'SOURCE', 'SPLIT',
                        'SPRSDE', 'TITLE ', 'ANISOU', 'CISPEP', 'CONECT',
                        'DBREF ', 'HELIX ', 'HET', 'LINK', 'MODRES',
                        'REVDAT', 'SEQADV', 'SHEET', 'SSBOND', 'FORMUL',
                        'HETNAM', 'HETSYN', 'SEQRES', 'SITE', 'ENDMDL', 'MODEL',
                        'JRNL', 'REMARK', 'TER', 'USER'):
                    continue
                elif line[:5] in ('ORIGX', 'SCALE', 'MTRIX'):
                    if line[5] not in '123':
                        return False
                elif words[0] in ('ATOM', 'HETATM'):
                    # Format is:
                    # rec atnum atname resname [chain] resnum x y z chg radius
                    # Where the chain ID is optional. rec must be ATOM or HETATM
                    if len(words) < 10:
                        return False
                    elif PDBFile.id_format(filename):
                        return False # It is a PDB file

                    if len(words) == 10:
                        offset = 0
                    elif len(words) >= 11:
                        offset = 1
                        try:
                            float(words[10])
                        except ValueError:
                            offset = 0
                    if not words[1].isdigit(): return False
                    if words[2].isdigit(): return False
                    if words[3].isdigit(): return False
                    if not words[4+offset].isdigit(): return False
                    try:
                        float(words[5+offset])
                        float(words[6+offset])
                        float(words[7+offset])
                        float(words[8+offset])
                        float(words[9+offset])
                    except ValueError:
                        return False
                    return True
                else:
                    return False
            return False

    #===================================================

    @staticmethod
    def parse(filename, skip_bonds=True):
        """ Read a PQR file and return a populated `Structure` class

        Parameters
        ----------
        filename : str or file-like
            Name of the PQR file to read, or a file-like object that can iterate
            over the lines of a PQR. Compressed file names can be specified and
            are determined by file-name extension (e.g., file.pqr.gz,
            file.pqr.bz2)
        skip_bonds : bool, optional
            If True, skip trying to assign bonds. This can save substantial time
            when parsing large files with non-standard residue names. However,
            no bonds are assigned. This is OK if, for instance, the PQR file is
            being parsed simply for its coordinates. Default is False.

        Returns
        -------
        structure : :class:`Structure`
            The Structure object initialized with all of the information from
            the PDB file.  No bonds or other topological features are added by
            default.
        """
        if isinstance(filename, str):
            own_handle = True
            fileobj = genopen(filename, 'r')
        else:
            own_handle = False
            fileobj = filename

        struct = Structure()
        # Add metadata fields
        modelno = 1 # For PDB files with multiple MODELs
        atomno = 0
        coordinates = []
        all_coordinates = []

        # Support hexadecimal numbering like that printed by VMD
        try:
            for line in fileobj:
                words = line.split()
                if words[0] in ('ATOM', 'HETATM'):
                    atomno += 1
                    if len(words) == 10:
                        _, num, nam, res, resn, x, y, z, chg, rad = words
                        chn = ''
                    elif len(words) >= 11:
                        _, num, nam, res, chn, resn, x, y, z, chg, rad = (
                            words[i] for i in range(11)
                        )
                        # If the radius is not a float (but rather a letter,
                        # like the element or something), then the chain might
                        # be missing. In this case, shift all tokens "back" one
                        # and empty the chn string
                        try:
                            float(rad)
                        except ValueError:
                            resn, x, y, z, chg, rad = chn, resn, x, y, z, chg
                    else:
                        raise ValueError('Illegal PQR record format: expected '
                                         '10 or 11 tokens on the atom line')
                    x, y, z = float(x), float(y), float(z)
                    chg, rad = float(chg), float(rad)
                    resn, num = int(resn), int(num)
                    elem = element_by_name(nam) # Yuck
                    atomic_number = AtomicNum[elem]
                    mass = Mass[elem]
                    if nam in ('EP', 'LP'): # lone pair
                        atom = ExtraPoint(atomic_number=atomic_number, name=nam,
                                          charge=chg, mass=mass, number=num,
                                          solvent_radius=rad)
                    else:
                        atom = Atom(atomic_number=atomic_number, name=nam,
                                    charge=chg, mass=mass, number=num,
                                    solvent_radius=rad)
                    atom.xx, atom.xy, atom.xz = float(x), float(y), float(z)
                    if modelno == 1:
                        struct.add_atom(atom, res, resn, chn)
                    else:
                        try:
                            orig_atom = struct.atoms[atomno-1]
                        except IndexError:
                            raise PDBError('Extra atom in MODEL %d' % modelno)
                        if (orig_atom.residue.name != res.strip()
                                or orig_atom.name != nam.strip()):
                            raise PDBError(
                                f'Atom {atomno} differs in MODEL {modelno} [{orig_atom.residue.name} '
                                f'{orig_atom.name} vs. {res} {nam}]'
                            )
                    coordinates.extend([atom.xx, atom.xy, atom.xz])
                elif words[0] == 'ENDMDL':
                    # End the current model
                    if len(struct.atoms) == 0:
                        raise PDBError('MODEL ended before any atoms read in')
                    modelno += 1
                    if len(struct.atoms)*3 != len(coordinates):
                        raise PDBError('Inconsistent atom numbers in some PDB models')
                    all_coordinates.append(coordinates)
                    atomno = 0
                    coordinates = []
                elif words[0] == 'MODEL':
                    if modelno == 1 and len(struct.atoms) == 0: continue
                    if len(coordinates) > 0:
                        if len(struct.atoms)*3 != len(coordinates):
                            raise PDBError('Inconsistent atom numbers in some PDB models')
                        warnings.warn('MODEL not explicitly ended', PDBWarning)
                        all_coordinates.append(coordinates)
                        coordinates = []
                    modelno += 1
                    atomno = 0
                elif words[0] == 'CRYST1':
                    a, b, c = (float(w) for w in words[1:4])
                    try:
                        A, B, C = (float(w) for w in words[4:7])
                    except ValueError:
                        A = B = C = 90.0
                    struct.box = [a, b, c, A, B, C]
        finally:
            if own_handle:
                fileobj.close()

        struct.unchange()
        if not skip_bonds:
            struct.assign_bonds()
        if coordinates:
            if len(coordinates) != 3 * len(struct.atoms):
                raise PDBError('bad number of atoms in some PQR models')
            all_coordinates.append(coordinates)
        struct._coordinates = np.array(all_coordinates).reshape((-1, len(struct.atoms), 3))
        return struct

    #===================================================

    @staticmethod
    def write(struct, dest, renumber=True, coordinates=None,
              standard_resnames=False):
        """ Write a PDB file from a Structure instance

        Parameters
        ----------
        struct : :class:`Structure`
            The structure from which to write the PDB file
        dest : str or file-like
            Either a file name or a file-like object containing a `write`
            method to which to write the PDB file. If it is a filename that
            ends with .gz or .bz2, a compressed version will be written using
            either gzip or bzip2, respectively.
        renumber : bool, optional
            If True, renumber the atoms and residues sequentially as they are
            stored in the structure.  If False, use the original numbering if
            it was assigned previously. Default is True
        coordinates : array-like of float, optional
            If provided, these coordinates will be written to the PDB file
            instead of the coordinates stored in the structure. These
            coordinates should line up with the atom order in the structure
            (not necessarily the order of the "original" PDB file if they
            differ)
        standard_resnames : bool, optional
            If True, common aliases for various amino and nucleic acid residues
            will be converted into the PDB-standard values. Default is False
        """
        own_handle = False
        if not hasattr(dest, 'write'):
            dest = genopen(dest, 'w')
            own_handle = True
        atomrec = ('ATOM  %5d %-3s  %-3s %1s %3d    %7.3f %7.3f %7.3f %8.4f '
                   '%8.4f\n')
        hetatomrec = atomrec.replace('ATOM  ', 'HETATM')
        if struct.box is not None:
            dest.write('CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f\n' % (
                    struct.box[0], struct.box[1], struct.box[2], struct.box[3],
                    struct.box[4], struct.box[5]))
        if coordinates is not None:
            coords = np.asanyarray(coordinates)
            try:
                coords = coords.reshape((-1, len(struct.atoms), 3))
            except ValueError:
                raise TypeError("Coordinates has unexpected shape")
        else:
            coords = struct.get_coordinates('all')
        # Create a function to process each atom and return which one we want
        # to print, based on our alternate location choice
        if standard_resnames:
            standardize = lambda x: _standardize_resname(x)
        else:
            standardize = lambda x: (x, _is_hetatm(x))
        last_number = 0
        last_rnumber = 0
        for model, coord in enumerate(coords):
            if coords.shape[0] > 1:
                dest.write('MODEL      %5d\n' % (model+1))
            for res in struct.residues:
                if renumber:
                    atoms = res.atoms
                else:
                    atoms = sorted(res.atoms, key=lambda atom: atom.number)
                for atom in atoms:
                    # Figure out the serial numbers we want to print
                    if renumber:
                        anum = (atom.idx + 1)
                        rnum = (res.idx + 1)
                    else:
                        anum = (atom.number or last_number + 1)
                        rnum = (atom.residue.number or last_rnumber + 1)
                    last_number = anum
                    last_rnumber = rnum
                    # Do any necessary name munging to respect the PDB spec
                    if (len(atom.name) < 4 and
                            len(Element[atom.atomic_number]) != 2):
                        aname = ' %-3s' % atom.name
                    else:
                        aname = atom.name
                    xyz = coord[atom.idx]
                    resname, hetatm = standardize(res.name)
                    if hetatm:
                        rec = hetatomrec
                    else:
                        rec = atomrec
                    dest.write(rec % (anum, aname, resname, res.chain, rnum,
                                      xyz[0], xyz[1], xyz[2], atom.charge,
                                      atom.solvent_radius))
            if coords.shape[0] > 1:
                dest.write('ENDMDL\n')

        dest.write("%-80s\n" % "END")
        if own_handle:
            dest.close()
