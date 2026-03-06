"""
This module contains functionality relevant to loading and parsing GROMACS GRO
(coordinate) files and building a stripped-down Structure from it
"""
from contextlib import closing
from ..constants import TINY
from ..exceptions import GromacsError
from ..formats.registry import FileFormatType
from ..geometry import (
    box_vectors_to_lengths_and_angles, box_lengths_and_angles_to_vectors, reduce_box_vectors
)
from ..periodic_table import AtomicNum, element_by_name, Mass
from ..structure import Structure
from ..topologyobjects import Atom, ExtraPoint
from .. import unit as u
from ..utils.io import genopen

class _AtomLineParser:
    """ Parses atom lines from GRO files """
    def __init__(self):
        self._digits = None
        self._pdeci = 0
        self._ndeci = 0

    def read(self, line):
        """ Reads a line

        Parameters
        ----------
        line : str
            A line with an atom record from a GRO file

        Returns
        -------
        atom, resname, resnum : Atom, str, int
            The Atom instance, residue name, and residue number containing the
            atom
        """
        resnum = int(line[:5])
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        elem = element_by_name(atomname)
        atomic_number = AtomicNum[elem]
        mass = Mass[elem]
        atnum = int(line[15:20])
        if atomic_number == 0:
            atom = ExtraPoint(name=atomname, number=atnum)
        else:
            atom = Atom(atomic_number=atomic_number, name=atomname,
                        number=atnum, mass=mass)
        if self._digits is None:
            self._pdeci = line.index('.', 20)
            self._ndeci = line.index('.', self._pdeci+1)
            self._digits = self._ndeci - self._pdeci
        atom.xx, atom.xy, atom.xz = (
            float(line[20+i*self._digits:20+(i+1)*self._digits])*10 for i in range(3)
        )
        wbeg = 20 + self._digits * 3
        wend = wbeg + self._digits
        if line[wbeg:wend].strip():
            atom.vx, atom.vy, atom.vz = (
                float(line[wbeg+i*self._digits:wend+i*self._digits])*10 for i in range(3)
            )
        return atom, resname, resnum

class GromacsGroFile(metaclass=FileFormatType):
    """ Parses and writes Gromacs GRO files """

    @staticmethod
    def id_format(filename):
        """ Identifies the file as a GROMACS GRO file

        Parameters
        ----------
        filename : str
            Name of the file to check if it is a Gromacs GRO file

        Returns
        -------
        is_fmt : bool
            If it is identified as a Gromacs GRO file, return True. False
            otherwise
        """
        with closing(genopen(filename)) as f:
            f.readline() # Title line
            try:
                int(f.readline().strip()) # number of atoms
            except ValueError:
                return False
            line = f.readline()
            try:
                int(line[:5])
                if not line[5:10].strip():
                    return False
                if not line[10:15].strip():
                    return False
                int(line[15:20])
                pdeci = [i for i, x in enumerate(line) if x == '.']
                ndeci = pdeci[1] - pdeci[0] - 5
                for i in range(1, 4):
                    wbeg = (pdeci[0]-4)+(5+ndeci)*(i-1)
                    wend = (pdeci[0]-4)+(5+ndeci)*i
                    float(line[wbeg:wend])
                i = 4
                wbeg = (pdeci[0]-4)+(5+ndeci)*(i-1)
                wend = (pdeci[0]-4)+(5+ndeci)*i
                if line[wbeg:wend].strip():
                    for i in range(4, 7):
                        wbeg = (pdeci[0]-4)+(5+ndeci)*(i-1)
                        wend = (pdeci[0]-4)+(5+ndeci)*i
                        float(line[wbeg:wend])
            except ValueError:
                return False
            return True

    @staticmethod
    def parse(filename, skip_bonds=False):
        """ Parses a Gromacs GRO file

        Parameters
        ----------
        filename : str or file-like
            Name of the file or the GRO file object
        skip_bonds : bool, optional
            If True, skip trying to assign bonds. This can save substantial time
            when parsing large files with non-standard residue names. However,
            no bonds are assigned. This is OK if, for instance, the GRO file is
            being parsed simply for its coordinates. This will also reduce the
            accuracy of assigned atomic numbers for typical ions. Default is
            False.

        Returns
        -------
        struct : :class:`Structure`
            The Structure instance instantiated with *just* residues and atoms
            populated (with coordinates)
        """
        struct = Structure()
        if isinstance(filename, str):
            fileobj = genopen(filename, 'r')
            own_handle = True
        else:
            fileobj = filename
            own_handle = False
        try:
            # Ignore the title line
            fileobj.readline()
            try:
                natom = int(fileobj.readline().strip())
            except ValueError:
                raise GromacsError(f'Could not parse {filename} as GRO file')
            line_parser = _AtomLineParser()
            for i, line in enumerate(fileobj):
                if i == natom: break
                try:
                    atom, resname, resnum = line_parser.read(line)
                except (ValueError, IndexError):
                    raise GromacsError(f'Could not parse the atom record of GRO file {filename}')
                struct.add_atom(atom, resname, resnum)
            else:
                # If no box exists, the break did not hit, so line still
                # contains the last atom (which cannot be interpreted as a box).
                # This wipes out line (IFF fileobj reached the line)
                line = fileobj.readline()
                if i+1 != natom:
                    raise GromacsError(f'Truncated GRO file. Found {i+1} of {natom} atoms')
            # Get the box from the last line if it's present
            if line.strip():
                try:
                    box = [float(x) for x in line.split()]
                except ValueError:
                    raise GromacsError(f'Could not understand box line of GRO file {filename}')
                if len(box) == 3:
                    struct.box = [box[0]*10, box[1]*10, box[2]*10, 90.0, 90.0, 90.0]
                elif len(box) == 9:
                    # Assume we have vectors
                    leng, ang = box_vectors_to_lengths_and_angles(
                        [box[0], box[3], box[4]] * u.nanometers,
                        [box[5], box[1], box[6]] * u.nanometers,
                        [box[7], box[8], box[2]] * u.nanometers,
                    )
                    a, b, c = leng.value_in_unit(u.angstroms)
                    alpha, beta, gamma = ang.value_in_unit(u.degrees)
                    struct.box = [a, b, c, alpha, beta, gamma]
        finally:
            if own_handle:
                fileobj.close()

        # Assign bonds (and improved element guesses)
        if not skip_bonds:
            struct.assign_bonds()

        return struct

    @staticmethod
    def write(struct, dest, precision=3, nobox=False, combine=False):
        """ Write a Gromacs Topology File from a Structure

        Parameters
        ----------
        struct : :class:`Structure`
            The structure to write to a Gromacs GRO file (must have coordinates)
        dest : str or file-like
            The name of a file or a file object to write the Gromacs topology to
        precision : int, optional
            The number of decimal places to print in the coordinates. Default 3
        nobox : bool, optional
            If the system does not have a periodic box defined, and this option
            is True, no box will be written. If False, the periodic box will be
            defined to enclose the solute with 0.5 nm clearance on all sides. If
            periodic box dimensions *are* defined, this variable has no effect.
        combine : 'all', None, or list of iterables, optional
            Equivalent to the combine argument of the GromacsTopologyFile.write
            method. If None, system atom order may be changed to meet the need
            for contiguously bonded groups of atoms to be part of a single
            moleculetype. All other values leave the atom order unchanged.
            Default is None.
        """

        def _write_atom_line(atom, atid, resid, has_vels, dest, precision):
            varwidth = 5 + precision
            crdfmt = f'{{:{varwidth}.{precision}f}}'
            velfmt = f'{{:{varwidth}.{precision + 1}f}}'
            dest.write(f'{resid:5d}{atom.residue.name[:5]:<5s}{atom.name[:5]:>5s}{atid:5d}')
            dest.write(crdfmt.format(atom.xx / 10)[:varwidth])
            dest.write(crdfmt.format(atom.xy / 10)[:varwidth])
            dest.write(crdfmt.format(atom.xz / 10)[:varwidth])
            if has_vels:
                dest.write(velfmt.format(atom.vx / 10)[:varwidth])
                dest.write(velfmt.format(atom.vy / 10)[:varwidth])
                dest.write(velfmt.format(atom.vz / 10)[:varwidth])
            dest.write('\n')

        own_handle = False
        if isinstance(dest, str):
            dest = genopen(dest, 'w')
            own_handle = True
        elif not hasattr(dest, 'write'):
            raise TypeError('dest must be a file name or file-like object')

        dest.write('GROningen MAchine for Chemical Simulation\n')
        dest.write(f'{len(struct.atoms):5d}\n')
        has_vels = all(hasattr(a, 'vx') for a in struct.atoms)
        if combine != 'all':
            resid, atid = 0, 0
            # use struct.split to get residue order as per topology file
            split_struct = struct.split()
            n_mols = sum(len(mol[1]) for mol in split_struct)
            unused_atoms = list(struct.atoms)
            for molid in range(n_mols):
                # loop through molids so we can get the correct molecule
                # according to the order they appear
                molecule = [mol[0] for mol in split_struct if molid in mol[1]][0]
                new_molecule = set()  # track atoms added
                last_found_atom = None  # track when gro and top diverge

                for residue in molecule.residues:
                    resid += 1
                    for atom in residue.atoms:
                        # for each atom in split topology get the first
                        # matching occurrence in the original structure
                        for original_atom in unused_atoms:
                            if atom.type == original_atom.type and \
                               atom.name == original_atom.name and \
                               atom.residue.name == original_atom.residue.name:

                                if last_found_atom is not None and \
                                   original_atom.idx != last_found_atom.idx + 1:
                                    # a rearrangement has occurred! Need to do
                                    # extra check that we've found the correct
                                    # original_atom
                                    if len(new_molecule.intersection(original_atom.bond_partners)) == 0:
                                        # original_atom must be bonded to at
                                        # least one atom in the molecule we
                                        # are currently writing otherwise find
                                        # next candidate
                                        continue

                                atid += 1
                                _write_atom_line(
                                    original_atom, atid % 100000,
                                    resid % 100000, has_vels, dest, precision
                                )
                                new_molecule.add(original_atom)
                                last_found_atom = original_atom
                                unused_atoms.remove(original_atom)
                                break
                        else:
                            raise RuntimeError(f"Could not find {atom}")
        else:
            for atom in struct.atoms:
                resid = (atom.residue.idx + 1) % 100000
                atid = (atom.idx + 1) % 100000
                _write_atom_line(atom, atid, resid, has_vels, dest, precision)
                
        # Box, in the weird format...
        if struct.box is not None:
            a, b, c = reduce_box_vectors(*box_lengths_and_angles_to_vectors(*struct.box))
            if all([abs(x-90) < TINY for x in struct.box[3:]]):
                dest.write(f'{a[0] / 10:10.5f}{b[1] / 10:10.5f}{c[2] / 10:10.5f}\n')
            else:
                dest.write(
                    f'{a[0] / 10:10.5f}{b[1] / 10:10.5f}{c[2] / 10:10.5f}{a[1] / 10:10.5f}{a[2] / 10:10.5f}'
                    f'{b[0] / 10:10.5f}{b[2] / 10:10.5f}{c[0] / 10:10.5f}{c[1] / 10:10.5f}\n'
                )
        elif not nobox and struct.atoms:
            # Find the extent of the molecule in all dimensions, and buffer it by 5 A
            crds = struct.coordinates
            diff = (crds.max(axis=0) - crds.min(axis=0)) / 10 + 0.5
            dest.write(f'{diff[0]:10.5f}{diff[1]:10.5f}{diff[2]:10.5f}\n')
        if own_handle:
            dest.close()
