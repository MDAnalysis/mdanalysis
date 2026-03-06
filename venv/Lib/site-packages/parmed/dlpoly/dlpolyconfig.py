"""
This module contains functionality relevant to loading and parsing DLPOLY CONFIG
(coordinate) files and building a stripped-down Structure from it
"""
from ..constants import TINY
from ..formats.registry import FileFormatType
from ..geometry import box_lengths_and_angles_to_vectors, reduce_box_vectors
from ..utils.io import genopen

class DlpolyConfigFile(metaclass=FileFormatType):
    """ Writes Dlpoly CONFIG files """
    #===================================================

    @staticmethod
    def write(struct, dest, precision=12, combine=False):
        """ Write a Dlpoly CONFIG File from a Structure

        Parameters
        ----------
        struct : :class:`Structure`
            The structure to write to a Dlpoly CONFIG file (must have coordinates)
        dest : str or file-like
            The name of a file or a file object to write the Dlpoly CONFIG file to
        precision : int, optional
            The number of decimal places to print in the coordinates. Default 9
        combine : 'all', None, or list of iterables, optional
            Equivalent to the combine argument of the DlpolyFieldFile.write
            method. If None, system atom order may be changed to meet the need
            for contiguously bonded groups of atoms to be part of a single
            moleculetype. All other values leave the atom order unchanged.
            Default is None.
        """

        def _write_atom_line(atom, atid, resid, has_vels, dest, precision):
            varwidth = 8 + precision
            crdfmt = '%%%d.%df' % (varwidth, precision)
            velfmt = '%%%d.%df' % (varwidth, precision)
            dest.write('%-8s %6d' % (atom.type[:5], atid))
            dest.write('\n')
            dest.write((crdfmt % (atom.xx))[:varwidth])
            dest.write((crdfmt % (atom.xy))[:varwidth])
            dest.write((crdfmt % (atom.xz))[:varwidth])
            dest.write('\n')
            if has_vels:
                dest.write((velfmt % (atom.vx))[:varwidth])
                dest.write((velfmt % (atom.vy))[:varwidth])
                dest.write((velfmt % (atom.vz))[:varwidth])
                dest.write('\n')

        own_handle = False
        if isinstance(dest, str):
            dest = genopen(dest, 'w')
            own_handle = True
        elif not hasattr(dest, 'write'):
            raise TypeError('dest must be a file name or file-like object')

        dest.write('CONFIGURATION\n')
        # CONFIG file key (do we have velocities as well?)
        has_vels = all(hasattr(a, 'vx') for a in struct.atoms)
        if has_vels:
            dest.write('         1')
        else:
            dest.write('         0')
        # Periodic Boundary Key
        if struct.box is not None:
            a, b, c = reduce_box_vectors(*box_lengths_and_angles_to_vectors(*struct.box))
            if all([abs(x-90) < TINY for x in struct.box[3:]]):
                if (a[0] == b[1] and b[1] == c[2]):
                    dest.write('         1')
                else:
                    dest.write('         2')
            else:
                dest.write('         3')
        else:
            dest.write('         0')
        dest.write('     %5d         0.0000000000E+00\n' % len(struct.atoms))
        # Box information
        if struct.box is not None:
            varwidth = 8 + precision
            boxfmt = '%%%d.%df%%%d.%df%%%d.%df\n' % (varwidth, precision, varwidth, precision, varwidth, precision)
            dest.write(boxfmt % (a[0], a[1], a[2]))
            dest.write(boxfmt % (b[0], b[1], b[2]))
            dest.write(boxfmt % (c[0], c[1], c[2]))
        if combine != 'all':
            resid, atid = 0, 0
            # use struct.split to get residue order as per topology file
            split_struct = struct.split()
            n_mols = sum(len(mol[1]) for mol in split_struct)
            unused_atoms = list(struct.atoms)
            for molid in range(n_mols):
                # loop through molids so we can get the correct molecule
                # according to the order they appear
                molecule = [
                    mol[0] for mol in split_struct if molid in mol[1]][0]
                new_molecule = set()  # track atoms added
                last_found_atom = None  # track when CONFIG and FIELD diverge

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
                                    if len(new_molecule.intersection(
                                            original_atom.bond_partners)) == 0:
                                        # original_atom must be bonded to at
                                        # least one atom in the molecule we
                                        # are currently writing otherwise find
                                        # next candidate
                                        continue

                                atid += 1
                                _write_atom_line(
                                    original_atom, atid % 100000,
                                    resid % 100000, has_vels, dest, precision)
                                new_molecule.add(original_atom)
                                last_found_atom = original_atom
                                unused_atoms.remove(original_atom)
                                break
                        else:
                            raise RuntimeError("Could not find %s" % atom)
        else:
            for atom in struct.atoms:
                resid = (atom.residue.idx + 1) % 100000
                atid = (atom.idx + 1) % 100000
                _write_atom_line(atom, atid, resid, has_vels, dest, precision)

        if own_handle:
            dest.close()
