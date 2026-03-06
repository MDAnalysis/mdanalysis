"""
This module contains classes for parsing and processing CHARMM parameter,
topology, and stream files. It only extracts atom properties from the
topology files and extracts all parameters from the parameter files

Author: Jason M. Swails
"""
import warnings
from collections import OrderedDict
from copy import copy
from itertools import permutations

from .exceptions import ParameterError, ParameterWarning
from .topologyobjects import (
    AtomType, DihedralType, DihedralTypeList, NoUreyBradley, UnassignedAtomType
)

class ParameterSet:
    """
    Stores a parameter set defining a force field

    Attributes
    ----------
    atom_types : dict(str:AtomType)
        Dictionary mapping the names of the atom types to the corresponding AtomType instances
    atom_types_int : dict(int:AtomType)
        Dictionary mapping the serial indexes of the atom types to the corresponding AtomType
        instances
    atom_types_tuple : dict((str,int):AtomType)
        Dictionary mapping the (name,number) tuple of the atom types to the corresponding AtomType
        instances
    bond_types : dict((str,str):AtomType)
        Dictionary mapping the 2-element tuple of the names of the two atom types involved in the
        bond to the BondType instances
    angle_types : dict((str,str,str):AngleType)
        Dictionary mapping the 3-element tuple of the names of the three atom types involved in the
        angle to the AngleType instances
    urey_bradley_types : dict((str,str,str):BondType)
        Dictionary mapping the 3-element tuple of the names of the three atom types involved in the
        angle to the BondType instances of the Urey-Bradley terms
    dihedral_types : dict((str,str,str,str):list(DihedralType))
        Dictionary mapping the 4-element tuple of the names of the four atom types involved in the
        dihedral to the DihedralType instances. Since each torsion term can be a multiterm
        expansion, each item corresponding to a key in this dict is a list of `DihedralType`s for
        each term in the expansion
    improper_types : dict((str,str,str,str):ImproperType)
        Dictionary mapping the 4-element tuple of the names of the four atom types involved in the
        improper torsion to the ImproperType instances
    improper_periodic_types : dict((str,str,str,str):DihedralType)
        Dictionary mapping the 4-element tuple of the names of the four atom types involved in the
        improper torsion (modeled as a Fourier series) to the DihedralType instances. Note, the
        central atom should always be put in the *third* position of the key
    rb_torsion_types : dict((str,str,str,str):RBTorsionType)
        Dictionary mapping the 4-element tuple of the names of the four atom types involved in the
        Ryckaert-Bellemans torsion to the RBTorsionType instances
    cmap_types : dict((str,str,str,str,str,str,str,str):CmapType)
        Dictionary mapping the 5-element tuple of the names of the five atom types involved in the
        correction map to the CmapType instances
    nbfix_types : dict((str,str):(float,float))
        Dictionary mapping the 2-element tuple of the names of the two atom types whose LJ terms are
        modified to the tuple of the (epsilon,rmin) terms for that off-diagonal term
    pair_types : dict((str,str):NonbondedExceptionType)
        Dictionary mapping the 2-element tuple of atom type names for which explicit exclusion rules
        should be applied
    nbthole_types : dict(())
    parametersets : list(str)
        List of parameter set names processed in the current ParameterSet
    residues : dict(str:ResidueTemplate|ResidueTemplateContainer)
        A library of ResidueTemplate objects mapped to the residue name defined in the force field
        library files
    """

    def __init__(self):
        # Instantiate the list types
        self.atom_types = self.atom_types_str = OrderedDict()
        self.atom_types_int = OrderedDict()
        self.atom_types_tuple = OrderedDict()
        self.bond_types = OrderedDict()
        self.angle_types = OrderedDict()
        self.urey_bradley_types = OrderedDict()
        self.dihedral_types = OrderedDict()
        self.improper_types = OrderedDict()
        self.improper_periodic_types = OrderedDict()
        self.rb_torsion_types = OrderedDict()
        self.cmap_types = OrderedDict()
        self.nbfix_types = OrderedDict()
        self.pair_types = OrderedDict()
        self.nbthole_types = OrderedDict()
        self.parametersets = []
        self._combining_rule = 'lorentz'
        self.residues = OrderedDict()
        self.patches = OrderedDict()
        self.default_scee = self.default_scnb = 1.0
        self._improper_key_map = OrderedDict()

    def __copy__(self):
        other = type(self)()
        for key, item in self.atom_types.items():
            other.atom_types[key] = copy(item)
        for key, item in self.atom_types_int.items():
            other.atom_types_int[key] = copy(item)
        for key, item in self.atom_types_tuple.items():
            other.atom_types_tuple[key] = copy(item)
        for key, item in self.bond_types.items():
            if key in other.bond_types: continue
            typ = copy(item)
            other.bond_types[key] = typ
            other.bond_types[tuple(reversed(key))] = typ
        for key, item in self.pair_types.items():
            if key in other.pair_types: continue
            typ = copy(item)
            other.pair_types[key] = typ
            other.pair_types[tuple(reversed(key))] = typ
        for key, item in self.angle_types.items():
            if key in other.angle_types: continue
            typ = copy(item)
            other.angle_types[key] = typ
            other.angle_types[tuple(reversed(key))] = typ
        for key, item in self.dihedral_types.items():
            if key in other.dihedral_types: continue
            typ = copy(item)
            other.dihedral_types[key] = typ
            other.dihedral_types[tuple(reversed(key))] = typ
        for key, item in self.rb_torsion_types.items():
            if key in other.rb_torsion_types: continue
            typ = copy(item)
            other.rb_torsion_types[key] = typ
            other.rb_torsion_types[tuple(reversed(key))] = typ
        for key, item in self.improper_types.items():
            if key in other.improper_types: continue
            other.improper_types[key] = copy(item)
        for key, item in self.improper_periodic_types.items():
            if key in other.improper_periodic_types: continue
            typ = copy(item)
            other.improper_periodic_types[key] = typ
            other.improper_periodic_types[tuple(reversed(key))] = typ
        for key, item in self.urey_bradley_types.items():
            if key in other.urey_bradley_types: continue
            typ = copy(item)
            other.urey_bradley_types[key] = typ
            other.urey_bradley_types[tuple(reversed(key))] = typ
        for key, item in self.cmap_types.items():
            if key in other.cmap_types: continue
            typ = copy(item)
            other.cmap_types[key] = typ
            other.cmap_types[tuple(reversed(key))] = typ
        for key, item in self.residues.items():
            other.residues[key] = copy(item)
        for key, item in self.patches.items():
            other.patches[key] = copy(item)
        other._improper_key_map = copy(self._improper_key_map)
        other.combining_rule = self.combining_rule

        return other

    @classmethod
    def from_structure(cls, struct, allow_unequal_duplicates=True):
        """ Extracts known parameters from a Structure instance

        Parameters
        ----------
        struct : :class:`parmed.structure.Structure`
            The parametrized ``Structure`` instance from which to extract
            parameters into a ParameterSet
        allow_unequal_duplicates : bool, optional
            If True, if two or more unequal parameter types are defined by the
            same atom types, the last one encountered will be assigned. If
            False, an exception will be raised. Default is True

        Returns
        -------
        params : :class:`ParameterSet`
            The parameter set with all parameters defined in the Structure

        Notes
        -----
        The parameters here are copies of the ones in the Structure, so
        modifying the generated ParameterSet will have no effect on ``struct``.
        Furthermore, the *first* occurrence of each parameter will be used. If
        future ones differ, they will be silently ignored, since this is
        expected behavior in some instances (like with Gromacs topologies in the
        ff99sb-ildn force field) unless ``allow_unequal_duplicates`` is set to
        ``False``

        Dihedrals are a little trickier. They can be multi-term, which can be
        represented either as a *single* entry in dihedrals with a type of
        DihedralTypeList or multiple entries in dihedrals with a DihedralType
        parameter type. In this case, the parameter is constructed from either
        the first DihedralTypeList found or the first DihedralType of each
        periodicity found if no matching DihedralTypeList is found.

        Raises
        ------
        :class:`parmed.exceptions.ParameterError` if allow_unequal_duplicates is
        False and 2+ unequal parameters are defined between the same atom types.

        `NotImplementedError` if any AMOEBA potential terms are defined in the
        input structure
        """
        params = cls()
        found_dihed_type_list = dict()
        for atom in struct.atoms:
            if atom.atom_type in (UnassignedAtomType, None):
                atom_type = AtomType(atom.type, None, atom.mass, atom.atomic_number)
                atom_type.set_lj_params(atom.epsilon, atom.rmin, atom.epsilon_14, atom.rmin_14)
                params.atom_types[atom.type] = atom_type
            else:
                atom_type = copy(atom.atom_type)
                params.atom_types[str(atom_type)] = atom_type
                if atom_type.number is not None:
                    params.atom_types_int[int(atom_type)] = atom_type
                    params.atom_types_tuple[(int(atom_type), str(atom_type))] = atom_type
        if struct.has_NBFIX():
            for atom in struct.atoms:
                if atom.atom_type.nbfix:
                    other_atoms = list(atom.atom_type.nbfix.keys())
                    for other_atom in other_atoms:
                        (rmin, epsilon, rmin14, epsilon14) = atom.atom_type.nbfix[other_atom]
                        if (other_atom, atom.type) in params.nbfix_types:
                            continue
                        params.nbfix_types[(atom.type, other_atom)] = (epsilon, rmin)
        for bond in struct.bonds:
            if bond.type is None: continue
            key = (bond.atom1.type, bond.atom2.type)
            if key in params.bond_types:
                if not allow_unequal_duplicates and params.bond_types[key] != bond.type:
                    raise ParameterError(f'Unequal bond types defined between {key[0]} and {key[1]}')
                continue # pragma: no cover
            typ = copy(bond.type)
            key = (bond.atom1.type, bond.atom2.type)
            params.bond_types[key] = typ
            params.bond_types[tuple(reversed(key))] = typ
        for angle in struct.angles:
            if angle.type is None: continue
            key = (angle.atom1.type, angle.atom2.type, angle.atom3.type)
            if key in params.angle_types:
                if not allow_unequal_duplicates and params.angle_types[key] != angle.type:
                    raise ParameterError('Unequal angle types defined between {}, {}, and {}'.format(*key))
                continue # pragma: no cover
            typ = copy(angle.type)
            key = (angle.atom1.type, angle.atom2.type, angle.atom3.type)
            params.angle_types[key] = typ
            params.angle_types[tuple(reversed(key))] = typ
            if angle.funct == 5:
                key = (angle.atom1.type, angle.atom2.type, angle.atom3.type)
                params.urey_bradley_types[key] = NoUreyBradley
                params.urey_bradley_types[tuple(reversed(key))] = NoUreyBradley
        for dihedral in struct.dihedrals:
            if dihedral.type is None: continue
            if dihedral.improper:
                # Make sure the central atom comes third, but add all permutations since it's not
                # clear which particular ordering the various programs require
                for key in _find_improper_keys(dihedral):
                    if key in params.improper_periodic_types:
                        if (not allow_unequal_duplicates and
                                params.improper_periodic_types[key] != dihedral.type):
                            raise ParameterError(
                                'Unequal dihedral types defined between {}, {}, {}, and {}'.format(*key)
                            )
                        continue # pragma: no cover
                    typ = copy(dihedral.type)
                    params.improper_periodic_types[key] = typ
            else:
                key = (dihedral.atom1.type, dihedral.atom2.type,
                       dihedral.atom3.type, dihedral.atom4.type)
                # Proper dihedral. Look out for multi-term forms
                if key in params.dihedral_types and found_dihed_type_list[key]:
                    # Already found a multi-term dihedral type list
                    if not allow_unequal_duplicates:
                        if isinstance(dihedral.type, DihedralTypeList):
                            if params.dihedral_types[key] != dihedral.type:
                                raise ParameterError(
                                    'Unequal dihedral types defined between {}, {}, {}, and {}'.format(*key)
                                )
                        elif isinstance(dihedral.type, DihedralType):
                            for dt in params.dihedral_types[key]:
                                if dt == dihedral.type:
                                    break
                            else:
                                raise ParameterError('Unequal dihedral types defined between %s, '
                                                     '%s, %s, and %s' % key)
                    continue # pragma: no cover
                elif key in params.dihedral_types:
                    # We have one term of a potentially multi-term dihedral.
                    if isinstance(dihedral.type, DihedralTypeList):
                        # This is a full Fourier series list
                        found_dihed_type_list[key] = True
                        found_dihed_type_list[tuple(reversed(key))] = True
                        typ = copy(dihedral.type)
                        params.dihedral_types[key] = typ
                        params.dihedral_types[tuple(reversed(key))] = typ
                    else:
                        # This *might* be another term. Make sure another term
                        # with its periodicity does not already exist
                        for t in params.dihedral_types[key]:
                            if t.per == dihedral.type.per:
                                if not allow_unequal_duplicates and t != dihedral.type:
                                    raise ParameterError(
                                        'Unequal dihedral types defined bewteen {}, {}, {}, and {}'.format(*key)
                                    )
                                break
                        else:
                            # If we got here, we did NOT find this periodicity.
                            # And since this is mutating a list in-place, it
                            # automatically propagates to the reversed key
                            typ = copy(dihedral.type)
                            params.dihedral_types[key].append(typ)
                else:
                    # New parameter. If it's a DihedralTypeList, assign it and
                    # be done with it. If it's a DihedralType, start a
                    # DihedralTypeList to be added to later.
                    if isinstance(dihedral.type, DihedralTypeList):
                        found_dihed_type_list[key] = True
                        found_dihed_type_list[tuple(reversed(key))] = True
                        typ = copy(dihedral.type)
                        params.dihedral_types[key] = typ
                        params.dihedral_types[tuple(reversed(key))] = typ
                    else:
                        found_dihed_type_list[key] = False
                        found_dihed_type_list[tuple(reversed(key))] = False
                        typ = DihedralTypeList()
                        typ.append(copy(dihedral.type))
                        params.dihedral_types[key] = typ
                        params.dihedral_types[tuple(reversed(key))] = typ
        for improper in struct.impropers:
            if improper.type is None: continue
            key = (improper.atom1.type, improper.atom2.type,
                   improper.atom3.type, improper.atom4.type)
            if key in params.improper_types:
                if not allow_unequal_duplicates and params.improper_types[key] != improper.type:
                    raise ParameterError('Unequal improper types defined between {}, {}, {}, and {}'.format(*key))
                continue # pragma: no cover
            params.improper_types[key] = copy(improper.type)
        for cmap in struct.cmaps:
            if cmap.type is None: continue
            key = (cmap.atom1.type, cmap.atom2.type, cmap.atom3.type, cmap.atom4.type,
                   cmap.atom2.type, cmap.atom3.type, cmap.atom4.type, cmap.atom5.type)
            if key in params.cmap_types:
                if not allow_unequal_duplicates and cmap.type != params.cmap_types[key]:
                    raise ParameterError(
                        f'Unequal CMAP types defined between {key[0]}, {key[1]}, {key[2]}, {key[3]}, and {key[7]}'
                    )
                continue # pragma: no cover
            typ = copy(cmap.type)
            params.cmap_types[key] = typ
            params.cmap_types[tuple(reversed(key))] = typ
        urey_brads_preassigned = len(params.urey_bradley_types) > 0
        for urey in struct.urey_bradleys:
            if urey.type is None or urey.type is NoUreyBradley: continue
            key = _find_ureybrad_key(urey)
            if key is None: continue
            if urey_brads_preassigned and key not in params.urey_bradley_types:
                warnings.warn('Angle corresponding to Urey-Bradley type not found', ParameterWarning)
            typ = copy(urey.type)
            params.urey_bradley_types[key] = typ
            params.urey_bradley_types[tuple(reversed(key))] = typ
        if not urey_brads_preassigned and len(params.urey_bradley_types) > 0:
            # Go through all of our angle parameters and make sure there is a
            # matching Urey-Bradley list. If there's not, that means there is no
            # Urey-Bradley term for that angle
            for key in params.angle_types:
                if key in params.urey_bradley_types: continue
                params.urey_bradley_types[key] = NoUreyBradley
        for adjust in struct.adjusts:
            if adjust.type is None: continue
            key = (adjust.atom1.type, adjust.atom2.type)
            if key in params.pair_types:
                if not allow_unequal_duplicates and params.pair_types[key] != adjust.type:
                    raise ParameterError('Unequal pair types defined between {} and {}'.format(*key))
                continue # pragma: no cover
            typ = copy(adjust.type)
            params.pair_types[key] = typ
            params.pair_types[tuple(reversed(key))] = typ
        # Trap for Amoeba potentials
        if (struct.trigonal_angles or struct.out_of_plane_bends or struct.torsion_torsions or
                struct.stretch_bends or struct.trigonal_angles or struct.pi_torsions):
            raise NotImplementedError('Cannot extract parameters from an Amoeba-parametrized system')
        return params

    def condense(self, do_dihedrals=True):
        """
        This function goes through each of the parameter type dicts and
        eliminates duplicate types. After calling this function, every unique
        bond, angle, dihedral, improper, or cmap type will pair with EVERY key
        in the type mapping dictionaries that points to the equivalent type

        Parameters
        ----------
        do_dihedrals : bool=True
            Dihedrals can take the longest time to compress since testing their
            equality takes the longest (this is complicated by the existence of
            multi-term torsions). This flag will allow you to *skip* condensing
            the dihedral parameter types (for large parameter sets, this can cut
            the compression time in half)

        Returns
        -------
        self
            The instance that is being condensed

        Notes
        -----
        The return value allows you to condense the types at construction time.

        Example
        -------
        >>> params = ParameterSet().condense()
        >>> params
        <parmed.parameters.ParameterSet at 0x7f88757de090>
        """
        # First scan through all of the bond types
        self._condense_types(self.bond_types)
        self._condense_types(self.angle_types)
        self._condense_types(self.urey_bradley_types)
        if do_dihedrals:
            self._condense_types(self.dihedral_types)
            self._condense_types(self.rb_torsion_types)
        self._condense_types(self.improper_periodic_types)
        self._condense_types(self.improper_types)
        self._condense_types(self.cmap_types)
        return self

    @staticmethod
    def _condense_types(typedict):
        """
        Loops through the given dict and condenses all types.

        Parameter
        ---------
        typedict : dict
            Type dictionary to condense
        """
        keylist = list(typedict.keys())
        for i in range(len(keylist) - 1):
            key1 = keylist[i]
            for j in range(i+1, len(keylist)):
                key2 = keylist[j]
                if typedict[key1] == typedict[key2]:
                    typedict[key2] = typedict[key1]

    @property
    def combining_rule(self):
        return self._combining_rule

    @combining_rule.setter
    def combining_rule(self, value):
        if value not in ('lorentz', 'geometric'):
            raise ValueError('combining_rule must be "lorentz" or "geometric"')
        self._combining_rule = value

    def typeify_templates(self):
        """ Assign atom types to atom names in templates """
        from .modeller import ResidueTemplateContainer, ResidueTemplate
        for residue in self.residues.values():
            if isinstance(residue, ResidueTemplateContainer):
                for res in residue:
                    for atom in res:
                        atom.atom_type = self.atom_types[atom.type]
            else:
                assert isinstance(residue, ResidueTemplate), 'Wrong type!'
                for atom in residue:
                    atom.atom_type = self.atom_types[atom.type]

def _find_ureybrad_key(urey):
    """
    Finds a key for a given Urey-Bradley by finding the middle atom in an angle.
    Raises a ParameterWarning if no middle atom found
    """
    a1, a2 = urey.atom1, urey.atom2
    shared_bond_partners = set(a1.bond_partners) & set(a2.bond_partners)
    if len({a.type for a in shared_bond_partners}) != 1:
        warnings.warn(f'Urey-Bradley {repr(urey)} shares multiple central atoms', ParameterWarning)
    return (a1.type, list(shared_bond_partners)[0].type, a2.type)

def _find_improper_keys(dih):
    """ Finds the central atom (i.e., that bonded to everything else) """
    assert dih.improper, 'Should not be called on non-improper!'
    all_atoms = {dih.atom1, dih.atom2, dih.atom3, dih.atom4}
    for atom in all_atoms:
        for oatom in all_atoms:
            if oatom is atom:
                continue
            if oatom not in atom.bond_partners:
                break
        else:
            # This *is* the central atom
            for key in permutations([a.type for a in all_atoms if a is not atom]):
                yield (key[0], key[1], atom.type, key[2])
            return # break out of the generator
    # If we got here, we found no central atom. *assume* it's the third spot already...
    for key in permutations([dih.atom1.type, dih.atom2.type, dih.atom4.type]):
        yield (key[0], key[1], dih.atom3.type, key[2])
