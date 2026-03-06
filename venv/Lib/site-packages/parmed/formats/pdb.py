"""
This package contains classes responsible for reading and writing both PDB and
PDBx/mmCIF files.
"""
from collections import OrderedDict, namedtuple, Counter
from contextlib import closing
from functools import lru_cache, reduce
from string import ascii_letters
import io
import urllib
import gzip
import logging
import urllib.error
import urllib.request
import numpy as np
from ..exceptions import PDBError, PDBWarning
from ..formats.pdbx import PdbxReader, PdbxWriter, containers
from ..formats.registry import FileFormatType
from ..periodic_table import AtomicNum, Mass, Element, element_by_name
from ..residue import AminoAcidResidue, RNAResidue, DNAResidue, WATER_NAMES
from ..modeller import StandardBiomolecularResidues, get_standard_residue_template_library, get_nonstandard_ccd_residues
from ..structure import Structure
from ..topologyobjects import Atom, ExtraPoint, Bond, Link, QualitativeBondType
from ..symmetry import Symmetry
from ..utils.io import genopen
import re
import warnings

LOGGER = logging.getLogger(__name__)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_ascii_letters_set = set(ascii_letters)

def _compare_atoms(old_atom, new_atom, resname, resid, chain, segid, inscode):
    """
    Compares two atom instances, along with the residue name, number, and chain
    identifier, to determine if two atoms are actually the *same* atom, but
    simply different conformations

    Parameters
    ----------
    old_atom : :class:`Atom`
        The original atom that has been added to the structure already
    new_atom : :class:`Atom`
        The new atom that we want to see if it is the same as the old atom
    resname : ``str``
        The name of the residue that the new atom would belong to
    resid : ``int``
        The number of the residue that the new atom would belong to
    chain : ``str``
        The chain identifier that the new atom would belong to
    segid : ``str``
        The segment identifier for the molecule
    inscode : ``str``
        The insertion code for the residue

    Returns
    -------
    True if they are the same atom, False otherwise
    """
    if old_atom.name != new_atom.name: return False
    if old_atom.residue.name != resname: return False
    if old_atom.residue.number != resid: return False
    if old_atom.residue.chain != chain.strip(): return False
    if old_atom.residue.segid != segid.strip(): return False
    if old_atom.residue.insertion_code != inscode.strip(): return False
    return True

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _standardize_resname(resname):
    """ Looks up a standardized residue name for the given resname """
    try:
        return AminoAcidResidue.get(resname, abbronly=True).abbr, False
    except KeyError:
        try:
            return RNAResidue.get(resname).abbr, False
        except KeyError:
            try:
                return DNAResidue.get(resname).abbr, False
            except KeyError:
                if resname.strip() in WATER_NAMES:
                    return 'HOH', True
                else:
                    return resname, True

def _is_hetatm(resname):
    """ Sees if residue name is "standard", otherwise, we need to use HETATM to
    print atom records instead of ATOM
    """
    if len(resname) != 3:
        return not (RNAResidue.has(resname) or DNAResidue.has(resname))
    return not (AminoAcidResidue.has(resname) or RNAResidue.has(resname) or DNAResidue.has(resname))

def _number_truncated_to_n_digits(num, digits):
    """ Truncates the given number to the specified number of digits """
    if num < 0:
        return int(-(-num % eval(f'1e{digits - 1:d}')))
    return int(num % eval(f'1e{digits:d}'))

_ATOM_NAME_TO_ATOM_CACHE = {}
def _atom_name_to_atom_map(residue_name, *reslibs):
    if residue_name in _ATOM_NAME_TO_ATOM_CACHE:
        return _ATOM_NAME_TO_ATOM_CACHE[residue_name]
    for reslib in reslibs:
        if residue_name in reslib:
            _ATOM_NAME_TO_ATOM_CACHE[residue_name] = {atom.name: atom for atom in reslib[residue_name].atoms}
    return _ATOM_NAME_TO_ATOM_CACHE.get(residue_name, {})

_ATOM_NAME_TO_BOND_CACHE = {}
def _atom_name_to_bond_map(residue_name, templates):
    if residue_name in _ATOM_NAME_TO_BOND_CACHE:
        return _ATOM_NAME_TO_BOND_CACHE[residue_name]
    if residue_name not in templates:
        return {}
    bond_type_map = {}
    for bond in templates[residue_name].bonds:
        bond_type_map[(bond.atom1.name, bond.atom2.name)] = bond.qualitative_type
        bond_type_map[(bond.atom2.name, bond.atom1.name)] = bond.qualitative_type
    _ATOM_NAME_TO_BOND_CACHE[residue_name] = bond_type_map
    return bond_type_map

def _replace_if_none(obj, template_obj, attribute):
    if getattr(obj, attribute) is None:
        setattr(obj, attribute, getattr(template_obj, attribute))

def _map_atomic_properties_from_templates(structure: Structure, *reslibs):
    """Modifies atoms in a structure in-place with attributes from the template"""
    for residue in structure.residues:
        template_atom_map = _atom_name_to_atom_map(residue.name, *reslibs)
        for atom in residue.atoms:
            if atom.name not in template_atom_map:
                continue
            template_atom = template_atom_map[atom.name]
            _replace_if_none(atom, template_atom, "formal_charge")
            _replace_if_none(atom, template_atom, "aromatic")
            _replace_if_none(atom, template_atom, "hybridization")

def _map_bond_types_from_templates(structure: Structure, *reslibs):
    templates = reduce(lambda d1, d2: {**d2, **d1}, reslibs, {})
    for bond in structure.bonds:
        # Only assign bond types of bonds we recognize
        if bond.atom1.residue.name not in templates or bond.atom2.residue.name not in templates:
            continue
        # Don't overwrite something already there
        if bond.qualitative_type is not None:
            continue
        # For RNA, DNA, and peptides all inter-residue bonds are single bonds.
        if bond.atom1.residue is not bond.atom2.residue:
            bond.qualitative_type = QualitativeBondType.SINGLE
            continue
        # Otherwise they're both in the same residue, and should be in the template
        bond_type_map = _atom_name_to_bond_map(bond.atom1.residue.name, templates)
        qualitative_type = bond_type_map.get((bond.atom1.name, bond.atom2.name), None)
        bond.qualitative_type = qualitative_type

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PDBFile(metaclass=FileFormatType):
    """ Standard PDB file format parser and writer """
    #===================================================

    AtomLookupKey = namedtuple(
        'AtomLookupKey', ('name', 'number', 'residue_name', 'residue_number', 'chain',
                          'insertion_code', 'segment_id', 'alternate_location')
    )

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as a PDB file

        Parameters
        ----------
        filename : str or file object
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is a PDB file
        """
        if isinstance(filename, str):
            own_handle = True
            fileobject = genopen(filename, 'r')
        elif hasattr(filename, 'read'):
            own_handle = False
            fileobject = filename

        try:
            for line in fileobject:
                if line[:6] in {'CRYST1', 'END   ', 'END', 'HEADER', 'NUMMDL', 'MASTER', 'AUTHOR', 'CAVEAT',
                        'COMPND', 'EXPDTA', 'MDLTYP', 'KEYWDS', 'OBSLTE', 'SOURCE', 'SPLIT ', 'SPRSDE', 'TITLE ',
                        'ANISOU', 'CISPEP', 'CONECT', 'DBREF ', 'HELIX ', 'HET   ', 'LINK  ', 'MODRES',
                        'REVDAT', 'SEQADV', 'SHEET ', 'SSBOND', 'FORMUL', 'HETNAM', 'HETSYN', 'SEQRES', 'SITE  ',
                        'ENDMDL', 'MODEL ', 'TER   ', 'JRNL  ', 'REMARK', 'TER', 'DBREF ', 'DBREF2', 'DBREF1',
                        'DBREF', 'HET', 'LINKR ', 'TURN  '}:
                    continue
                # Hack to support reduce-added flags
                elif line[:6] == 'USER  ' and line[6:9] == 'MOD':
                    continue
                elif line[:5] in ('ORIGX', 'SCALE', 'MTRIX'):
                    if line[5] not in '123':
                        return False
                elif line[:6] in ('ATOM  ', 'HETATM'):
                    atnum, atname = line[6:11], line[12:16]
                    resname, resid = line[17:20], line[22:26]
                    x, y, z = line[30:38], line[38:46], line[46:54]
                    occupancy, bfactor = line[54:60], line[60:66]
                    elem = line[76:78]
                    # Check for various attributes. This is the first atom, so
                    # we can assume we haven't gotten into the regime of "weird"
                    # yet, like hexadecimal atom/residue indices.
                    if not atnum.strip().lstrip('-').isdigit(): return False
                    if atname.strip().isdigit(): return False
                    if not resname.strip(): return False
                    if not resid.strip().lstrip('-').isdigit(): return False
                    try:
                        float(x), float(y), float(z)
                    except ValueError:
                        return False
                    if occupancy.strip():
                        try:
                            float(occupancy)
                        except ValueError:
                            return False
                    if bfactor.strip():
                        try:
                            float(bfactor)
                        except ValueError:
                            return False
                    if elem.strip():
                        if any(x.isdigit() for x in elem):
                            return False
                    return True
                else:
                    return False
            return False
        finally:
            if own_handle:
                fileobject.close()

    #===================================================

    def __init__(self, fileobj):
        # Open file object that we are parsing
        self.fileobj = fileobj

        self._atom_map_from_attributes = OrderedDict()
        self._atom_map_from_all_attributes = OrderedDict()
        self._current_model_number = 1
        self._residue_indices_overflow = False
        self._atom_indices_overflow = False
        self.struct = Structure()
        self._symmetry_lines = []
        self._link_lines = []
        self._coordinates = [[]]
        self._model_atom_counts = [0]
        self._anisou_records = dict()
        self._atom_map_from_atom_number = dict()
        self._atom_map_to_parent = dict()
        self._model1_atoms_in_structure = set()
        self._model_open = True
        self._insertion_codes = set() # Only permitted for compliant PDB files
        self._last_atom = None
        self._last_residue = None
        # Some writers use additional fields to hold extra digits for the residue number if it goes
        # more than 4 digits
        self._residue_number_field_extended_by = 0
        # Fallback if parser fails to process residue numbers
        self._last_residue_number_label = None

    #===================================================

    @staticmethod
    def download(pdb_id, timeout=10, saveto=None):
        """
        Goes to the wwPDB website and downloads the requested PDB, loading it
        as a :class:`Structure` instance

        Parameters
        ----------
        pdb_id : str
            The 4-letter PDB ID to try and download from the RCSB PDB database
        timeout : float, optional
            The number of seconds to wait before raising a timeout error.
            Default is 10 seconds
        saveto : str, optional
            If provided, this will be treated as a file name to which the PDB
            file will be saved. If None (default), no PDB file will be written.
            This will be a verbatim copy of the downloaded PDB file, unlike the
            somewhat-stripped version you would get by using
            :meth:`Structure.write_pdb <parmed.structure.Structure.write_pdb>`

        Returns
        -------
        struct : :class:`Structure <parmed.structure.Structure>`
            Structure instance populated by the requested PDB

        Raises
        ------
        socket.timeout if the connection times out while trying to contact the
        FTP server

        IOError if there is a problem retrieving the requested PDB or writing a
        requested ``saveto`` file

        ImportError if the gzip module is not available

        TypeError if pdb_id is not a 4-character string
        """
        if not isinstance(pdb_id, str) or len(pdb_id) != 4:
            raise ValueError('pdb_id must be the 4-letter PDB code')

        pdb_id = pdb_id.lower()
        
        url = (
            'https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/%s/pdb%s.ent.gz'
        ) % (pdb_id[1:3], pdb_id)

        fileobj = io.BytesIO()
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                fileobj.write(response.read())
        except (urllib.error.URLError, urllib.error.HTTPError) as err:
            raise IOError(f"Could not retrieve PDB ID {pdb_id}; {err}") from err

        # Rewind, wrap it in a GzipFile and send it to parse
        fileobj.seek(0)
        fileobj = io.TextIOWrapper(gzip.GzipFile(fileobj=fileobj, mode='r'))
        if saveto is not None:
            with closing(genopen(saveto, 'w')) as f:
                f.write(fileobj.read())
            fileobj.seek(0)
        return PDBFile.parse(fileobj)

    #===================================================

    _relatere = re.compile(r'RELATED ID: *(\w+) *RELATED DB: *(\w+)', re.I)

    def _add_metadata_fields(self):
        self.struct.experimental = ''
        self.struct.journal = ''
        self.struct.authors = ''
        self.struct.keywords = ''
        self.struct.doi = ''
        self.struct.pmid = ''
        self.struct.journal_authors = ''
        self.struct.volume = ''
        self.struct.title = ''
        self.struct.year = None
        self.struct.resolution = None
        self.struct.related_entries = []

    @classmethod
    def parse(cls, filename, skip_bonds=False, expanded_residue_template_match=False, all_residue_template_match=False):
        """ Read a PDB file and return a populated `Structure` class

        Parameters
        ----------
        filename : str or file-like
            Name of the PDB file to read, or a file-like object that can iterate
            over the lines of a PDB. Compressed file names can be specified and
            are determined by file-name extension (e.g., file.pdb.gz,
            file.pdb.bz2)
        skip_bonds : bool, optional
            If True, skip trying to assign bonds. This can save substantial time
            when parsing large files with non-standard residue names. However,
            no bonds are assigned. This is OK if, for instance, the PDB file is
            being parsed simply for its coordinates. This may also reduce
            element assignment if element information is not present in the PDB
            file already. Default is False.
        expanded_residue_template_match : bool, optional
            If True, an expanded set of amino and nucleic acid residue templates taken from the
            chemical component database is used to assign bonds and atomic properties.
        all_residue_template_match : bool, optional
            If True, all residue templates taken from the chemical component database is used
            to assign bonds and atomic properties. Note that loading the library for all known
            3-letter residues is slow and memory-intensive. The library is cached after the
            first use.

        Metadata
        --------
        The PDB parser also adds metadata to the returned Structure object that
        may be present in the PDB file

        experimental : ``str``
            EXPDTA record
        journal : ``str``
            JRNL record
        authors : ``str``
            AUTHOR records
        keywords : ``str``
            KEYWDS records
        doi : ``str``
            DOI from the JRNL record
        pmid : ``str``
            PMID from the JRNL record
        journal_authors : ``str``
            Author info from the JRNL record
        volume : ``str``
            Volume of the published article from the JRNL record
        page : ``str``
            Page of the published article from the JRNL record
        title : ``str``
            TITL section of the JRNL record
        year : ``int``
            Year that the article was published, from the JRNL record
        resolution : ``float``
            The X-RAY resolution in Angstroms, or None if not found
        related_entries : ``list of (str, str)``
            List of entries in other databases

        Returns
        -------
        structure : :class:`Structure`
            The Structure object initialized with all the information from the PDB
            file.  No bonds or other topological features are added by default.
        """
        if all_residue_template_match and not expanded_residue_template_match:
            raise ValueError("all_residue_template_match cannot be True if expanded_residue_template_match is False")
        if expanded_residue_template_match and skip_bonds:
            raise ValueError("extended_residue_template_match cannot be True if skip_bonds is also True")
        if isinstance(filename, str):
            own_handle = True
            fileobj = genopen(filename, 'r')
        else:
            own_handle = False
            fileobj = filename

        inst = cls(fileobj)
        inst._add_metadata_fields()

        # Support hexadecimal numbering like that printed by VMD
        try:
            inst._parse_open_file(fileobj)
        finally:
            if own_handle:
                fileobj.close()

        # Assign bonds based on standard templates and simple distances
        if not skip_bonds:
            if expanded_residue_template_match:
                residue_libraries = [get_standard_residue_template_library()]
                if all_residue_template_match:
                    residue_libraries.append(get_nonstandard_ccd_residues())
                inst.struct.assign_bonds(*residue_libraries)
                _map_atomic_properties_from_templates(inst.struct, *residue_libraries)
                _map_bond_types_from_templates(inst.struct, *residue_libraries)
            else:
                inst.struct.assign_bonds()

        inst._postprocess_metadata()
        inst.struct.unchange()

        try:
            inst.struct.coordinates = inst._coordinates
        except ValueError:
            raise PDBError(
                'Coordinate shape mismatch. Probably caused by different atom counts in some of the PDB models'
            )

        # Postprocess features of the PDB file we couldn't resolve until the end
        inst._process_structure_symmetry()
        inst._process_link_records()
        inst._assign_anisou_to_atoms()
        if inst._residue_indices_overflow:
            # If we have overflows in our residue numbers, wipe out all insertion codes. CHARMM
            # abuses the PDB format and usurps the insertion code (and later fields) to extend the
            # residue number. I don't feel bad discarding insertion code information when the PDB
            # format is abused like this. It's really never used for overflowing PDBs in my
            # experience anyway
            for residue in inst.struct.residues:
                residue.insertion_code = ''

        return inst.struct

    def _parse_open_file(self, fileobj):
        method_dispatch = {
            'REMARK': self._parse_remarks,
            'ATOM': self._parse_atom_record,
            'ANISOU': self._parse_anisou_record,
            'TER': self._parse_ter_record,
            'HETATM': self._parse_atom_record,
            'EXPDTA': self._parse_experimental,
            'AUTHOR': self._parse_author,
            'JRNL': self._parse_journal,
            'KEYWDS': self._parse_keywords,
            'CONECT': self._parse_connect_record,
            'LINK': self._parse_links,
            'CRYST1': self._parse_cryst1,
            'MODEL': self._new_model,
            'ENDMDL': self._end_model,
        }

        for line in self.fileobj:
            rec = line[:6].strip()
            if rec in method_dispatch:
                method_dispatch[rec](line)

    def _parse_remarks(self, line):
        """ Parse the various remarks, which may contain various metadata """
        if line[6:10] == ' 290':
            self._symmetry_lines.append(line)
        elif line[6:10] == ' 900':
            rematch = self._relatere.match(line[11:])
            if rematch:
                self.struct.related_entries.append(rematch.groups())
        elif line[6:10] == '   2': # Resolution
            if (not line[11:22].strip() or self.struct.resolution is not None
                or line[23:38] == 'NOT APPLICABLE.'):
                return
            elif line[11:22] !=  'RESOLUTION.':
                warnings.warn('Unrecognized RESOLUTION record in PDB file: %s' % line.strip())
            else:
                try:
                    self.struct.resolution = float(line[23:30])
                except ValueError:
                    warnings.warn('Trouble converting resolution (%s) to float' % line[23:30])

    def _parse_links(self, line):
        self._link_lines.append(line)

    def _parse_keywords(self, line):
        self.struct.keywords += '%s,' % line[10:].strip()

    def _parse_ter_record(self, *args):
        if self._current_model_number == 1:
            self._last_atom.residue.ter = True

    def _parse_experimental(self, line):
        if self.struct.experimental:
            self.struct.experimental += ' %s' % line[6:].strip()
        else:
            self.struct.experimental += line[6:].strip()

    def _parse_author(self, line):
        if self.struct.authors:
            self.struct.authors += ' %s' % line[10:].strip()
        else:
            self.struct.authors = line[10:].strip()

    def _parse_journal(self, line):
        part = line[12:16]
        if part == 'AUTH':
            self.struct.journal_authors += line[19:].strip()
        elif part == 'TITL':
            self.struct.title += ' %s' % line[19:].strip()
        elif part == 'REF ':
            self.struct.journal += ' %s' % line[19:47].strip()
            if not line[16:18].strip():
                self.struct.volume = line[51:55].strip()
                self.struct.page = line[56:61].strip()
                try:
                    self.struct.year = int(line[62:66])
                except ValueError:
                    # Shouldn't happen, but don't throw a fit
                    pass
        elif part == 'PMID':
            self.struct.pmid = line[19:].strip()
        elif part == 'DOI ':
            self.struct.doi = line[19:].strip()

    def _parse_cryst1(self, line):
        """ Parses unit cell information from the CRYST1 record """
        a = float(line[6:15])
        b = float(line[15:24])
        c = float(line[24:33])
        try:
            A = float(line[33:40])
            B = float(line[40:47])
            C = float(line[47:54])
        except (IndexError, ValueError):
            A = B = C = 90.0
        self.struct.box = [a, b, c, A, B, C]
        self.struct.space_group = line[55:66].strip()

    @staticmethod
    def _parse_atom_parts_1(line):
        return dict(
            number=line[6:11], name=line[12:16].strip(), alternate_location=line[16].strip(),
            residue_name=line[17:21].strip(), chain=line[21].strip(),
            residue_number=line[22:26].strip(), insertion_code=line[26].strip(),
        )

    @classmethod
    def _parse_atom_parts(cls, line):
        """ Pulls out atom attributes from the line and packs it into a dict """
        # segment_id is CHARMM-specific
        parts = cls._parse_atom_parts_1(line)
        parts.update(dict(
            x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
            occupancy=line[54:60], bfactor=line[60:66], element=line[76:78],
            charge=line[78:80], segment_id=line[72:76].strip(),
        ))

        elem = '%-2s' % parts['element']
        # Make sure the space is at the end
        elem = elem[1] + ' ' if elem[0] == ' ' else elem
        try:
            elem = (elem[0].upper() + elem[1].lower()).strip()
            atomic_number = AtomicNum[elem]
        except KeyError:
            elem = element_by_name(parts['name'])
            atomic_number = AtomicNum[elem]
        parts['atomic_number'] = atomic_number
        parts['mass'] = Mass[elem]
        parts['bfactor'] = try_convert(parts['bfactor'], float, 0.0)
        parts['occupancy'] = try_convert(parts['occupancy'], float, 0.0)
        parts['formal_charge'] = try_convert(parts['charge'], int, None)
        parts['charge'] = try_convert(parts['charge'], float, 0.0)

        return parts

    def _determine_residue_number(self, residue_number, line):
        if self._last_atom is not None:
            last_residue_number = self._last_atom.residue.number
        else:
            last_residue_number = 0
        extended_by = self._residue_number_field_extended_by
        if self._residue_indices_overflow:
            if extended_by:
                try:
                    self._last_residue_number_label = residue_number + line[26:26+extended_by]
                    residue_number = int(residue_number + line[26:26+extended_by])
                    return residue_number
                except ValueError:
                    pass
            else:
                if residue_number == self._last_residue_number_label:
                    return last_residue_number
                else:
                    self._last_residue_number_label = residue_number
                    return last_residue_number + 1
        if last_residue_number >= 9999 and residue_number == '1000' + '0' * extended_by and \
                line[26+extended_by] == '0':
            self._residue_indices_overflow = True
            self._residue_number_field_extended_by += 1
            self._last_residue_number_label = '1000' + '0' * self._residue_number_field_extended_by
            return int(self._last_residue_number_label)
        elif last_residue_number == 9999 and residue_number != '9999':
            # This is the first time we notice residue indices overflowing
            self._residue_indices_overflow = True
            self._last_residue_number_label = residue_number
            return last_residue_number + 1
        else:
            self._last_residue_number_label = residue_number
            return int(residue_number)

    def _parse_anisou_record(self, line):
        parts = self._parse_atom_parts_1(line)
        try:
            u11 = int(line[28:35])
            u22 = int(line[35:42])
            u33 = int(line[42:49])
            u12 = int(line[49:56])
            u13 = int(line[56:63])
            u23 = int(line[63:70])
        except ValueError:
            warnings.warn('Problem parsing anisotropic factors from ANISOU record from: ' + line,
                          PDBWarning)
        else:
            anisou = np.array([u11/1e4, u22/1e4, u33/1e4, u12/1e4, u13/1e4, u23/1e4])
            key = self._make_atom_key_from_parts(parts, all_parts=True)
            self._anisou_records[key] = (anisou, line)

    def _determine_atom_number(self, atom_number):
        if self._last_atom is not None:
            self._atom_indices_overflow = (self._atom_indices_overflow or
                                           self._last_atom.number >= 99999)
        if self._atom_indices_overflow and self._last_atom is not None:
            return self._last_atom.number + 1
        try:
            return int(atom_number)
        except ValueError:
            return self._last_atom.number + 1 if self._last_atom is not None else 1

    @classmethod
    def _make_atom_key_from_parts(cls, atom_parts, all_parts=False):
        """ If all_parts is False, only key from the values that determine alternate locations """
        return cls.AtomLookupKey(
            name=atom_parts['name'],
            number=atom_parts['number'] if all_parts else None,
            residue_name=atom_parts['residue_name'],
            residue_number=atom_parts['residue_number'],
            chain=atom_parts['chain'],
            insertion_code=atom_parts['insertion_code'],
            segment_id=atom_parts.get('segment_id', ''),
            alternate_location=atom_parts['alternate_location'] if all_parts else None,
        )

    def _parse_atom_record(self, line):
        """ Parses an atom record from a PDB file """
        atom_parts = self._parse_atom_parts(line)
        residue_number = self._determine_residue_number(atom_parts['residue_number'], line)
        atom_number = self._determine_atom_number(atom_parts['number'])

        AtomClass = ExtraPoint if atom_parts['name'] in ('EP', 'LP') else Atom
        atom = AtomClass(atomic_number=atom_parts['atomic_number'], name=atom_parts['name'],
                         charge=atom_parts['charge'], mass=atom_parts['mass'],
                         occupancy=atom_parts['occupancy'], bfactor=atom_parts['bfactor'],
                         altloc=atom_parts['alternate_location'], number=atom_number,
                         formal_charge=atom_parts['formal_charge'])
        atom.xx = atom_parts['x']
        atom.xy = atom_parts['y']
        atom.xz = atom_parts['z']
        attribute_key = self._make_atom_key_from_parts(atom_parts)
        all_attribute_key = self._make_atom_key_from_parts(atom_parts, all_parts=True)
        current_atom = self._atom_map_from_attributes.get(attribute_key, None)
        if (current_atom is not None and atom_parts['alternate_location'] in _ascii_letters_set and
            not self._atom_indices_overflow and not self._residue_indices_overflow):
            if self._current_model_number == 1:
                current_atom.other_locations[atom_parts['alternate_location']] = atom
                self._atom_map_to_parent[atom] = current_atom
                if atom_number not in self._atom_map_from_atom_number:
                    self._atom_map_from_atom_number[atom_number] = atom
                # altloc atoms should be reachable in the all_attributes map
                self._atom_map_from_all_attributes[all_attribute_key] = atom
                return
            elif not current_atom in self._model1_atoms_in_structure:
                # This is if the atom is not in the structure, but in the alternate locations. We
                # don't currently store coordinates for those atoms beyond the first frame. Note
                # that this should be incredibly rare, since alt-locs are common in static structure
                # determination, like X-Ray crystallography. Ensemble methods like NMR won't have
                # alternate locations in addition to multiple frames, so this is not expected to be
                # an impactful limitation
                return
        current_atom_index = self._model_atom_counts[-1]
        self._model_atom_counts[-1] += 1
        if self._current_model_number == 1:
            self._atom_map_from_all_attributes[all_attribute_key] = atom
            self._atom_map_from_attributes[attribute_key] = atom
            self._atom_map_from_atom_number[atom_number] = atom
            self.struct.add_atom(
                atom,
                atom_parts['residue_name'],
                residue_number,
                atom_parts['chain'],
                atom_parts['insertion_code'],
                atom_parts['segment_id'],
            )
            self._model1_atoms_in_structure.add(atom)
            self._last_atom = atom
        else:
            try:
                atom_from_first_model = self.struct.atoms[current_atom_index]
            except IndexError:
                raise PDBError('Atom number mismatch between models')
            if (atom_from_first_model.residue.name != atom_parts['residue_name'] or
                atom_from_first_model.name != atom_parts['name']):
                raise PDBError(
                    f'Atom/residue name mismatch in different models in model {self._current_model_number} [{line.strip()}]!'
                )
        if self._current_model_number == 1 or current_atom in self._model1_atoms_in_structure:
            self._coordinates[-1].extend([atom.xx, atom.xy, atom.xz])

    def _new_model(self, line):
        if self._current_model_number == 1 and len(self.struct.atoms) == 0:
            return # MODEL 1
        if self._model_open:
            warnings.warn(
                f'{line.strip()} begun before last model ended. Assuming it is ending',
                PDBWarning,
            )
            self._end_model(line)
        self._coordinates.append([])
        self._model_atom_counts.append(0)
        self._model_open = True
        self._current_model_number += 1

    def _end_model(self, line):
        """
        Ends the current model and validates the processed model is the same size as the models
        that came before this one
        """
        if not self._model_open:
            raise PDBError('No model was begun before ENDMDL was encountered')
        self._model_open = False
        if len(self._atom_map_from_attributes) == 0:
            raise PDBError('No atoms found in model')
        if len(self._coordinates[-1]) != 3 * len(self._atom_map_from_attributes):
            raise PDBError(f'Coordinate mismatch in model {self._current_model_number}')

    def _parse_connect_record(self, line):
        """
        Parses the CONECT records and creates the bond. According to the format spec, the first two
        atom indexes are required. The final 3 are optional.
        """
        origin_index = try_convert(line[6:11], int)
        index_1 = try_convert(line[11:16], int)
        index_2 = try_convert(line[16:21], int)
        index_3 = try_convert(line[21:26], int)
        index_4 = try_convert(line[26:31], int)
        if origin_index is None or index_1 is None:
            warnings.warn(
                f'Bad CONECT record -- not enough atom indexes in line: {line}', PDBWarning
            )
            return
        origin_atom = self._atom_map_from_atom_number.get(origin_index, None)
        atom_1 = self._atom_map_from_atom_number.get(index_1, None)
        atom_2 = self._atom_map_from_atom_number.get(index_2, None)
        atom_3 = self._atom_map_from_atom_number.get(index_3, None)
        atom_4 = self._atom_map_from_atom_number.get(index_4, None)
        if origin_atom is None or atom_1 is None:
            warnings.warn(
                f'CONECT record - could not find atoms {origin_index} and/or {index_1} to connect. Line: {line}',
                PDBWarning,
            )
            return
        origin_atom = self._atom_map_to_parent.get(origin_atom, origin_atom)
        for partner in (atom_1, atom_2, atom_3, atom_4):
            partner = self._atom_map_to_parent.get(partner, partner)
            if partner is None or partner in origin_atom.bond_partners:
                continue
            self.struct.bonds.append(Bond(origin_atom, partner))

    def _process_structure_symmetry(self):
        if self._symmetry_lines:
            data = []
            for line in self._symmetry_lines:
                if line.strip().startswith('REMARK 290   SMTRY'):
                    data.append(line.split()[4:])
            tensor = np.asarray(data, dtype='f8')
            self.struct.symmetry = Symmetry(tensor)

    def _process_link_records(self):
        for line in self._link_lines:
            atom_1_parts = self._parse_atom_parts_1(line)
            atom_2_parts = self._parse_atom_parts_1(line[30:])

            symop1 = line[59:65].strip()
            symop2 = line[66:72].strip()
            try:
                length = float(line[73:78])
            except ValueError:
                warnings.warn(f'Malformed LINK line (bad distance): {line}', PDBWarning)
                continue

            key1 = self._make_atom_key_from_parts(atom_1_parts)
            key2 = self._make_atom_key_from_parts(atom_2_parts)
            try:
                a1 = self._atom_map_from_attributes[key1]
                a2 = self._atom_map_from_attributes[key2]
            except KeyError:
                warnings.warn(f'Could not find link atoms {key1} and {key2}', PDBWarning)
            else:
                self.struct.links.append(Link(a1, a2, length, symop1, symop2))

    def _assign_anisou_to_atoms(self):
        """ Assigns the ANISOU tensors to the atoms they belong to """
        for key, (anisou_tensor, line) in self._anisou_records.items():
            try:
                self._atom_map_from_all_attributes[key].anisou = anisou_tensor
            except KeyError:
                warnings.warn(
                    f'Could not find atom belonging to anisou tensor with key {key}. Line: {line}',
                    PDBWarning
                )

    def _postprocess_metadata(self):
        self.struct.keywords = [s.strip() for s in self.struct.keywords.split(',') if s.strip()]
        self.struct.journal = self.struct.journal.strip()
        self.struct.title = self.struct.title.strip()

    @staticmethod
    def write(struct, dest, renumber=True, coordinates=None, altlocs='all',
              write_anisou=False, charmm=False, use_hetatoms=True,
              standard_resnames=False, increase_tercount=True, write_links=False):
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
        renumber : bool, optional, default True
            If True, renumber the atoms and residues sequentially as they are
            stored in the structure.  If False, use the original numbering if
            it was assigned previously.
        coordinates : array-like of float, optional
            If provided, these coordinates will be written to the PDB file
            instead of the coordinates stored in the structure. These
            coordinates should line up with the atom order in the structure
            (not necessarily the order of the "original" PDB file if they
            differ)
        altlocs : str, optional, default 'all'
            Keyword controlling which alternate locations are printed to the
            resulting PDB file. Allowable options are:

                - 'all' : print all alternate locations
                - 'first' : print only the first alternate locations
                - 'occupancy' : print the one with the largest occupancy. If two
                  conformers have the same occupancy, the first one to occur is
                  printed

            Input is case-insensitive, and partial strings are permitted as long
            as it is a substring of one of the above options that uniquely
            identifies the choice.
        write_anisou : bool, optional, default False
            If True, an ANISOU record is written for every atom that has one. If
            False, ANISOU records are not written.
        charmm : bool, optional, default False
            If True, SEGID will be written in columns 73 to 76 of the PDB file
            in the typical CHARMM-style PDB output. This will be omitted for any
            atom that does not contain a SEGID identifier.
        use_hetatoms: bool, optional, default True
            If True, certain atoms will have the HETATM tag instead of ATOM
            as per the PDB-standard. 
        standard_resnames : bool, optional, default False
            If True, common aliases for various amino and nucleic acid residues
            will be converted into the PDB-standard values.
        increase_tercount : bool, optional, default True
            If True, the TER atom number field increased by one compared to
            atom card preceding it; this conforms to PDB standard.
        write_links : bool, optional, default False
            If True, any LINK records stored in the Structure will be written to
            the LINK records near the top of the PDB file. If this is True, then
            renumber *must* be False or a ValueError will be thrown

        Notes
        -----
        If multiple coordinate frames are present, these will be written as
        separate models (but only the unit cell from the first model will be
        written, as the PDB standard dictates that only one set of unit cells
        shall be present).
        """
        # Determine if we have *any* atom or residue numbers set. If none of
        # them are set, force renumbering
        no_atom_numbers_assigned = {a.number for a in struct.atoms} == {-1}
        no_residue_numbers_assigned = {r.number for r in struct.residues} == {-1}
        renumber = renumber or (no_atom_numbers_assigned and no_residue_numbers_assigned)
        if renumber and write_links:
            raise ValueError('write_links requires renumber=False AND original (not implied) '
                             'numbers to be assigned')
        if altlocs.lower() == 'all'[:len(altlocs)]:
            altlocs = 'all'
        elif altlocs.lower() == 'first'[:len(altlocs)]:
            altlocs = 'first'
        elif altlocs.lower() == 'occupancy'[:len(altlocs)]:
            altlocs = 'occupancy'
        else:
            raise ValueError("Illegal value of occupancy [%s]; expected 'all', "
                             "'first', or 'occupancy'" % altlocs)
        own_handle = False
        if not hasattr(dest, 'write'):
            dest = genopen(dest, 'w')
            own_handle = True
        if charmm:
            atomrec = ('ATOM  %5d %-4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f'
                       '%6.2f      %-4s%2s%-2s\n')
            anisourec = 'ANISOU%5d %-4s%1s%-4s%1s%4d%1s %7d%7d%7d%7d%7d%7d      %2s%-2s\n'
            terrec = 'TER   %5d      %-4s%1s%4d\n'
            reslen = 4
        else:
            atomrec = ('ATOM  %5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f'
                       '%6.2f      %-4s%2s%-2s\n')
            anisourec = 'ANISOU%5d %-4s%1s%-3s %1s%4d%1s %7d%7d%7d%7d%7d%7d      %2s%-2s\n'
            terrec = ('TER   %5d      %-3s %1s%4d\n')
            reslen = 3
        linkrec = ('LINK        %-4s%1s%-3s %1s%4d%1s               '
                   '%-4s%1s%-3s %1s%4d%1s  %6s %6s %5.2f\n')
        hetatomrec = atomrec.replace('ATOM  ', 'HETATM') if use_hetatoms else atomrec
        if struct.box is not None:
            dest.write('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4s\n' % (
                    struct.box[0], struct.box[1], struct.box[2], struct.box[3],
                    struct.box[4], struct.box[5], struct.space_group, ''))
        if struct.symmetry is not None:
            fmt = '%d%4d%10.6f%10.6f%10.6f%15.5f\n'
            for index, arr in enumerate(struct.symmetry.data):
                arr_list = [1 + index % 3, 1 + index//3] + arr.tolist()
                symm_line = "REMARK 290   SMTRY" + fmt % tuple(arr_list)
                dest.write(symm_line)
        if coordinates is not None:
            coords = np.asanyarray(coordinates)
            try:
                coords = coords.reshape((-1, len(struct.atoms), 3))
            except ValueError:
                raise TypeError("Coordinates has unexpected shape")
        else:
            coords = struct.get_coordinates('all')
            if coords is None:
                raise ValueError('Cannot write PDB file with no coordinates')
        # Create a function to process each atom and return which one we want
        # to print, based on our alternate location choice
        if altlocs == 'all':
            def print_atoms(atom, coords):
                return atom, atom.other_locations, coords[atom.idx]
        elif altlocs == 'first':
            def print_atoms(atom, coords):
                return atom, dict(), coords[atom.idx]
        elif altlocs == 'occupancy':
            def print_atoms(atom, coords):
                occ = atom.occupancy
                a = atom
                for key, item in atom.other_locations.items():
                    if item.occupancy > occ:
                        occ = item.occupancy
                        a = item
                return a, dict(), [a.xx, a.xy, a.xz]
        else:
            assert False, 'Should not be here'
        if standard_resnames:
            standardize = lambda x: _standardize_resname(x)[:reslen]
        else:
            standardize = lambda x: (x[:reslen], _is_hetatm(x))
        nmore = 0 # how many *extra* atoms have been added?
        last_number = 0
        if write_links:
            for link in struct.links:
                rec = (
                    _format_atom_name_for_pdb(link.atom1),
                    link.atom1.altloc,
                    link.atom1.residue.name,
                    link.atom1.residue.chain,
                    link.atom1.residue.number,
                    link.atom1.residue.insertion_code,

                    _format_atom_name_for_pdb(link.atom2),
                    link.atom2.altloc,
                    link.atom2.residue.name,
                    link.atom2.residue.chain,
                    link.atom2.residue.number,
                    link.atom2.residue.insertion_code,

                    link.symmetry_op1,
                    link.symmetry_op2,
                    link.length,
                )
                dest.write(linkrec % rec)
        for model, coord in enumerate(coords):
            if coords.shape[0] > 1:
                dest.write('MODEL      %5d\n' % (model+1))
            for res in struct.residues:
                if renumber:
                    atoms = res.atoms
                else:
                    atoms = sorted(res.atoms, key=lambda atom: atom.number)
                if charmm:
                    segid = (res.segid or res.chain)[:4]
                else:
                    segid = ''
                for atom in atoms:
                    pa, others, (x, y, z) = print_atoms(atom, coord)
                    # Figure out the serial numbers we want to print
                    if renumber:
                        anum = _number_truncated_to_n_digits(atom.idx + 1 + nmore, 5)
                        rnum = _number_truncated_to_n_digits(res.idx + 1, 4)
                    else:
                        anum = _number_truncated_to_n_digits(pa.number, 5)
                        rnum = _number_truncated_to_n_digits(res.number, 4)
                    last_number = anum
                    # Do any necessary name munging to respect the PDB spec
                    aname = _format_atom_name_for_pdb(pa)
                    resname, hetatom = standardize(res.name)
                    if hetatom:
                        rec = hetatomrec
                    else:
                        rec = atomrec
                    dest.write(rec % (anum, aname, pa.altloc, resname,
                               res.chain[:1], rnum, res.insertion_code[:1],
                               x, y, z, pa.occupancy, pa.bfactor, segid,
                               Element[pa.atomic_number].upper(), ''))
                    if write_anisou and pa.anisou is not None:
                        anisou = [int(ani*1e4) for ani in pa.anisou]
                        dest.write(anisourec % (anum, aname, pa.altloc,
                                   resname, res.chain[:1], rnum,
                                   res.insertion_code[:1], anisou[0], anisou[1],
                                   anisou[2], anisou[3], anisou[4], anisou[5],
                                   Element[pa.atomic_number].upper(), ''))
                    for key in sorted(others.keys()):
                        oatom = others[key]
                        x, y, z = oatom.xx, oatom.xy, oatom.xz
                        if renumber:
                            nmore += 1
                            anum = (pa.idx + 1 + nmore)
                        else:
                            anum = oatom.number or last_number + 1
                        anum = anum - anum // 100000 * 100000
                        last_number = anum
                        aname = _format_atom_name_for_pdb(oatom)
                        dest.write(rec % (anum, aname, key, resname,
                                   res.chain[:1], rnum, res.insertion_code[:1],
                                   x, y, z, oatom.occupancy, oatom.bfactor, segid,
                                   Element[oatom.atomic_number].upper(), ''))
                        if write_anisou and oatom.anisou is not None:
                            anisou = [int(ani*1e4) for ani in oatom.anisou]
                            el = Element[oatom.atomic_number].upper()
                            dest.write(anisourec % (anum, aname,
                                oatom.altloc[:1], resname, res.chain[:1],
                                rnum, res.insertion_code[:1], anisou[0],
                                anisou[1], anisou[2], anisou[3],
                                anisou[4], anisou[5], el, ''))
                if res.ter or (len(struct.bonds) > 0 and _needs_ter_card(res)):
                    if increase_tercount:
                        dest.write(terrec % (anum+1, resname, res.chain, rnum))
                        if renumber:
                            nmore += 1
                        else:
                            last_number += 1
                    else:
                        dest.write(terrec % (anum, resname, res.chain, rnum))
            if coords.shape[0] > 1:
                dest.write('ENDMDL\n')

        dest.write("%-80s\n" % "END")
        if own_handle:
            dest.close()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CIFFile(metaclass=FileFormatType):
    """ Standard PDBx/mmCIF file format parser and writer """
    #===================================================

    @staticmethod
    def id_format(filename):
        """ Identifies the file type as a PDBx/mmCIF file

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is a PDBx/mmCIF file
        """
        f = genopen(filename)
        try:
            for line in f:
                if line.startswith('#'): continue
                if line[:5] == 'data_' and len(line.split()) == 1:
                    return True
                else:
                    return False
            return False
        finally:
            f.close()

    #===================================================

    @staticmethod
    def download(pdb_id, timeout=10, saveto=None):
        """
        Goes to the wwPDB website and downloads the requested PDBx/mmCIF,
        loading it as a :class:`Structure` instance

        Parameters
        ----------
        pdb_id : str
            The 4-letter PDB ID to try and download from the RCSB PDB database
        timeout : float, optional
            The number of seconds to wait before raising a timeout error.
            Default is 10 seconds
        saveto : str, optional
            If provided, this will be treated as a file name to which the PDB
            file will be saved. If None (default), no CIF file will be written.
            This will be a verbatim copy of the downloaded CIF file, unlike the
            somewhat-stripped version you would get by using
            :meth:`Structure.write_cif <parmed.structure.Structure.write_cif>`

        Returns
        -------
        struct : :class:`Structure <parmed.structure.Structure>`
            Structure instance populated by the requested PDBx/mmCIF

        Raises
        ------
        socket.timeout if the connection times out while trying to contact the
        FTP server

        IOError if there is a problem retrieving the requested PDB

        ImportError if the gzip module is not available

        TypeError if pdb_id is not a 4-character string
        """
        if not isinstance(pdb_id, str) or len(pdb_id) != 4:
            raise ValueError('pdb_id must be the 4-letter PDB code')

        pdb_id = pdb_id.lower()
                
        url = (
            'https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/%s/%s.cif.gz'
        ) % (pdb_id[1:3], pdb_id)
        
        fileobj = io.BytesIO()
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                fileobj.write(response.read())
        except (urllib.error.URLError, urllib.error.HTTPError) as err:
            raise IOError(f"Could not retrieve PDB ID {pdb_id}; {err}") from err
        
        fileobj.seek(0)
        fileobj = io.TextIOWrapper(gzip.GzipFile(fileobj=fileobj, mode='r'))
        if saveto is not None:
            with closing(genopen(saveto, 'w')) as f:
                f.write(fileobj.read())
            fileobj.seek(0)
        return CIFFile.parse(fileobj)

    #===================================================

    @classmethod
    def parse(cls, filename, skip_bonds=False, expanded_residue_template_match=False, all_residue_template_match=False):
        """
        Read a PDBx or mmCIF file and return a populated `Structure` class

        Parameters
        ----------
        filename : ``str or file-like``
            Name of PDB file to read, or a file-like object that can iterate
            over the lines of a PDB. Compressed file names can be specified and
            are determined by file-name extension (e.g., file.pdb.gz,
            file.pdb.bz2)
        skip_bonds : bool, optional
            If True, skip trying to assign bonds. This can save substantial time
            when parsing large files with non-standard residue names. However,
            no bonds are assigned. This is OK if, for instance, the CIF file is
            being parsed simply for its coordinates. Default is False.
        expanded_residue_template_match : bool, optional
            If True, an expanded set of residue templates taken from the chemical component database is used
            to assign bonds and atomic properties.
        all_residue_template_match : bool, optional
            If True, all residue templates taken from the chemical component database is used
            to assign bonds and atomic properties.

        Metadata
        --------
        The PDB parser also adds metadata to the returned Structure object that
        may be present in the PDB file

        experimental : ``str``
            EXPDTA record
        journal : ``str``
            JRNL record
        authors : ``str``
            AUTHOR records
        keywords : ``str``
            KEYWDS records
        doi : ``str``
            DOI from the JRNL record
        pmid : ``str``
            PMID from the JRNL record
        journal_authors : ``str``
            Author info from the JRNL record
        volume : ``str``
            Volume of the published article from the JRNL record
        page : ``str``
            Page of the published article from the JRNL record
        title : ``str``
            TITL section of the JRNL record
        year : ``str``
            Year that the article was published, from the JRNL record
        resolution : ``float``
            The X-RAY resolution in Angstroms, or None if not found
        related_entries : ``list of (str, str)``
            List of entries in other databases

        Returns
        -------
        structure1 [, structure2 [, structure3 [, ...] ] ]

        structure# : :class:`Structure`
            The Structure object initialized with all of the information from
            the PDBx/mmCIF file.  No bonds or other topological features are
            added by default. If multiple structures are defined in the CIF
            file, multiple Structure instances will be returned as a tuple.

        Raises
        ------
        ValueError if the file severely violates the PDB format specification.
        If this occurs, check the formatting on each line and make sure it
        matches the others.
        """
        if all_residue_template_match and not expanded_residue_template_match:
            raise ValueError("all_residue_template_match cannot be True if expanded_residue_template_match is False")
        if expanded_residue_template_match and skip_bonds:
            raise ValueError("extended_residue_template_match cannot be True if skip_bonds is also True")
        if isinstance(filename, str):
            own_handle = True
            fileobj = genopen(filename, 'r')
        else:
            own_handle = False
            fileobj = filename

        try:
            cifobj = PdbxReader(fileobj)
            data = []
            cifobj.read(data)
        finally:
            if own_handle: fileobj.close()

        structures = []
        for cont in data:
            struct = Structure()
            structures.append(struct)
            # Add metadata fields
            struct.experimental = struct.journal = struct.authors = ''
            struct.keywords = struct.doi = struct.pmid = ''
            struct.journal_authors = struct.volume = struct.title = ''
            struct.year = struct.resolution = None
            struct.related_entries = []

            # Now we have the data. First get the metadata if it exists
            exptl = cont.getObj('exptl')
            if exptl is not None:
                struct.experimental = exptl.getValue('method')
            auth = cont.getObj('audit_author')
            if auth is not None:
                nameidx = auth.getAttributeIndex('name')
                if nameidx != -1:
                    struct.authors = ', '.join([t[nameidx] for t in
                                               auth.getRowList()])
            reflns = cont.getObj('reflns')
            if reflns is not None:
                res = reflns.getValue('d_resolution_high')
                if res != '?':
                    try:
                        struct.resolution = float(res)
                    except ValueError:
                        warnings.warn('Could not convert resolution (%s) to float' % res)
            cite = cont.getObj('citation_author')
            if cite is not None:
                nameidx = cite.getAttributeIndex('name')
                if nameidx != -1:
                    journal_authors = []
                    for i in range(cite.getRowCount()):
                        a = cite.getRow(i)[nameidx]
                        if a not in journal_authors:
                            journal_authors.append(a)
                    struct.journal_authors = ', '.join(journal_authors)
            cite = cont.getObj('citation')
            if cite is not None:
                doiid = cite.getAttributeIndex('pdbx_database_id_DOI')
                pmiid = cite.getAttributeIndex('pdbx_database_id_PubMed')
                titlid = cite.getAttributeIndex('title')
                yearid = cite.getAttributeIndex('year')
                pageid = cite.getAttributeIndex('page_first')
                jrnlid = cite.getAttributeIndex('journal_abbrev')
                volid = cite.getAttributeIndex('journal_volume')
                rows = cite.getRowList()
                if doiid != -1:
                    struct.doi = ', '.join([row[doiid] for row in rows if row[doiid] != '?'])
                if pmiid != -1:
                    struct.pmid = ', '.join([row[pmiid] for row in rows if row[pmiid] != '?'])
                if titlid != -1:
                    struct.title = '; '.join([row[titlid] for row in rows])
                if yearid != -1:
                    struct.year = ', '.join([row[yearid] for row in rows])
                if pageid != -1:
                    struct.page = ', '.join([row[pageid] for row in rows])
                if jrnlid != -1:
                    struct.journal = '; '.join([row[jrnlid] for row in rows])
                if volid != -1:
                    struct.volume = ', '.join([row[volid] for row in rows])
            keywds = cont.getObj('struct_keywords')
            if keywds is not None:
                textid = keywds.getAttributeIndex('text')
                if textid != -1:
                    rows = keywds.getRowList()
                    struct.keywords = ', '.join([row[textid] for row in rows])
                    struct.keywords = [key.strip() for key in struct.keywords.split(',') if key.strip()]
            dbase = cont.getObj('pdbx_database_related')
            if dbase is not None:
                dbid = dbase.getAttributeIndex('db_id')
                nameid = dbase.getAttributeIndex('db_name')
                if dbid != -1 and nameid != -1:
                    rows = dbase.getRowList()
                    struct.related_entries = [(r[dbid],r[nameid]) for r in rows]
            # Now go through all the atoms. Any items that do *not* exist are
            # given an index of -1. So we append an empty string on the end of
            # each row so that the default value for any un-specified value is
            # the empty string. This avoids needing any conditionals inside the
            # loop
            atoms = cont.getObj('atom_site')
            atnumid = atoms.getAttributeIndex('id')
            elemid = atoms.getAttributeIndex('type_symbol')
            atnameid = atoms.getAttributeIndex('auth_atom_id')
            altlocid = atoms.getAttributeIndex('label_alt_id')
            resnameid = atoms.getAttributeIndex('auth_comp_id')
            chainid = atoms.getAttributeIndex('auth_asym_id')
            resnumid = atoms.getAttributeIndex('auth_seq_id')
            inscodeid = atoms.getAttributeIndex('pdbx_PDB_ins_code')
            formal_charge_id = atoms.getAttributeIndex("pdbx_formal_charge")
            xid = atoms.getAttributeIndex('Cartn_x')
            yid = atoms.getAttributeIndex('Cartn_y')
            zid = atoms.getAttributeIndex('Cartn_z')
            occupid = atoms.getAttributeIndex('occupancy')
            bfactorid = atoms.getAttributeIndex('B_iso_or_equiv')
            modelid = atoms.getAttributeIndex('pdbx_PDB_model_num')
            origmodel = None
            lastmodel = None
            all_coords = []
            xyz = []
            atommap = dict()
            last_atom = Atom()
            for i in range(atoms.getRowCount()):
                row = atoms.getRow(i) + ['']
                atnum = int(row[atnumid])
                elem = row[elemid]
                atname = row[atnameid]
                altloc = row[altlocid]
                if altloc == '.': altloc = ''
                resname = row[resnameid]
                chain = row[chainid]
                resnum = int(row[resnumid])
                inscode = row[inscodeid]
                if inscode in '?.': inscode = ''
                model = int(row[modelid])
                if origmodel is None:
                    origmodel = lastmodel = model
                x, y, z = float(row[xid]), float(row[yid]), float(row[zid])
                occup = float(row[occupid])
                bfactor = float(row[bfactorid])
                try:
                    formal_charge = int(row[formal_charge_id]) if formal_charge_id != -1 else None
                except ValueError:
                    formal_charge = None
                try:
                    atsym = elem.strip().title()
                    atomic_number = AtomicNum[atsym]
                    mass = Mass[atsym]
                except KeyError:
                    # Now try based on the atom name... but don't try too hard
                    # (e.g., don't try to differentiate b/w Ca and C)
                    try:
                        atomic_number = AtomicNum[atname.strip()[0].upper()]
                        mass = Mass[atname.strip()[0].upper()]
                    except KeyError:
                        try:
                            sym = atname.strip()[:2]
                            sym = '%s%s' % (sym[0].upper(), sym[1].lower())
                            atomic_number = AtomicNum[sym]
                            mass = Mass[sym]
                        except KeyError:
                            atomic_number = 0 # give up
                            mass = 0.0
                if atname.startswith('EP') or atname.startswith('LP'):
                    atom = ExtraPoint(atomic_number=atomic_number, name=atname,
                                mass=mass, occupancy=occup, bfactor=bfactor,
                                altloc=altloc, number=atnum)
                else:
                    atom = Atom(atomic_number=atomic_number, name=atname,
                                mass=mass, occupancy=occup, bfactor=bfactor,
                                altloc=altloc, number=atnum, formal_charge=formal_charge)
                atom.xx, atom.xy, atom.xz = x, y, z
                if (_compare_atoms(last_atom, atom, resname, resnum, chain, '', inscode) and altloc):
                    atom.residue = last_atom.residue
                    last_atom.other_locations[altloc] = atom
                else:
                    if model == origmodel:
                        # Only add the atoms once
                        struct.add_atom(atom, resname, resnum, chain, inscode)
                    last_atom = atom
                    if model == lastmodel:
                        xyz.extend([x, y, z])
                    else:
                        if all_coords and len(xyz) != len(all_coords[-1]):
                            raise ValueError('All frames must have same number of atoms')
                        all_coords.append(xyz)
                        xyz = [x, y, z]
                        lastmodel = model
                # Keep a mapping in case we need to go back and add attributes,
                # like anisotropic b-factors
                if model == origmodel:
                    key = (resnum,resname,inscode,chain,atnum,altloc,atname)
                    atommap[key] = atom
            # Check for unit cell parameters
            cell = cont.getObj('cell')
            if cell is not None:
                aid = cell.getAttributeIndex('length_a')
                bid = cell.getAttributeIndex('length_b')
                cid = cell.getAttributeIndex('length_c')
                alphaid = cell.getAttributeIndex('angle_alpha')
                betaid = cell.getAttributeIndex('angle_beta')
                gammaid = cell.getAttributeIndex('angle_gamma')
                row = cell.getRow(0)
                struct.box = np.array(
                    [float(row[aid]), float(row[bid]), float(row[cid]),
                     float(row[alphaid]), float(row[betaid]), float(row[gammaid])]
                )
            symmetry = cont.getObj('symmetry')
            if symmetry is not None:
                spaceid = symmetry.getAttributeIndex('space_group_name_H-M')
                row = symmetry.getRow(0)
                if spaceid != -1:
                    struct.space_group = row[spaceid]
            # Check for anisotropic B-factors
            anisou = cont.getObj('atom_site_anisotrop')
            if anisou is not None:
                atnumid = anisou.getAttributeIndex('id')
                atnameid = anisou.getAttributeIndex('pdbx_auth_atom_id')
                altlocid = anisou.getAttributeIndex('pdbx_label_alt_id')
                resnameid = anisou.getAttributeIndex('pdbx_auth_comp_id')
                chainid = anisou.getAttributeIndex('pdbx_auth_asym_id')
                resnumid = anisou.getAttributeIndex('pdbx_auth_seq_id')
                inscodeid = anisou.getAttributeIndex('pdbx_PDB_ins_code')
                u11id = anisou.getAttributeIndex('U[1][1]')
                u22id = anisou.getAttributeIndex('U[2][2]')
                u33id = anisou.getAttributeIndex('U[3][3]')
                u12id = anisou.getAttributeIndex('U[1][2]')
                u13id = anisou.getAttributeIndex('U[1][3]')
                u23id = anisou.getAttributeIndex('U[2][3]')
                if -1 in (atnumid, atnameid, altlocid, resnameid, chainid,
                          resnumid, u11id, u22id, u33id, u12id, u13id, u23id):
                    warnings.warn('Incomplete anisotropic B-factor CIF section. Skipping', PDBWarning)
                else:
                    try:
                        for i in range(anisou.getRowCount()):
                            row = anisou.getRow(i) + ['']
                            atnum = int(row[atnumid])
                            atname = row[atnameid]
                            altloc = row[altlocid]
                            resname = row[resnameid]
                            chain = row[chainid]
                            resnum = int(row[resnumid])
                            inscode = row[inscodeid]
                            u11 = float(row[u11id])
                            u22 = float(row[u22id])
                            u33 = float(row[u33id])
                            u12 = float(row[u12id])
                            u13 = float(row[u13id])
                            u23 = float(row[u23id])
                            altloc = '' if altloc == '.' else altloc
                            inscode = '' if inscode in '?.' else inscode
                            key = (resnum, resname, inscode, chain, atnum, altloc, atname)
                            atommap[key].anisou = np.array([u11, u22, u33, u12, u13, u23])
                    except (ValueError, KeyError):
                        # If at least one went wrong, set them all to None
                        for key, atom in atommap.items():
                            atom.anisou = None
                        warnings.warn('Problem processing anisotropic B-factors. Skipping', PDBWarning)
            conn = cont.getObj("struct_conn")
            if conn is not None:
                try:
                    cls._process_struct_conn(struct, conn)
                except IndexError:
                    pass
            if xyz:
                if len(xyz) != len(struct.atoms) * 3:
                    raise ValueError('Corrupt CIF; all models must have the same atoms')
                all_coords.append(xyz)
            if all_coords:
                struct._coordinates = np.array(all_coords).reshape((-1, len(struct.atoms), 3))

        # Make sure we assign bonds for all the structures we parsed
        if not skip_bonds:
            for struct in structures:
                if expanded_residue_template_match:
                    residue_libraries = [get_standard_residue_template_library()]
                    if all_residue_template_match:
                        residue_libraries.append(get_nonstandard_ccd_residues())
                    struct.assign_bonds(*residue_libraries)
                    _map_atomic_properties_from_templates(struct, *residue_libraries)
                    _map_bond_types_from_templates(struct, *residue_libraries)
                else:
                    struct.assign_bonds()
        # Build the return value
        if len(structures) == 1:
            return structures[0]
        return tuple(structures)

    #===================================================

    @staticmethod
    def write(struct, dest, renumber=True, coordinates=None,
              altlocs='all', write_anisou=False, standard_resnames=False):
        """
        Write a PDB file from the current Structure instance

        Parameters
        ----------
        struct : :class:`Structure`
            The structure from which to write the PDBx/mmCIF file
        dest : ``str or file-like``
            Either a file name or a file-like object containing a `write`
            method to which to write the PDB file. If it is a filename that
            ends with .gz or .bz2, a compressed version will be written using
            either gzip or bzip2, respectively.
        renumber : ``bool``
            If True, renumber the atoms and residues sequentially as they are
            stored in the structure.  If False, use the original numbering if
            it was assigned previously
        coordinates : ``array-like of float``
            If provided, these coordinates will be written to the PDB file
            instead of the coordinates stored in the structure. These
            coordinates should line up with the atom order in the structure
            (not necessarily the order of the "original" PDB file if they
            differ)
        altlocs : ``str``
            Keyword controlling which alternate locations are printed to the
            resulting PDB file. Allowable options are:

                - 'all' : (default) print all alternate locations
                - 'first' : print only the first alternate locations
                - 'occupancy' : print the one with the largest occupancy. If two
                  conformers have the same occupancy, the first one to occur is
                  printed

            Input is case-insensitive, and partial strings are permitted as long
            as it is a substring of one of the above options that uniquely
            identifies the choice.
        write_anisou : ``bool``
            If True, an ANISOU record is written for every atom that has one. If
            False, ANISOU records are not written
        standard_resnames : bool, optional
            If True, common aliases for various amino and nucleic acid residues
            will be converted into the PDB-standard values. Default is False

        Notes
        -----
        If multiple coordinate frames are present, these will be written as
        separate models (but only the unit cell from the first model will be
        written, as the PDBx standard dictates that only one set of unit cells
        shall be present).
        """
        if altlocs.lower() == 'all'[:len(altlocs)]:
            altlocs = 'all'
        elif altlocs.lower() == 'first'[:len(altlocs)]:
            altlocs = 'first'
        elif altlocs.lower() == 'occupancy'[:len(altlocs)]:
            altlocs = 'occupancy'
        else:
            raise ValueError(f"Illegal value of occupancy [{altlocs}]; expected 'all', 'first', or 'occupancy'")
        own_handle = False
        if not hasattr(dest, 'write'):
            dest = genopen(dest, 'w')
            own_handle = True
        # Make the main container
        cont = containers.DataContainer('cell')
        # Add cell info if applicable
        if struct.box is not None:
            cell = containers.DataCategory('cell')
            cell.appendAttribute('length_a')
            cell.appendAttribute('length_b')
            cell.appendAttribute('length_c')
            cell.appendAttribute('angle_alpha')
            cell.appendAttribute('angle_beta')
            cell.appendAttribute('angle_gamma')
            cell.append(struct.box[:])
            cont.append(cell)
        # symmetry
        sym = containers.DataCategory('symmetry')
        sym.appendAttribute('space_group_name_H-M')
        sym.append([struct.space_group])
        cont.append(sym)
        if coordinates is not None:
            coords = np.asanyarray(coordinates)
            try:
                coords = coords.reshape((-1, len(struct.atoms), 3))
            except ValueError:
                raise TypeError("Coordinates has unexpected shape")
        else:
            coords = struct.get_coordinates('all')
            if coords is None:
                raise ValueError('Cannot write CIF file with no coordinates')
        # Create a function to process each atom and return which one we want
        # to print, based on our alternate location choice
        if altlocs == 'all':
            def print_atoms(atom, coords):
                return atom, atom.other_locations, coords[atom.idx]
        elif altlocs == 'first':
            def print_atoms(atom, coords):
                return atom, dict(), coords[atom.idx]
        elif altlocs == 'occupancy':
            def print_atoms(atom, coords):
                occ = atom.occupancy
                a = atom
                for key, item in atom.other_locations.items():
                    if item.occupancy > occ:
                        occ = item.occupancy
                        a = item
                return a, dict(), [a.xx, a.xy, a.xz]
        else:
            assert False, 'Should not be here'
        if standard_resnames:
            standardize = lambda x: _standardize_resname(x)
        else:
            standardize = lambda x: (x, _is_hetatm(x))
        # Now add the atom section. Include all names that the CIF standard
        # usually includes, but put '?' in sections that contain data we don't
        # store in the Structure, Residue, or Atom classes
        cifatoms = containers.DataCategory('atom_site')
        cont.append(cifatoms)
        cifatoms.setAttributeNameList(
            ['group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id', 'label_asym_id',
             'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy',
             'B_iso_or_equiv', 'Cartn_x_esd', 'Cartn_y_esd', 'Cartn_z_esd', 'occupancy_esd', 'B_iso_or_equiv_esd',
             'pdbx_formal_charge', 'auth_seq_id', 'auth_comp_id', 'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num']
        )
        write_anisou = write_anisou and any(atom.anisou is not None
                                            for atom in struct.atoms)
        if write_anisou:
            cifanisou = containers.DataCategory('atom_site_anisotrop')
            cont.append(cifanisou)
            cifanisou.setAttributeNameList(
                ['id', 'type_symbol', 'pdbx_label_atom_id', 'pdbx_label_alt_id', 'pdbx_label_comp_id',
                 'pdbx_label_asym_id', 'pdbx_label_seq_id', 'U[1][1]', 'U[2][2]', 'U[3][3]', 'U[1][2]', 'U[1][3]',
                 'U[2][3]', 'U[1][1]_esd', 'U[2][2]_esd', 'U[3][3]_esd', 'U[1][2]_esd', 'U[1][3]_esd', 'U[2][3]_esd',
                 'pdbx_auth_seq_id', 'pdbx_auth_comp_id', 'pdbx_auth_asym_id', 'pdbx_auth_atom_id']
            )
        nmore = 0 # how many *extra* atoms have been added?
        last_number = 0
        last_rnumber = 0
        for model, coord in enumerate(coords):
            for res in struct.residues:
                if renumber:
                    atoms = res.atoms
                else:
                    atoms = sorted(res.atoms, key=lambda atom: atom.number)
                resname, hetatom = standardize(res.name)
                if hetatom:
                    atomrec = 'HETATM'
                else:
                    atomrec = 'ATOM  '
                for atom in atoms:
                    pa, others, (x, y, z) = print_atoms(atom, coord)
                    # Figure out the serial numbers we want to print
                    if renumber:
                        anum = (atom.idx + 1 + nmore)
                        rnum = (res.idx + 1)
                    else:
                        anum = (pa.number or last_number + 1)
                        rnum = (atom.residue.number or last_rnumber + 1)
                    last_number = anum
                    last_rnumber = rnum
                    cifatoms.append(
                        [atomrec, anum, Element[pa.atomic_number].upper(), pa.name, pa.altloc, resname, res.chain, '?',
                         rnum, res.insertion_code, x, y, z, pa.occupancy, pa.bfactor, '?', '?', '?', '?', '?', '', rnum,
                         resname, res.chain, pa.name, str(model+1)]
                    )
                    if write_anisou and pa.anisou is not None:
                        cifanisou.append(
                            [anum, Element[pa.atomic_number].upper(), pa.name, pa.altloc, resname, res.chain, rnum,
                             pa.anisou[0], pa.anisou[1], pa.anisou[2], pa.anisou[3], pa.anisou[4], pa.anisou[5], '?',
                             '?', '?', '?', '?', '?', rnum, resname, res.chain, pa.name]
                        )
                    for key in sorted(others.keys()):
                        oatom = others[key]
                        x, y, z = oatom.xx, oatom.xy, oatom.xz
                        if renumber:
                            nmore += 1
                            anum = (pa.idx + 1 + nmore)
                        else:
                            anum = oatom.number or last_number + 1
                        last_number = anum
                        el = Element[oatom.atomic_number].upper()
                        cifatoms.append(
                            [atomrec, anum, el, oatom.name, oatom.altloc, resname, res.chain, '?', rnum,
                             res.insertion_code, x, y, z, oatom.occupancy, oatom.bfactor, '?', '?', '?', '?', '?', '',
                             rnum, resname, res.chain, oatom.name, '1']
                        )
                        if write_anisou and oatom.anisou is not None:
                            cifanisou.append(
                                [anum, Element[oatom.atomic_number].upper(), oatom.name, oatom.altloc, resname,
                                 res.chain, rnum, oatom.anisou[0], oatom.anisou[1], oatom.anisou[2],
                                 oatom.anisou[3], oatom.anisou[4], oatom.anisou[5], '?', '?', '?', '?', '?',
                                 '?', rnum, resname, res.chain, oatom.name]
                            )
        # Now write the PDBx file
        writer = PdbxWriter(dest)
        writer.write([cont])
        if own_handle:
            dest.close()

    @staticmethod
    def _process_struct_conn(struct, conn):
        """Adds bonds from the connect record to the structure

        May raise:
         IndexError (badly formatted rows that do not have the right number of elements)
        """
        atom_by_residue_identifiers_map = {
            (atom.residue.chain, atom.residue.name, atom.residue.number, atom.name): atom for atom in struct.atoms
        }

        chain_1_id = conn.getAttributeIndex("ptnr1_auth_asym_id")
        residue_1_name_id = conn.getAttributeIndex("ptnr1_auth_comp_id")
        residue_1_seq_id = conn.getAttributeIndex("ptnr1_auth_seq_id")
        atom_1_name_id = conn.getAttributeIndex("ptnr1_label_atom_id")

        chain_2_id = conn.getAttributeIndex("ptnr2_auth_asym_id")
        residue_2_name_id = conn.getAttributeIndex("ptnr2_auth_comp_id")
        residue_2_seq_id = conn.getAttributeIndex("ptnr2_auth_seq_id")
        atom_2_name_id = conn.getAttributeIndex("ptnr2_label_atom_id")

        order_id = conn.getAttributeIndex("pdbx_value_order")

        bond_type_order_map = {
            "sing": QualitativeBondType.SINGLE, "doub": QualitativeBondType.DOUBLE,
            "trip": QualitativeBondType.TRIPLE, "quad": QualitativeBondType.QUADRUPLE,
        }

        # NOTE - we are not currently using the *_auth_* keys that are present in the atom_site container.
        # They seem to be redundant, and the chain, residue name, and residue index seem to be enough (with the
        # atom name) to identify each atom.

        for i in range(conn.getRowCount()):
            row = conn.getRow(i)
            chain_1 = row[chain_1_id].strip().replace("?", "")
            residue_1_name = row[residue_1_name_id].strip()
            atom_1_name = row[atom_1_name_id].strip()
            residue_1_seq = try_convert(row[residue_1_seq_id].strip(), int)
            atom_1_key = (chain_1, residue_1_name, residue_1_seq, atom_1_name)


            chain_2 = row[chain_2_id].strip().replace("?", "")
            residue_2_name = row[residue_2_name_id].strip()
            atom_2_name = row[atom_2_name_id].strip()
            residue_2_seq = try_convert(row[residue_2_seq_id].strip(), int)
            atom_2_key = (chain_2, residue_2_name, residue_2_seq, atom_2_name)

            order = row[order_id].lower() if order_id != -1 else None

            qualitative_type = bond_type_order_map.get(order, None)

            a1 = atom_by_residue_identifiers_map.get(atom_1_key, None)
            a2 = atom_by_residue_identifiers_map.get(atom_2_key, None)

            if a1 is None or a2 is None:
                LOGGER.warning(
                    f"Cannot find 1 or more atoms in the bond {chain_1} {residue_1_name} {residue_1_seq} {atom_1_name}"
                    f" -- {chain_2} {residue_2_name} {residue_2_seq} {atom_2_name}"
                )
                continue
            struct.bonds.append(Bond(a1, a2, qualitative_type=qualitative_type))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _find_atom_index(struct, idx):
    """
    Returns the atom with the given index in the structure. This is required
    because atom indices may not start from 1 and may contain gaps in a PDB
    file. This tries to find atoms quickly, assuming that indices *do* start
    from 1 and have no gaps. It then looks up or down, depending on whether we
    hit an index too high or too low. So it *assumes* that the sequence is
    monotonically increasing. If the atom can't be found, None is returned
    """
    idx0 = min(max(idx - 1, 0), len(struct.atoms)-1)
    if struct[idx0].number == idx:
        return struct[idx0]
    if struct[idx0].number < idx:
        idx0 += 1
        while idx0 < len(struct.atoms):
            if struct[idx0].number == idx:
                return struct[idx0]
            idx0 += 1
        return None # not found
    else:
        idx0 -= 1
        while idx0 > 0:
            if struct[idx0].number == idx:
                return struct[idx0]
            idx0 -= 1
        return None # not found

def _needs_ter_card(res):
    """ Determines if a TER card is needed by seeing if the residue is a
    polymeric residue that is *not* bonded to the next residue
    """
    # First see if it's in the list of standard biomolecular residues. If so,
    # and it has no tail, no TER is needed
    std_resname = _standardize_resname(res.name)[0]
    if std_resname in StandardBiomolecularResidues:
        is_std_res = True
        if StandardBiomolecularResidues[std_resname].tail is None:
            return False
    else:
        is_std_res = False
    my_res_idx = res.idx
    residxs = set()
    for atom in res.atoms:
        for bond in atom.bonds:
            residxs |= {bond.atom1.residue.idx, bond.atom2.residue.idx}
    if my_res_idx + 1 in residxs:
        return False # It's connected to next residue
    elif is_std_res:
        return True
    else:
        # Heuristic -- add a TER if it's bonded to the previous residue, which
        # indicates it's polymeric. Otherwise don't.
        return my_res_idx - 1 in residxs

def _format_atom_name_for_pdb(atom):
    if len(atom.name) < 4 and len(Element[atom.atomic_number]) != 2:
        return ' %-3s' % atom.name
    return atom.name[:4]

def try_convert(value, cast_type, default=None):
    try:
        return cast_type(value)
    except ValueError:
        return default
