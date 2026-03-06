"""
This module contains basic information and functionality related to individual
residues in typical biopolymers.
"""
from abc import ABC

__all__ = ['AminoAcidResidue', 'RNAResidue', 'DNAResidue', 'ALA', 'ARG', 'ASN',
           'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'HYP', 'ILE', 'LEU', 'LYS',
           'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'DA', 'DT',
           'DG', 'DC', 'A', 'U', 'G', 'C', 'SOLVENT_NAMES', 'EXTRA_POINT_NAMES',
           'CATION_NAMES', 'ANION_NAMES', 'ALLION_NAMES']

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BiomolecularResidue(ABC):
    """ Base class for different classes of biopolymer residues """
    _all_residues_by_name = dict()
    _all_residues_by_abbr = dict()
    _all_residues_by_symbol = dict()
    all_residues = []

    def __init_(self, *args, **kwargs):
        super().__init__()

    @classmethod
    def get(cls, key):
        raise NotImplementedError()

    def __str__(self):
        return self.name

    @classmethod
    def has(cls, thing):
        """
        Determines if a particular BiomolecularResidue or residue name is
        present in this classification of biomolecular residues

        Parameters
        ----------
        thing : str or :class:`BiomolecularResidue`

        Returns
        -------
        contains : bool
            If the residue or residue name *is* of this type, True. Otherwise,
            False.
        """
        if isinstance(thing, BiomolecularResidue):
            return thing in cls.all_residues
        try:
            cls.get(thing)
        except KeyError:
            return False
        else:
            return True

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AminoAcidResidue(BiomolecularResidue):
    """
    An individual amino acid residue.

    Parameters
    ----------
    name : str
        The name of the residue
    abbr : str
        The 3-letter abbreviation of the amino acid residue
    symbol : str
        The 1-letter symbol of the amino acid
    aliases : list of str, optional
        A list of other abbreviations that *also* refer to this residue

    Raises
    ------
    ValueError
        If any aliases have the same abbreviation as *other*
    """
    _all_residues_by_name = dict()
    _all_residues_by_abbr = dict()
    _all_residues_by_symbol = dict()
    all_residues = []

    def __init__(self, name, abbr, symbol, aliases=None):
        self.name = name
        self.abbr = abbr
        self.symbol = symbol
        type(self)._all_residues_by_name[name.upper()] = self
        type(self)._all_residues_by_abbr[abbr.upper()] = self
        if symbol is not None:
            type(self)._all_residues_by_symbol[symbol.upper()] = self
        type(self).all_residues.append(self)
        if aliases is not None:
            for alias in aliases:
                alias = alias.upper()
                if alias in type(self)._all_residues_by_abbr:
                    raise ValueError('%s is already an abbreviation' % alias)
                type(self)._all_residues_by_abbr[alias] = self

    def __repr__(self):
        return '<Amino Acid Residue %s: %s [%s]>' % (self.name, self.abbr,
                self.symbol)

    @classmethod
    def get(cls, key, abbronly=False):
        """
        Gets the amino acid corresponding to either the residue name, 3-letter
        abbreviation or 1-letter symbol. It is case-insensitive.

        Parameters
        ----------
        key : str
            1-letter symbol, 3-letter abbreviation, or residue name
        abbronly : bool
            If True, only look for the 3-letter abbreviation (not the 1-letter
            symbol)

        Returns
        -------
        residue : :class:`AminoAcidResidue`
            The residue corresponding to the given key

        Raises
        ------
        KeyError if ``key`` is not a symbol, abbreviation, or case-insensitive
        name of an amino acid residue, or any of its abbreviations.
        """
        if len(key) == 1 and not abbronly:
            return cls._all_residues_by_symbol[key.upper()]
        if len(key) == 3:
            return cls._all_residues_by_abbr[key.upper()]
        # Handle C- and N-termini that may be prepended with C or N
        if len(key) == 4 and key[0].upper() in 'CN':
            return cls._all_residues_by_abbr[key[1:].upper()]
        return cls._all_residues_by_name[key.upper()]

ALA = AminoAcidResidue('Alanine', 'ALA', 'A')
ARG = AminoAcidResidue('Arginine', 'ARG', 'R')
ASN = AminoAcidResidue('Asparagine', 'ASN', 'N')
ASP = AminoAcidResidue('Aspartate' ,'ASP', 'D', ['ASH', 'AS4'])
CYS = AminoAcidResidue('Cysteine', 'CYS', 'C', ['CYM', 'CYX'])
GLU = AminoAcidResidue('Glutamate', 'GLU', 'E', ['GLH', 'GL4'])
GLN = AminoAcidResidue('Glutamine', 'GLN', 'Q')
GLY = AminoAcidResidue('Glycine', 'GLY', 'G')
HIS = AminoAcidResidue('Histidine', 'HIS', 'H', ['HIP', 'HIE', 'HID'])
HYP = AminoAcidResidue('Hydroxyproline', 'HYP', None)
ILE = AminoAcidResidue('Isoleucine', 'ILE', 'I')
LEU = AminoAcidResidue('Leucine', 'LEU', 'L')
LYS = AminoAcidResidue('Lysine', 'LYS', 'K', ['LYN'])
MET = AminoAcidResidue('Methionine', 'MET', 'M')
PHE = AminoAcidResidue('Phenylalanine', 'PHE', 'F')
PRO = AminoAcidResidue('Proline', 'PRO', 'P')
SER = AminoAcidResidue('Serine', 'SER', 'S')
THR = AminoAcidResidue('Threonine', 'THR', 'T')
TRP = AminoAcidResidue('Tryptophan', 'TRP', 'W')
TYR = AminoAcidResidue('Tyrosine', 'TYR', 'Y')
VAL = AminoAcidResidue('Valine', 'VAL', 'V')

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DNAResidue(BiomolecularResidue):
    """ An individual DNA residue

    Parameters
    ----------
    name : str
        The name of the residue
    abbr : str
        The abbreviation of the nucleic acid residue
    aliases : list of str, optional
        A list of other abbreviations that *also* refer to this residue
    """
    _all_residues_by_name = dict()
    _all_residues_by_abbr = dict()
    _all_residues_by_symbol = dict()
    all_residues = []

    def __init__(self, name, abbr, aliases=None):
        self.name = name
        self.abbr = abbr
        type(self)._all_residues_by_name[name.upper()] = self
        type(self)._all_residues_by_abbr[abbr.upper()] = self
        type(self).all_residues.append(self)
        if aliases is not None:
            for alias in aliases:
                alias = alias.upper()
                if alias in type(self)._all_residues_by_abbr:
                    raise ValueError('%s is already an abbreviation' % alias)
                type(self)._all_residues_by_abbr[alias] = self

    def __repr__(self):
        return '<DNA Residue %s: %s>' % (self.name, self.abbr)

    @classmethod
    def get(cls, key):
        """
        Gets the nucleic acid corresponding to either the residue name or
        abbreviation. It is case-insensitive.

        Parameters
        ----------
        key : str
            abbreviation or residue name

        Returns
        -------
        residue : :class:`DNAResidue`
            The residue corresponding to the given key

        Raises
        ------
        KeyError if ``key`` is not a recognized residue name or abbreviation for
        an DNA residue.
        """
        try:
            if key[-1] in '35':
                return cls._all_residues_by_abbr[key[:-1].upper()]
            return cls._all_residues_by_abbr[key.upper()]
        except KeyError:
            return cls._all_residues_by_name[key.upper()]

class RNAResidue(DNAResidue):
    """ An individual RNA residue

    Parameters
    ----------
    name : str
        The name of the residue
    abbr : str
        The abbreviation of the nucleic acid residue
    aliases : list of str, optional
        A list of other abbreviations that *also* refer to this residue
    """
    _all_residues_by_name = dict()
    _all_residues_by_abbr = dict()
    _all_residues_by_symbol = dict()
    all_residues = []

    def __repr__(self):
        return '<RNA Residue %s: %s>' % (self.name, self.abbr)

    @classmethod
    def get(cls, key):
        """
        Gets the nucleic acid corresponding to either the residue name or
        abbreviation. It is case-insensitive.

        Parameters
        ----------
        key : str
            abbreviation or residue name

        Returns
        -------
        residue : :class:`RNAResidue`
            The residue corresponding to the given key

        Raises
        ------
        KeyError if ``key`` is not a recognized residue name or abbreviation for
        an RNA residue.
        """
        try:
            if key[-1] in '35':
                return cls._all_residues_by_abbr[key[:-1].upper()]
            return cls._all_residues_by_abbr[key.upper()]
        except KeyError:
            return cls._all_residues_by_name[key.upper()]

DG = DNAResidue('Guanine', 'DG', ['GUA', 'DG5', 'DG3', 'DGN'])
DC = DNAResidue('Cytosine', 'DC', ['CYT', 'DC5', 'DC3', 'DCN', 'DCP'])
DA = DNAResidue('Adenine', 'DA', ['ADE', 'DA5', 'DA3', 'DAN', 'DAP'])
DT = DNAResidue('Thymine', 'DT', ['THY', 'DT5', 'DT3'])
G = RNAResidue('Guanine', 'G', ['GUA', 'G5', 'G3', 'GN', 'RG', 'RG3', 'RG5', 'RGN'])
C = RNAResidue('Cytosine', 'C', ['CYT', 'CP', 'C5', 'C3', 'CN', 'RC', 'RC5', 'RC3', 'RCN'])
A = RNAResidue('Adenine', 'A', ['ADE', 'AP', 'A5', 'A3', 'AN', 'RA', 'RA3', 'RA5'])
U = RNAResidue('Uracil', 'U', ['URA', 'U3', 'U5', 'UN', 'RU', 'RU3', 'RU5', 'RUN'])
T = RNAResidue('Thymine', 'T', ['THY', 'T3', 'T5', 'TN', 'RT', 'RT3', 'RT5', 'RTN'])

WATER_NAMES = {'WAT', 'HOH', 'TIP3', 'TIP4', 'TIP5', 'SPCE', 'SPC', 'SWM4', 'SWM6'}
SOLVENT_NAMES = WATER_NAMES | {'SOL'}
EXTRA_POINT_NAMES = {'EP', 'LP'}
CATION_NAMES = {'Na+', 'Li+', 'Mg+', 'Rb+', 'MG', 'Cs+', 'POT', 'SOD', 'MG2',
                'CAL', 'RUB', 'LIT', 'ZN2', 'CD2', 'NA', 'K+', 'K', 'NA+'}
ANION_NAMES = {'Cl-', 'Br-', 'F-', 'I-', 'CLA', 'CL', 'BR', 'CL-'}
ALLION_NAMES = CATION_NAMES | ANION_NAMES
