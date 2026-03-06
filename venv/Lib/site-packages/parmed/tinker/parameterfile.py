"""
This module contains classes regarding the Amoeba potential and loading in a
TINKER-based parameter file.
"""
import logging
import re
import warnings
from collections import OrderedDict
from ..exceptions import TinkerError, TinkerWarning

LOGGER = logging.getLogger(__name__)

class BookmarkedFile(object):
    """ Allows setting a bookmark and rewinding to that bookmark """

    def __init__(self, *args, **kwargs):
        self._stream = open(*args, **kwargs)

    def mark(self):
        self.bookmark = self.tell()

    def rewind(self):
        self.seek(self.bookmark)

    def __getattr__(self, attr):
        return getattr(self._stream, attr)

    def __del__(self):
        try:
            self._stream.close()
        except:
            LOGGER.exception("Unexpected error closing BookmarkedFile")

def _IS_INT(thing):
    try:
        int(thing)
        return True
    except ValueError:
        return False

def _IS_FLOAT(thing):
    try:
        float(thing)
        return True
    except ValueError:
        return False

class _ParamType:
    " All parameter types. This caches the list of parameters for easy access "
    TypeList = dict()
    _param_type = ''

    def __init__(self, *args, **kwargs):
        """ Instantiates the parameter type """
        raise NotImplementedError('virtual method')

    @classmethod
    def register(cls, obj, key):
        if key in cls.TypeList:
            warnmsg = f'Duplicate {cls._param_type} type found: {key}'
            if cls.TypeList[key] == obj:
                warnmsg += ' [same parameters]'
            else:
                warnmsg += ' [different parameters]'
            warnings.warn(warnmsg, TinkerWarning)
        cls.TypeList[key] = obj

    @classmethod
    def reset(cls):
        """
        Resets the cached list -- allows us to read in multiple parameter sets
        """
        cls.TypeList = dict()

    def __eq__(self, other):
        """ Make sure all attributes are the same """
        cls = type(self)
        if not isinstance(other, type(self)):
            raise TypeError(f'cannot compare type {type(self)} to type {type(other)}')
        for prop in dir(cls):
            if prop.startswith('_') or prop == 'TypeList':
                continue
            # Skip all callable attributes
            if hasattr(getattr(self, prop), '__call__'):
                continue
            if getattr(self, prop) != getattr(other, prop):
                return False

        return True

class _BondType(_ParamType):
    """ Bond parameter type """

    TypeList = dict()
    _param_type = 'bond'

    def __init__(self, idx1, idx2, k, req):
        idx1, idx2, = sorted((int(idx1), int(idx2)))
        self.k, self.req = float(k), float(req)
        key = f'{idx1}-{idx2}'
        self.register(self, key)

    def __repr__(self):
        return f'<_BondType: k={self.k}; req={self.req}>'

def get_angle_type(typecode, *args, **kwargs):
    """ Factory that returns the appropriate angle type """
    if typecode in 'fF':
        return _FourierAngleType(*args, **kwargs)
    return _AngleType(*args, **kwargs)

class _AngleType(_ParamType):

    TypeList = dict()
    _param_type = 'angle'

    def __init__(self, idx1, idx2, idx3, k, theteq, theteq2=None, theteq3=None):
        idx1, idx2, idx3 = int(idx1), int(idx2), int(idx3)
        idx1, idx3 = sorted((idx1, idx3))
        key = f'{idx1}-{idx2}-{idx3}'
        self.k, self.theteq = float(k), float(theteq)

        if theteq2 is not None:
            self.theteq2 = float(theteq2)
        else:
            self.theteq2 = None
        if theteq3 is not None:
            self.theteq3 = float(theteq3)
        else:
            self.theteq3 = None

        self.register(self, key)

    def __repr__(self):
        retval = f"<_AngleType: k={self.k}; theteq={self.theteq}"
        if self.theteq2 is not None:
            retval += f'; theteq2={self.theteq2}'
        if self.theteq3 is not None:
            retval += f'; theteq3={self.theteq3}'
        return retval + '>'

class _FourierAngleType(_AngleType):

    def __init__(self, idx1, idx2, idx3, k, theteq, periodicity):
        idx1, idx2, idx3 = int(idx1), int(idx2), int(idx3)
        idx1, idx3 = sorted((idx1, idx3))
        key = f'{idx1}-{idx2}-{idx3}'
        self.k = float(k)
        self.theteq = float(theteq)
        self.periodicity = float(periodicity)
        self.register(self, key)

    def __repr__(self):
        return f'<_FourierAngleType: k={self.k}; theteq={self.theteq}; periodicity={self.periodicity}>'

class _StretchBendType(_AngleType):
   
    TypeList = dict()
    _param_type = 'stretch-bend'

    def __init__(self, idx1, idx2, idx3, k1, k2):
        idx1, idx2, idx3 = int(idx1), int(idx2), int(idx3)
        idx1, idx3 = sorted((idx1, idx3))
        key = f'{idx1}-{idx2}-{idx3}'
        self.k1, self.k2 = float(k1), float(k2)
        self.register(self, key)

    def __repr__(self):
        return f'<_StretchBendType: k1={self.k1}; k2={self.k2}>'

class _UreyBradleyType(_AngleType):
   
    TypeList = dict()
    _param_type = 'urey-bradley'

    def __init__(self, idx1, idx2, idx3, k, req):
        idx1, idx2, idx3 = int(idx1), int(idx2), int(idx3)
        idx1, idx3 = sorted((idx1, idx3))
        key = f'{idx1}-{idx2}-{idx3}'
        self.k, self.req = float(k), float(req)
        self.register(self, key)

    def __repr__(self):
        return f'<_UreyBradleyType: k={self.k}; req={self.req}>'

class _OPBendType(_ParamType):

    TypeList = dict()
    _param_type = 'out-of-plane bending'

    def __init__(self, idx1, idx2, idx3, idx4, k):
        idx1, idx2, idx3, idx4 = int(idx1), int(idx2), int(idx3), int(idx4)
        idx3, idx4 = sorted((idx3, idx4))
        self.k = float(k)
        key = f'{idx1}-{idx2}-{idx3}-{idx4}'
        self.register(self, key)

    def __repr__(self):
        return f'<_OPBendType: k={self.k}>'

class _DihedralType(_ParamType):

    TypeList = dict()
    _param_type = 'dihedral'

    def __init__(self, idx1, idx2, idx3, idx4, *args):
        idx1, idx2, idx3, idx4 = int(idx1), int(idx2), int(idx3), int(idx4)
        if idx2 < idx3 or (idx2 == idx3 and (idx1 < idx4 or idx1 == idx4)):
            key = f'{idx1}-{idx2}-{idx3}-{idx4}'
        elif idx2 > idx3 or (idx2 == idx3 and idx1 > idx4):
            key = f'{idx4}-{idx3}-{idx2}-{idx1}'
        self.k, self.phase, self.periodicity = [], [], []
        for i in range(len(args)//3):
            self.k.append(float(args[i*3]))
            self.phase.append(float(args[i*3+1]))
            self.periodicity.append(float(args[i*3+2]))
        self.register(self, key)

    def __repr__(self):
        return f'<_DihedralType: k={self.k}; phase={self.phase}; per={self.periodicity}>'

class _PiTorsionType(_ParamType):

    TypeList = dict()
    _param_type = 'pi-torsion'

    def __init__(self, idx1, idx2, k):
        idx1, idx2 = sorted((int(idx1), int(idx2)))
        key = f'{idx1}-{idx2}'
        self.k = float(k)
        self.register(self, key)

    def __repr__(self):
        return f'<_PiTorsionType: k={self.k}>'

class _TorsionTorsionType(_ParamType):

    TypeList = dict()
    _param_type = 'torsion-torsion'

    def __init__(self, indexes, nx, ny):
        indexes = (int(i) for i in indexes)
        key = '{}-{}-{}-{}-{}'.format(*indexes)
        self.nx, self.ny = int(nx), int(ny)
        self.potential_grid = OrderedDict()
        self.register(self, key)

    def add_point(self, x, y, potential):
        self.potential_grid[(float(x), float(y))] = float(potential)

    def __repr__(self):
        return f'<_TorsionTorsion: {self.nx} x {self.ny} potential grid>'

class _MultipoleType(_ParamType):

    TypeList = dict()
    _param_type = 'multipole'

    def __init__(self, indexes, p1):
        indexes = (str(i) for i in indexes)
        key = '-'.join(indexes)
        self.potential_terms = [float(p1)]
        self.register(self, key)

    def __repr__(self):
        return f'<_MultipoleType: terms={self.potential_terms}>'

    def add_terms(self, terms):
        for term in terms:
            self.potential_terms.append(float(term))

def get_atom_type(index, atomic_number, mass, valence):
    """
    Factory for getting an _AtomType, but making sure that only one instance of
    a particular type is created
    """
    index = int(index)
    try:
        return _AtomType.TypeList[index]
    except KeyError:
        return _AtomType(index, atomic_number, mass, valence)

class _AtomType(object):
    """ An atom type """
    TypeList = dict() # All cached types -- ensures atom types are unique

    def __init__(self, index, atomic_number, mass, valence):
        self.index = int(index)
        self.atomic_number = int(atomic_number)
        self.mass = float(mass)
        self.valence = int(valence)
        _AtomType.TypeList[index] = self # cache this type

    @staticmethod
    def set_vdw_params(index, size, epsilon, reduction=None):
        inst = _AtomType.TypeList[int(index)]
        inst.size = float(size)
        inst.epsilon = float(epsilon)
        if reduction is not None:
            inst.reduction = float(reduction)
        else:
            inst.reduction = None

    def __repr__(self):
        retval = '<_AtomType: idx=%d; elem=%d; mass=%s; val=%d' % (
                       self.index, self.atomic_number, self.mass, self.valence)
        if hasattr(self, 'size'):
            retval += '; size=%s; eps=%s; red=%s' % (self.size, self.epsilon,
                self.reduction)
        return retval + '>'

    @classmethod
    def reset(cls):
        cls.TypeList = dict()

class _Atom(object):
    """ An atom in a parameter set """
    AtomList=dict()
    def __init__(self, index,typeindex, name, descrip, atomic_number, mass, val):
        self.name = name
        self.description = descrip
        self.type = get_atom_type(typeindex, atomic_number, mass, val)
        _Atom.AtomList[index] = self # cache this type

    # Allow _Atom instances to access (but not modify) type properties
    def _typeindex(self): return self.type.index
    def _atomic_number(self): return self.type.atomic_number
    def _element(self): return self.type.atomic_number
    def _valence(self): return self.type.valence
    def _size(self): return self.type.size
    def _epsilon(self): return self.type.epsilon
    def _reduction(self): return self.type.reduction
    def _blocked(self): raise NotImplementedError('Cannot set this attribute')
    # Now set the above as properties
    typeindex = property(fget=_typeindex, fset=_blocked)
    atomic_number = property(fget=_atomic_number, fset=_blocked)
    element = property(fget=_element, fset=_blocked)
    valence = property(fget=_valence, fset=_blocked)
    size = property(fget=_size, fset=_blocked)
    epsilon = property(fget=_epsilon, fset=_blocked)
    reduction = property(fget=_reduction, fset=_blocked)

    def set_polarizability(self, polarizability, thole, connected_types):
        self.polarizability = float(polarizability)
        self.thole = float(thole)
        self.connected_types = [int(i) for i in connected_types]

    def __repr__(self):
        retstr = f'<_Atom "{self.description}": name={self.name}; type={self.typeindex}'
        if hasattr(self, 'polarizability'):
            retstr += f'; dipole pol={self.polarizability}; thole={self.thole}; connected atoms={self.connected_types}'
        return retstr+'>'

    @classmethod
    def reset(cls):
        cls.TypeList = dict()

#==============================================================================

def reset():
    """
    Resets all of the TypeList instances (without destroying the data inside
    them) so we can load multiple parameter sets
    """
    _BondType.reset()
    _AngleType.reset()
    _FourierAngleType.reset()
    _StretchBendType.reset()
    _UreyBradleyType.reset()
    _OPBendType.reset()
    _DihedralType.reset()
    _PiTorsionType.reset()
    _TorsionTorsionType.reset()
    _MultipoleType.reset()
    _AtomType.reset()

#==============================================================================

class AmoebaParameterSet(object):
    """
    Contains all of the parameters found in an Amoeba parameter file from TINKER
    """
    atomre = re.compile(r'atom *(\d+) *(\d+) *([A-Za-z\-\+\*0-9]+) *"(.+)" *'
                        r'(\d+) *(\d+\.\d+) *(\d+)', re.I)
    anglere = re.compile(r'angle([ 345fF])')

    def __init__(self, fname=None):
        self.atoms = dict()
        self.atom_types = _AtomType.TypeList # For easy access
        self.bonds = _BondType.TypeList
        self.angles = _AngleType.TypeList # includes FourierAngleTypes
        self.stretch_bends = _StretchBendType.TypeList
        self.urey_bradleys = _UreyBradleyType.TypeList
        self.opbends = _OPBendType.TypeList
        self.dihedrals = _DihedralType.TypeList
        self.torsion_torsions = _TorsionTorsionType.TypeList
        self.multipoles = _MultipoleType.TypeList
        if fname is not None:
            self.load_parameter_file(fname)

    def load_parameter_file(self, fname):
        """
        Parses a parameter file and loads all of the parameters found into data
        structures.
        """
        self.attributes = dict()

        # First load the attributes from the header
        f = BookmarkedFile(fname, 'r')
        done_with_attributes = False
        line = f.readline().replace('\t', ' ')
        while line and not done_with_attributes:
            if 'Literature References' in line:
                done_with_attributes = True
                break
            line = (line + '#').strip()
            line = line[:line.index('#')].strip()
            if not line:
                line = f.readline().replace('\t', ' ')
                continue
            words = line.split()
            if len(words) != 2: continue
            # Now extract the property
            if _IS_INT(words[1]):
                self.attributes[words[0].lower()] = int(words[1])
            elif _IS_FLOAT(words[1]):
                self.attributes[words[0].lower()] = float(words[1])
            else:
                self.attributes[words[0].lower()] = words[1]
            line = f.readline().replace('\t', ' ')
        if not done_with_attributes:
            raise TinkerError('Could not find force field attributes.')
        # Now get the atom types
        while line.lstrip()[:5].lower() != 'atom ':
            line = f.readline().replace('\t', ' ')
        # Now loop through all atoms
        while line.lstrip()[:5].lower() == 'atom ':
            rematch = self.atomre.match(line)
            num, typenum, name, descrip, anum, mass, val = rematch.groups()
            self.atoms[int(num)] = _Atom(int(num),typenum, name, descrip, anum, mass, val)
            line = f.readline().replace('\t', ' ')
        # Now parse out the van der waals terms
        while line.lstrip()[:4].lower() != 'vdw ':
            line = f.readline().replace('\t', ' ')
        while line.lstrip()[:4].lower() == 'vdw ':
            _AtomType.set_vdw_params(*line.split()[1:])
            line = f.readline().replace('\t', ' ')
        # Now parse out the bonds
        while line.lstrip()[:5].lower() != 'bond ':
            line = f.readline().replace('\t', ' ')
        while line.lstrip()[:5].lower() == 'bond ':
            _BondType(*line.split()[1:])
            line = f.readline().replace('\t', ' ')
        # Now parse out the angles. Handle iring and Fourier terms
        rematch = self.anglere.match(line)
        while not rematch:
            line = f.readline().replace('\t', ' ')
            rematch = self.anglere.match(line)
        while rematch:
            try:
                get_angle_type(rematch.groups()[0], *line.split()[1:])
            except TypeError:
                LOGGER.debug('%s, %s', repr(rematch.groups()[0]), line.split()[1:])
                raise
            line = f.readline().replace('\t', ' ')
            rematch = self.anglere.match(line)
        # Now parse out the stretch-bend parameters. From here on out, some of
        # the terms may not exist in all versions of the force field, so
        # protect for EOF and make sure we rewind to avoid missing any terms.
        f.mark()
        while line.lstrip()[:7].lower() != 'strbnd ' and line:
            line = f.readline().replace('\t', ' ')
        while line.lstrip()[:7].lower() == 'strbnd ' and line:
            f.mark()
            _StretchBendType(*line.split()[1:])
            line = f.readline().replace('\t', ' ')
        # Get the Urey-Bradley term(s)
        f.rewind(); line = f.readline().replace('\t', ' ')
        while line.lstrip()[:9].lower() != 'ureybrad ' and line:
            line = f.readline().replace('\t', ' ')
        while line.lstrip()[:9].lower() == 'ureybrad ' and line:
            f.mark()
            _UreyBradleyType(*line.split()[1:])
            line = f.readline().replace('\t', ' ')
        # Get the out-of-plane bending
        f.rewind(); line = f.readline().replace('\t', ' ')
        while line.lstrip()[:7].lower() != 'opbend ' and line:
            line = f.readline().replace('\t', ' ')
        while line.lstrip()[:7].lower() == 'opbend ' and line:
            f.mark()
            _OPBendType(*line.split()[1:])
            line = f.readline().replace('\t', ' ')
        # Get the torsion parameters
        f.rewind(); line = f.readline().replace('\t', ' ')
        while line.lstrip()[:8].lower() != 'torsion ' and line:
            line = f.readline().replace('\t', ' ')
        while line.lstrip()[:8].lower() == 'torsion ' and line:
            f.mark()
            _DihedralType(*line.split()[1:])
            line = f.readline().replace('\t', ' ')
        # Get the pitorsions
        f.rewind(); line = f.readline().replace('\t', ' ')
        while line.lstrip()[:7] != 'pitors ' and line:
            line = f.readline().replace('\t', ' ')
        while line.lstrip()[:7] == 'pitors ' and line:
            f.mark()
            _PiTorsionType(*line.split()[1:])
            line = f.readline().replace('\t', ' ')
        # Get the coupled torsions
        f.rewind(); line = f.readline().replace('\t', ' ')
        while line.lstrip()[:8] != 'tortors ' and line:
            line = f.readline().replace('\t', ' ')
        while line.lstrip()[:8] == 'tortors ' and line:
            words = line.split()
            tortor = _TorsionTorsionType(words[1:6], words[6], words[7])
            line = f.readline().replace('\t', ' ')
            while line.strip():
                tortor.add_point(*line.split())
                line = f.readline().replace('\t', ' ')
            line = f.readline().replace('\t', ' ')
            f.mark()
        # Get the multipole terms
        f.rewind(); line = f.readline().replace('\t', ' ')
        while line.lstrip()[:10].lower() != 'multipole ' and line:
            line = f.readline().replace('\t', ' ')
        while line.lstrip()[:10].lower() == 'multipole ' and line:
            words = line.split()
            multipole = _MultipoleType(words[1:-1], words[-1])
            multipole.add_terms(f.readline().split())
            multipole.add_terms(f.readline().split())
            multipole.add_terms(f.readline().split())
            multipole.add_terms(f.readline().split())
            line = f.readline().replace('\t', ' ')
            f.mark()
        # Get the dipole polarizabilities
        f.rewind(); line = f.readline().replace('\t', ' ')
        while line.lstrip()[:9] != 'polarize ' and line:
            line = f.readline().replace('\t', ' ')
        while line.lstrip()[:9] == 'polarize ' and line:
            words = line.split()
            index = int(words[1])
            try:
                self.atoms[index].set_polarizability(
                            words[2], words[3], words[4:]
                )
            except IndexError:
                self.atoms[index].set_polarizability(words[2], words[3], [])
            line = f.readline().replace('\t', ' ')
        f.close()
        # Now clean up so we can load another parameter set
        reset()
        return
