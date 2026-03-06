"""
Contains all of the class objects for the Tinker topology
"""
from collections import OrderedDict

class Atom(object):
    """ An atom in the system """

    def __init__(self, symbol, type_, class_, 
                 atomic_number, mass, valence, desc):
        self.symbol = str(symbol).strip()
        self.type = int(type_)
        self.class_ = int(class_)
        self.atomic_number = int(atomic_number)
        self.mass = float(mass)
        self.valence = int(valence)
        self.desc = str(desc).strip()

    def add_vdw(self, radius, epsilon, radius14, epsilon14, reduction):
        """ Adds van der Waals terms to the atom """
        self.radius = float(radius)
        self.epsilon = float(epsilon)
        # These values might not exist
        if not str(radius14).strip():
            self.radius14 = None
        else:
            self.radius14 = float(radius14)
        if not str(epsilon14).strip():
            self.epsilon14 = None
        else:
            self.epsilon14 = float(epsilon14)
        if not str(reduction).strip():
            self.reduction = None
        else:
            self.reduction = float(reduction)

    def __repr__(self):
        return f'<Atom {self.symbol} [#{self.atomic_number}]; Mass={self.mass:.2f}; "{self.desc}">'

    def __str__(self):
        retstr =  (f'Atom {self.symbol} [#{self.atomic_number}] {self.desc}\n'
                   f'   Type:        {self.type}\n'
                   f'   Class:       {self.class_}\n'
                   f'   Mass:        {self.mass:.2f}\n'
                   f'   Valence:     {self.valence}')

        if hasattr(self, 'radius'):
            retstr += (f'\n   vdW Rad:     {self.radius:.2f}\n'
                       f'   vdW Eps:     {self.epsilon:.2f}')
            if self.radius14 is not None or self.epsilon14 is not None:
                retstr += f'\n   vdW 1-4 Rad: {self.radius14:.2f}\n   vdW 1-4 Eps: {self.epsilon14:.2f}'

        if self.reduction is not None:
            retstr += f'\n   vdW Reduct:  {self.reduction:.2f}'

        return retstr

class _ParamTypeList(list):
    """
    A list for all types of parameters. This should be subclassed for each
    parameter type that's necessary. It performs type checking as well as other
    tasks.
    """

    typeclass = None

    def __init__(self):
        super().__init__()

    def add(self, *args, **kwargs):
        self.append(self.typeclass(*args, **kwargs))

    def append(self, thing):
        if not isinstance(thing, self.typeclass):
            raise TypeError(f'Can only append "{self.typeclass.__name__}" objects to {self.__class__.__name__}')
        super().append(thing)

    def extend(self, things):
        for thing in things:
            self.append(thing)

class AtomList(_ParamTypeList):
    """ List for Atom instances """
    typeclass = Atom

class BondStretch(object):
    """ A bond-stretching term """
    def __init__(self, atom1, atom2, k, req):
        self.atom1 = atom1
        self.atom2 = atom2
        self.k = float(k)
        self.req = float(req)

    def __repr__(self):
        return f'<{self.__class__.__name__} [{self.atom1.symbol}-{self.atom2.symbol}]; k={self.k:.2f}; Req={self.req:.2f}>'

    def __str__(self):
        return (f'{self.__class__.name__} {repr(self.atom1)} --- {repr(self.atom2)}\n'
                f'     K = {self.k:.4f}\n'
                f'   Req = {self.req:.4f}')

class BondList(_ParamTypeList):
    typeclass = BondStretch

class AngleBend(object):
    """ An angle-bending term """
    def __init__(self, atom1, atom2, atom3, k, theteq, fold, type_):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.k = float(k)
        self.theteq = float(theteq)
        if str(fold).strip():
            self.fold = float(fold)
        else:
            self.fold = None
        self.type = type_.strip() or 'Harmonic'

    def __repr__(self):
        retstr = (
            f'<AngleBend [{self.atom1.symbol}-{self.atom2.symbol}-{self.atom3.symbol}]; '
            f'k={self.k:.2f}; Theta_eq={self.theteq:.2f}'
        )
        if self.fold is not None:
            retstr += f'; fold={self.fold:.1f}'
        return retstr + f'; {self.type}>'

    def __str__(self):
        retstr = (f'AngleBend {repr(self.atom1)} --- {repr(self.atom2)} --- {repr(self.atom3)}\n'
                  f'     K = {self.k:.4f}\n'
                  f' Theta = {self.theteq:.4f}\n')
        if self.fold is not None:
            retstr += f'  Fold = {self.fold:.1f}\n'
        return retstr + f'  Type = {self.type}'

class AngleList(_ParamTypeList):
    typeclass = AngleBend

class StretchBend(object):
    """ A stretch-bending term """
    def __init__(self, atom1, atom2, atom3, k, theteq, r1eq, r2eq):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.k = float(k)
        self.theteq = float(theteq)
        self.r1eq = float(r1eq)
        self.r2eq = float(r2eq)

    def __repr__(self):
        return (
            f'<StretchBend [{self.atom1.symbol}-{self.atom2.symbol}-{self.atom3.symbol}]; '
            f'k={self.k:.2f}; Theta_eq={self.theteq:.2f}; R1eq={self.r1eq:.2f}; R2eq={self.r2eq:.2f}>'
        )

    def __str__(self):
        return (f'StretchBend {repr(self.atom1)} --- {repr(self.atom2)} --- {repr(self.atom3)}\n'
                f'     K = {self.k:.4f}\n'
                f' Theta = {self.theteq:.4f}\n'
                f'  R1eq = {self.r1eq:.4f}\n'
                f'  R2eq = {self.r2eq:.4f}')

class StretchBendList(_ParamTypeList):
    typeclass = StretchBend

class UreyBradley(BondStretch):
    """ A urey-bradley term -- functionally identical to a bond-stretch """

class UreyBradleyList(_ParamTypeList):
    typeclass = UreyBradley

class OutOfPlaneBend:
    """ An out-of-plane bending term """
    def __init__(self, atom1, atom2, atom3, atom4, k):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.k = float(k)

    def __repr__(self):
        return f'<{self.__class__.__name__} [{self.atom1.symbol}-{self.atom2.symbol}-{self.atom3.symbol}-{self.atom4.symbol}]; k={self.k:.2f}>'

    def __str__(self):
        return (f'{self.__class__.__name__} {repr(self.atom1)} --- {repr(self.atom2)} --- {repr(self.atom3)} --- {repr(self.atom4)}\n'
                f'   K = {self.k:.2f}')

class OutOfPlaneBendList(_ParamTypeList):
    typeclass = OutOfPlaneBend

class OutOfPlaneDist(OutOfPlaneBend):
    """ An out-of-plane distance (functionally equivalent to OOP Bend """

class OutOfPlaneDistList(_ParamTypeList):
    typeclass = OutOfPlaneDist

class TorsionAngle(object):
    """ Torsional Angle parameter """
    def __init__(self, atom1, atom2, atom3, atom4, args):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        # We don't know how many terms to expect, but it will be a multiple of 3
        if len(args) % 3 != 0:
            raise TypeError('TorsionAngle expects an equal number of '
                            'Amplitudes, phases, and periodicities.')
        nterms = len(args) // 3
        self.amplitude = tuple([float(args[3*i]) for i in range(nterms)])
        self.phase = tuple([float(args[3*i+1]) for i in range(nterms)])
        self.periodicity = tuple([int(args[3*i+2]) for i in range(nterms)])
        if (len(self.amplitude) != len(self.phase) or
                len(self.amplitude) != len(self.periodicity)):
            raise RuntimeError('BUGBUG!! Inconsistent # of terms in torsion')

    def __repr__(self):
        if len(self.amplitude) == 0:
            return f'<TorsionAngle [{self.atom1.symbol}-{self.atom1.symbol}-{self.atom1.symbol}-{self.atom1.symbol}]; No Terms>'
        return (f'<TorsionAngle [{self.atom1.symbol}-{self.atom1.symbol}-{self.atom1.symbol}-{self.atom1.symbol}]; '
                f'Ampl={self.amplitude}; Phase={self.phase}; Per={self.periodicity}>')

    def __str__(self):
        if len(self.amplitude) == 0:
            return (f'TorsionAngle {repr(self.atom1)} --- {repr(self.atom2)} --- {repr(self.atom3)} --- {repr(self.atom4)}\n'
                    '   No Terms')
        retstr = f'TorsionAngle {repr(self.atom1)} --- {repr(self.atom2)} --- {repr(self.atom3)} --- {repr(self.atom4)}'
        seq = range(len(self.amplitude))
        for i, amp, phase, per in enumerate(zip(seq, self.amplitude, self.phase, self.periodicity)):
            retstr += (f'\n   Term {i + 1}\n'
                       f'   Amplitude = {amp:.4f}\n'
                       f'       Phase = {phase:.4f}\n'
                       f' Periodicity = {per:.4f}')
        return retstr

class TorsionAngleList(_ParamTypeList):
    typeclass = TorsionAngle

class PiTorsion(object):
    """ A Pi-Orbital Torsion parameter """
    def __init__(self, atom1, atom2, amplitude):
        self.atom1, self.atom2 = atom1, atom2
        self.amplitude = float(amplitude)

    def __repr__(self):
        return f'<PiTorsion [{self.atom1.symbol}-{self.atom2.symbol}]; Ampl={self.amplitude:.2f}>'

    def __str__(self):
        return f'PiTorsion {repr(self.atom1)} --- {repr(self.atom1)}\n   Amplitude = {self.amplitude:.4f}'

class PiTorsionList(_ParamTypeList):
    typeclass = PiTorsion

class TorsionTorsion(object):
    """ A coupled-torsion parameter (like CMAP) """
    def __init__(self, atom1, atom2, atom3, atom4, atom5, spline1, spline2):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.atom5 = atom5
        self.spline1 = int(spline1)
        self.spline2 = int(spline2)

    def __repr__(self):
        return (
            f'<TorsionTorsion [{self.atom1.symbol}-{self.atom2.symbol}-{self.atom3.symbol}-{self.atom4.symbol}-{self.atom5.symbol}]; '
            f'Spline=({self.spline1},{self.spline2})>'
        )

    def __str__(self):
        return (f'TorsionTorsion {repr(self.atom1)} --- {repr(self.atom2)} --- {repr(self.atom3)} ---\n'
                f'                      {repr(self.atom4)} --- {repr(self.atom5)}\n'
                f'   Spline Grid = ({self.spline1} x {self.spline2})')

class TorsionTorsionList(_ParamTypeList):
    typeclass = TorsionTorsion

class TorsionTorsionGrid(OrderedDict):
    """
    The interpolation grid of a coupled-torsion correction map. Every unique
    grid is cached and if a duplicate grid is instantiated, a reference to the
    original grid is returned. As a result, all unique TorsionTorsionGrid
    instances are singletons and should be compared for equality with "is"
    """
    # Because these grids are space-consuming, we only hold onto unique grids
    _typelist = list()

    @classmethod
    def new(cls, data):
        inst = cls()
        for d in data:
            inst[tuple(d[:2])] = d[2]
        # Potentially expensive comparison of all grids.
        for ttg in TorsionTorsionGrid._typelist:
            if ttg == inst:
                return ttg
        TorsionTorsionGrid._typelist.append(inst)
        return inst
   
    def __eq__(self, other):
        if self.keys() != other.keys(): return False
        for key in self:
            if abs(self[key] - other[key]) > 1e-8: return False
        return True
   
    def __ne__(self, other):
        return not self == other

    # No comparisons are implemented
    def __gt__(self, other):
        raise NotImplemented('TorsionTorsionGrid instances are not well-ordered')
    __lt__ = __ge__ = __le__ = __gt__

class AtomicMultipole(object):
    """ Atomic multipole parameter """
    def __init__(self, atom, frame, definition, moments):
        self.atom = atom
        self.frame = [int(i) for i in frame if str(i).strip()]
        self.definition = definition.strip()
        self.moment = [float(x) for x in moments]

    def __repr__(self):
        return f'<AtomicMultipole [{self.atom.symbol}] {self.definition}; Frame={self.frame}; Moments={self.moment}>'

    def __str__(self):
        return (f'AtomicMultipole {repr(self.atom)} "{self.definition}"\n'
                f'     Frame = {self.frame}\n'
                f'   Moments = {self.moment[0]:8.5f}\n'
                f'             {self.moment[1]:8.5f} {self.moment[2]:8.5f} {self.moment[3]:8.5f}\n'
                f'             {self.moment[4]:8.5f}\n'
                f'             {self.moment[5]:8.5f} {self.moment[6]:8.5f}\n'
                f'             {self.moment[7]:8.5f} {self.moment[8]:8.5f} {self.moment[9]:8.5f}')

class AtomicMultipoleList(_ParamTypeList):
    typeclass = AtomicMultipole

class DipolePolarizability(object):
    """ A dipole polarizability parameter """
    def __init__(self, atom, alpha, group):
        self.atom = atom
        self.alpha = float(alpha)
        self.group = [int(g) for g in group]

    def __repr__(self):
        return f'<DipolePolarizability [{self.atom.symbol}]; Alpha={self.alpha:.2f}; Group={self.group}>'

    def __str__(self):
        return (f'DipolePolarizability {repr(self.atom)}\n'
                f'    Alpha = {self.alpha:.4f}\n'
                f'    Group = {self.group}')

class DipolePolarizabilityList(_ParamTypeList):
    typeclass = DipolePolarizability
