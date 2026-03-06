from ._exceptions import *

def get_edit_descriptor_obj(name):
    '''Returns a new object instance from a string'''
    name = name.upper()
    if name == 'A':
        return A()
    elif name == 'B':
        return B()
    elif name == 'BN':
        return BN()
    elif name == 'BZ':
        return BZ()
    elif name == ':':
        return Colon()
    elif name == 'D':
        return D()
    elif name == 'E':
        return E()
    elif name == 'EN':
        return EN()
    elif name == 'ES':
        return ES()
    elif name == 'F':
        return F()
    elif name == 'G':
        return G()
    elif name == 'H':
        return H()
    elif name == 'I':
        return I()
    elif name == 'L':
        return L()
    elif name == 'O':
        return O()
    elif name == 'P':
        return P()
    elif name =='S':
        return S()
    elif name == '/':
        return Slash()
    elif name == 'SP':
        return SP()
    elif name == 'SS':
        return SS()
    elif name == 'T':
        return T()
    elif name == 'TL':
        return TL()
    elif name == 'TR':
        return TR()
    elif name == 'X':
        return X()
    elif name == 'Z':
        return Z()
    else:
        raise InvalidFormat('Expected an edit descriptor, got %s' % name)

# All the tokens defined in the F77 specification unless specified

class A(object):
    def __init__(self):
        self.repeat = None
        self.width = None
    def __repr__(self):
        return '<A repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + '>'

class QuotedString(object):
    def __init__(self, char_string=None):
        self.char_string = char_string
    def get_width(self):
        return len(self.char_string)
    width = property(get_width)
    def __repr__(self):
        return '<QuotedString char_string=' + str(self.char_string) + '>'

# Only in F95
class B(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.min_digits = None
    def __repr__(self):
        return '<B repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' min_digits=' + str(self.min_digits) + '>'
    
class BN(object):
    def __init__(self):
        pass
    def __repr__(self):
        return '<BN>'

class BZ(object):
    def __init__(self):
        pass
    def __repr__(self):
        return '<BZ>'

class Colon(object):
    def __init__(self):
        pass
    def __repr__(self):
        return '<Colon>'
    
class D(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.decimal_places = None
    def __repr__(self):
        return '<D repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' decimal_places=' + str(self.decimal_places) + '>'

class E(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.decimal_places = None
        self.exponent = None
    def __repr__(self):
        return '<E repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' decimal_places=' + str(self.decimal_places) + \
                ' exponent=' + str(self.exponent) + '>'
    
# Only in F95
class EN(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.decimal_places = None
        self.exponent = None
    def __repr__(self):
        return '<EN repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' decimal_places=' + str(self.decimal_places) + \
                ' exponent=' + str(self.exponent) + '>'

# Only in F95
class ES(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.decimal_places = None
        self.exponent = None
    def __repr__(self):
        return '<ES repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' decimal_places=' + str(self.decimal_places) + \
                ' exponent=' + str(self.exponent) + '>'

class F(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.decimal_places = None
    def __repr__(self):
        return '<F repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' decimal_places=' + str(self.decimal_places) + '>'
    
class FormatGroup(object):
    pass

class G(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.decimal_places = None
        self.exponent = None
    def __repr__(self):
        return '<G repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' decimal_places=' + str(self.decimal_places) + \
                ' exponent=' + str(self.exponent) + '>'

# Only in F77
class H(object):
    def __init__(self):
        self.num_chars = None
        self.char_string = None
    def __repr__(self):
        return '<H num_chars=' + str(self.num_chars) + \
                ' char_string=' + str(self.char_string) + '>'
    
class I(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.min_digits = None
    def __repr__(self):
        return '<I repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' min_digits=' + str(self.min_digits) + '>'
    
class L(object):
    def __init__(self):
        self.repeat = None
        self.width = None
    def __repr__(self):
        return '<L repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + '>'

# Only in F95
class O(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.min_digits = None
    def __repr__(self):
        return '<O repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' min_digits=' + str(self.min_digits) + '>'

class P(object):
    def __init__(self):
        self.scale = None
    def __repr__(self):
        return '<P scale=' + str(self.scale) + '>'
    
class S(object):
    def __init__(self):
        pass
    def __repr__(self):
        return '<S>'
    
class Slash(object):
    def __init__(self):
        self.repeat = None
        pass
    def __repr__(self):
        return '<Slash repeat=' + str(self.repeat) + '>'
    
class SP(object):
    def __init__(self):
        pass
    def __repr__(self):
        return '<SP>'
    
class SS(object):
    def __init__(self):
        pass
    def __repr__(self):
        return '<SS>'
    
class T(object):
    def __init__(self):
        self.num_chars = None
    def __repr__(self):
        return '<T num_chars=' + str(self.num_chars) + '>'
    
class TL(object):
    def __init__(self):
        self.num_chars = None
    def __repr__(self):
        return '<TL num_chars=' + str(self.num_chars) + '>'
    
class TR(object):
    def __init__(self):
        self.num_chars = None
    def __repr__(self):
        return '<TR num_chars=' + str(self.num_chars) + '>'

class X(object):
    def __init__(self):
        self.num_chars = None
    def __repr__(self):
        return '<X num_chars=' + str(self.num_chars) + '>'

# Only in F95
class Z(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.min_digits = None
    def __repr__(self):
        return '<Z repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' min_digits=' + str(self.min_digits) + '>'

# Categorise the edit descriptors depnding on how they should be parsed

ED1 = ['BN', 'BZ', 'SP', 'SS', 'S'] # Of form X only
ED2 = ['X'] # Of form nX only
ED3 = ['T', 'TR', 'TL', 'L'] # Of form Xn only
ED4 = ['A'] # Of form X or Xn
ED5 = ['D', 'F'] # Of form Xn.m only
ED6 = ['B', 'I', 'O', 'Z'] # Of form Xn or Xn.m
ED7 = ['E', 'EN', 'ES', 'G'] # Of form Xn.m or Xn.mEe
ED8 = ['P'] # Of form kX only, where k is a signed integer, may omit comma if followed by Type 5 or 7 edit descriptor
ED9 = [':'] # Of form X only, may omit comma either side
ED10 = ['/'] # Of form X only, may omit following comma and leading comma if no repeat
REPEATABLE_EDS = ['L', 'A', 'D', 'F', 'B', 'I', 'O', 'Z', 'E', 'EN', 'ES', 'G', '/']
OUTPUT_EDS = (L, A, D, F, B, I, O, Z, E, EN, ES, G)
CONTROL_EDS = (BN, BZ, P, SP, SS, S, X, T, TR, TL, Colon, Slash)
NON_REVERSION_EDS = (P, S, SP, SS, BN, BZ)
ALL_ED = ED1 + ED2 + ED3 + ED4 + ED5 + ED6 + ED7 + ED8 + ED9 + ED10

