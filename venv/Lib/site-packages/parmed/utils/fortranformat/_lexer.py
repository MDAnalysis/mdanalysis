from ._edit_descriptors import *
from ._exceptions import *

# Some lexer tokens to look out for
DIGITS = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
SIGNS = ['+', '-']
COMMA = [',']
DOT = ['.']
WHITESPACE = [' ', '\t', '\n']
QUOTE_CHARS = ['"', "'"]
DOUBLE_EDIT_DESCRIPTORS = ['EN', 'ES', 'TR', 'TL', 'BN', 'BZ', 'SP', 'SS']
SINGLE_EDIT_DESCRIPTORS = ['A', 'B', 'D', 'E', 'F', 'G', 'I', 'L', 'O', 'P', 'S', 'T', 'X', 'Z', ':', '/']
H_EDIT_DESCRIPTOR = ['H']
LEFT_PARENS = ['(']
RIGHT_PARENS = [')']
COLON = ':'
SLASH = '/'


def lexer(format):
    '''Lex the FORTRAN format statement into tokens'''
    tokens = []
    s = -1
    h_chars = None
    while True:
        # Get the next set of characters
        s = s + 1
        c0, c1, c2 = _get_chars(format, s)
        # If at end of format, end it all - aieee!
        if c0 is None:
            break
        # Read in an H edit descriptor string
        elif h_chars is not None:
            buff = format[s:s+h_chars]
            tokens.append(Token('QUOTED_STRING', buff))
            s = s + (h_chars - 1)
            h_chars = None
        # Skip whitespace
        elif c0 in WHITESPACE:
            continue
        # Read in a quoted string
        elif c0 in QUOTE_CHARS:
            buff = ''
            delim = c0
            while True:
                s = s + 1
                c0, c1, c2 = _get_chars(format, s)
                # Check if an escaped delimiter
                if (c0 == delim) and (c1 == delim):
                    s = s + 1
                    buff = buff + delim
                elif (c0 == delim):
                    break
                elif c0 is None:
                    # Premature end of format
                    raise InvalidFormat('Premature end of quoted string in format')
                else:
                    buff = buff + c0
            tokens.append(Token('QUOTED_STRING', buff))
        # Read in an integer
        elif c0 in DIGITS + SIGNS:
            # Check sign followed by at least one digit
            if (c0 in SIGNS) and (c1 not in DIGITS):
                # TODO: Is whitesapce allowed between sign and digit?
                raise InvalidFormat("Orphaned sign '%s' with no digits at position %d" % (c0, s))
            buff = c0
            while True:
                s = s + 1
                c0, c1, c2 = _get_chars(format, s)
                if (c0 not in DIGITS) or (c0 is None):
                    break
                else:
                    buff = buff + c0
            s = s - 1
            val = int(buff)
            if buff[0] in SIGNS:
                tokens.append(Token('INT', val))
            elif val == 0:
                tokens.append(Token('UINT', val))
            else:
                tokens.append(Token('NZUINT', val))
        # Read in a comma
        elif c0 in COMMA:
            tokens.append(Token('COMMA', None))
        # Read in a dot
        elif c0 in DOT:
            tokens.append(Token('DOT', None))
        # Read in double lettered edit descriptors
        elif (c1 is not None) and ((c0 + c1).upper() in DOUBLE_EDIT_DESCRIPTORS):
            ed_type = _get_ed_type((c0 + c1).upper())
            tokens.append(Token(ed_type, (c0 + c1).upper()))
            s = s + 1
        # Read in an H edit descriptor
        elif c0.upper() in H_EDIT_DESCRIPTOR:
            if (len(tokens) > 0) and (tokens[-1].type in ('NZUINT', 'UINT')):
                h_chars = tokens[-1].value
                tokens = tokens[:-1]
            else:
                raise InvalidFormat("Missing H descriptor number argument at position %d" % s)
        # Read in single lettered edit descriptors
        elif c0.upper() in SINGLE_EDIT_DESCRIPTORS:
            ed_type = _get_ed_type(c0.upper())
            tokens.append(Token(ed_type, c0.upper()))
        # Read in left parens
        elif c0 in LEFT_PARENS:
            tokens.append(Token('LEFT_PARENS', None))
        # Read in right parens
        elif c0 in RIGHT_PARENS:
            tokens.append(Token('RIGHT_PARENS', None))
        else:
            raise InvalidFormat('Character %s not recognised at position %d' % (c0, s))
    return tokens

def _get_ed_type(ed_string):
    if ed_string in ED1:
        ed_type = 'ED1'
    elif ed_string in ED2:
        ed_type = 'ED2'
    elif ed_string in ED3:
        ed_type = 'ED3'
    elif ed_string in ED4:
        ed_type = 'ED4'
    elif ed_string in ED5:
        ed_type = 'ED5'
    elif ed_string in ED6:
        ed_type = 'ED6'
    elif ed_string in ED7:
        ed_type = 'ED7'
    elif ed_string in ED8:
        ed_type = 'ED8'
    elif ed_string in ED9:
        ed_type = 'ED9'
    elif ed_string in ED10:
        ed_type = 'ED10'
    else:
        ed_type = None
    return ed_type

def _get_chars(format, s):
    try:
        c0 = format[s]
    except IndexError:
        c0 = None
    try:
        c1 = format[s+1]
    except IndexError:
        c1 = None
    try:
        c2 = format[s+2]
    except IndexError:
        c2 = None
    return (c0, c1, c2)


class InvalidFormat(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class Token(object):
    def __init__(self, type, value):
        self.type = type
        self.value = value
    def __repr__(self):
        return "\n  Token: type=%s,\tvalue=%s" % (self.type, str(self.value))

# Do some testing when run as a module

#if __name__ == '__main__':
#    import doctest
#    import os
#    globs = {'lexer' : lexer}
#    # Need to normalize whitespace since pasting into VIM converts tabs to
#    # spaces
#    doctest.testfile(os.path.join('tests', 'lexer_test.txt'), \
#        globs=globs, optionflags=doctest.NORMALIZE_WHITESPACE)
