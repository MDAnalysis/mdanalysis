from ._lexer import Token
from ._edit_descriptors import *
from ._exceptions import *
from . import config

def parser(tokens, version=None):
    # Parse the full edit descriptors
    eds = _parse_tokens(tokens, reversion=False, version=None)
    # Parse the edit descriptors used for the reversion of format control 
    # (F95 format 12.2.2)
    reversion_eds = _parse_tokens(tokens, reversion=True, version=None)
    return eds, reversion_eds

def _parse_tokens(tokens, reversion=False, version=None):
    # Remove outer parens is there are any
    tokens = _remove_outer_parens(tokens)
    # Get only the reversion tokens
    if reversion == True:
        tokens = _get_reversion_tokens(tokens)
    # First expand the parentheses
    tokens = _expand_parens(tokens)
    # Split on commas
    token_sets = _split_on_commas(tokens)
    # Split on ED9 (i.e. :)
    token_sets = _split_on_ed9(token_sets)
    # Split the ED10 (i.e. /)
    token_sets = _split_on_ed10(token_sets)
    # Split the ED8 (i.e. P edit descriptors)
    token_sets = _split_on_ed8(token_sets)
    # Process each set of edit descriptors
    eds = []
    for token_set in token_sets:
        # Assume first edit descriptor is the one to process
        ed_type = None
        ed_value = None
        for token in token_set:
            if token.type in ['ED1', 'ED2', 'ED3', 'ED4', 'ED5', 'ED6', 'ED7', 'ED8', 'ED9', 'ED10', 'QUOTED_STRING']:
                ed_type = token.type
                ed_value = token.value
                break
        # TODO: something more responsible here ...
        if ed_type is None:
            continue
        # If repeatable and first token is repeat number then cache
        repeat = None
        if ed_value in REPEATABLE_EDS and (token_set[0].type in ['NZUINT', 'UINT']):
            repeat = token_set[0].value
            token_set = token_set[1:]
        # Process the edit descriptor
        if ed_type == 'QUOTED_STRING':
            ed = _read_quoted_string(token_set)
        elif ed_type == 'ED1':
            ed = _read_ed1(token_set)
        elif ed_type == 'ED2':
            ed = _read_ed2(token_set)
        elif ed_type == 'ED3':
            ed = _read_ed3(token_set)
        elif ed_type == 'ED4':
            ed = _read_ed4(token_set)
        elif ed_type == 'ED5':
            ed = _read_ed5(token_set)
        elif ed_type == 'ED6':
            ed = _read_ed6(token_set)
        elif ed_type == 'ED7':
            ed = _read_ed7(token_set)
        elif ed_type == 'ED8':
            ed = _read_ed8(token_set)
        elif ed_type == 'ED9':
            ed = _read_ed9(token_set)
        elif ed_type == 'ED10':
            ed = _read_ed10(token_set)
        else:
            raise InvalidFormat('Could not identify edit descriptor in sequence $s' % str(token_set))
        # If there is a repeat number cached, then apply
        if repeat is not None:
            ed.repeat = repeat
        # Add to the list
        eds.append(ed)
    return eds


# Functions that group the tokens into sets

def _expand_parens(tokens):
    new_tokens = []
    get_tokens = iter(tokens)
    for t0 in get_tokens:
        if t0.type != 'LEFT_PARENS':
            new_tokens.append(t0)
        else:
            # Read in all tokens in subformat and recurse back to self
            paren_tokens = []
            nesting = 1
            while nesting > 0:
                try:
                    t1 = next(get_tokens)
                except StopIteration:
                    raise InvalidFormat('Open parens in format')
                if t1.type == 'LEFT_PARENS':
                    nesting = nesting + 1
                elif t1.type == 'RIGHT_PARENS':
                    nesting = nesting - 1
                paren_tokens.append(t1)
            # Remove last right paren
            paren_tokens = paren_tokens[:-1]
            # If repeated, then repeat the tokens accordingly
            if (len(new_tokens) > 0) and (new_tokens[-1].type in ['NZUINT', 'UINT']):
                repeat = new_tokens[-1].value
                # Remove the repeat NZUINT, UINT
                new_tokens = new_tokens[:-1]
                new_tokens.extend(repeat * (_expand_parens(paren_tokens) + [Token('COMMA', None)]))
            else:
                new_tokens.extend(_expand_parens(paren_tokens))
    return new_tokens


def _split_on_commas(tokens):
    token_sets = []
    set_buff = []
    for t0 in tokens:
        if t0.type == 'COMMA':
            token_sets.append(set_buff)
            set_buff = []
        else:
            set_buff.append(t0)
    token_sets.append(set_buff)
    return token_sets


def _split_on_ed9(token_sets):
    '''Splits on :'''
    new_token_sets = []
    for token_set in token_sets:
        if 'ED9' not in [t.type for t in token_set]:
            new_token_sets.append(token_set)
        else:
            buff = []
            for token in token_set:
                if token.type == 'ED9':
                    if len(buff) > 0:
                        new_token_sets.append(buff)
                        buff = []
                    new_token_sets.append([token])
                else:
                    buff.append(token)
            if len(buff) > 0:
                new_token_sets.append([token])
    return new_token_sets


def _split_on_ed10(token_sets):
    '''Splits on /'''
    new_token_sets = []
    for token_set in token_sets:
        # May have a repeat on the slash if preceded by a comma
        if (len(token_set) > 2) and ((token_set[0].type in ['UINT', 'NZUINT']) and (token_set[1].type == 'ED10')):
            new_token_sets.append(token_set[:2])
            token_set = token_set[2:]
        buff = []
        for token in token_set:
            if token.type == 'ED10':
                if len(buff) > 0:
                    new_token_sets.append(buff)
                    buff = []
                new_token_sets.append([token])
            else:
                buff.append(token)
        if len(buff) > 0:
            new_token_sets.append(buff)
    return new_token_sets


def _split_on_ed8(token_sets):
    '''Splits on ED8 (i.e. P edit descriptors)'''
    new_token_sets = []
    for token_set in token_sets:
        # Append to new list if no ED8
        if 'ED8' not in [t.type for t in token_set]:
            new_token_sets.append(token_set)
        # Otherwise split string on ED9
        elif (token_set[0].type in ['INT', 'UINT', 'NZUINT']) and (token_set[1].type == 'ED8'):
            new_token_sets.append(token_set[:2])
            new_token_sets.append(token_set[2:])
        else:
            raise InvalidFormat('P edit descriptor in invalid position')
    return new_token_sets

# Function to extract only the tokens for the reversion of control

def _get_reversion_tokens(tokens):
    reversion_tokens = []
    # Easier to work backwards
    nesting = None
    for token in tokens[::-1]:
        # End of loop condition
        if (nesting is not None) and (nesting < 1):
            # Parens may have a repeat number
            if token.type in ['UINT', 'NZUINT']:
                reversion_tokens.append(token)
            break
        # Read till the first right parens
        if token.type == 'RIGHT_PARENS':
            if nesting is None:
                nesting = 1
            else:
                nesting = nesting + 1
        elif token.type == 'LEFT_PARENS':
            if nesting is None:
                raise InvalidFormat('Unbalanced parens in format')
            else:
                nesting = nesting - 1
        reversion_tokens.append(token)
    # Tokens are in reverse order
    reversion_tokens.reverse()
    return reversion_tokens

# The functions that read particular edit descriptors sequences

def _read_quoted_string(tokens):
    # Of form X only
    type_string = ",".join([t.type for t in tokens])
    if type_string != "QUOTED_STRING":
        raise InvalidFormat('Token %s has invalid neighbouring token' % tokens[0])
    ed = QuotedString()  
    ed.char_string = tokens[0].value
    return ed

def _read_ed1(tokens):
    # Of form X only
    type_string = ",".join([t.type for t in tokens])
    if type_string != "ED1":
        raise InvalidFormat('Token %s has invalid neighbouring token' % tokens[0])
    ed = get_edit_descriptor_obj(tokens[0].value)   
    return ed

def _read_ed2(tokens):
    # Of form nX only
    type_string = ",".join([t.type for t in tokens])
    if type_string != "NZUINT,ED2":
        raise InvalidFormat('Token %s has invalid neighbouring token' % tokens[0])
    ed = get_edit_descriptor_obj(tokens[1].value)
    ed.num_chars = tokens[0].value
    return ed

def _read_ed3(tokens):
    # Of form Xn only
    type_string = ",".join([t.type for t in tokens])
    if type_string != "ED3,NZUINT":
        raise InvalidFormat('Token %s has invalid neighbouring token' % tokens[0])
    ed = get_edit_descriptor_obj(tokens[0].value)
    # L edit descriptor has a width rather than num_chars
    if hasattr(ed, 'width'):
        ed.width = tokens[1].value
    else:
        ed.num_chars = tokens[1].value
    return ed

def _read_ed4(tokens):
    # Of form X or Xn
    type_string = ",".join([t.type for t in tokens])
    if type_string in ["ED4", "ED4,NZUINT"] or \
        (config.ALLOW_ZERO_WIDTH_EDS and (type_string == "ED4,UINT")):
        ed = get_edit_descriptor_obj(tokens[0].value)
        if len(tokens) > 1:
            ed.width = tokens[1].value
    else:
        raise InvalidFormat('Token %s has invalid neighbouring token' % tokens[0])
    return ed

def _read_ed5(tokens):
    # Of form Xn.m only
    type_string = ",".join([t.type for t in tokens])
    if type_string in ["ED5,NZUINT,DOT,UINT", "ED5,NZUINT,DOT,NZUINT"] or \
      (config.ALLOW_ZERO_WIDTH_EDS and (type_string in \
      ["ED5,UINT,DOT,UINT", "ED5,UINT,DOT,NZUINT"])):
        ed = get_edit_descriptor_obj(tokens[0].value)
        ed.width = tokens[1].value
        ed.decimal_places = tokens[3].value
    else:
        raise InvalidFormat('%s has invalid neighbouring token' % tokens[0])
    return ed

def _read_ed6(tokens):
    # Of form Xn or Xn.m
    type_string = ",".join([t.type for t in tokens])
    if type_string == "ED6,NZUINT" or \
        (config.ALLOW_ZERO_WIDTH_EDS and (type_string == "ED6,UINT")):
        ed = get_edit_descriptor_obj(tokens[0].value)
        ed.width = tokens[1].value
        ed.min_digits = None
    elif type_string in ["ED6,NZUINT,DOT,UINT", "ED6,NZUINT,DOT,NZUINT"] or \
      (config.ALLOW_ZERO_WIDTH_EDS and (type_string in \
      ["ED6,UINT,DOT,UINT", "ED6,UINT,DOT,NZUINT"])):
        ed = get_edit_descriptor_obj(tokens[0].value)
        ed.width = tokens[1].value
        ed.min_digits = tokens[3].value
    else:
        raise InvalidFormat('%s has invalid neighbouring token' % tokens[0])
    return ed

def _read_ed7(tokens):
    # Of form Xn.m or Xn.mEe
    type_string = ",".join([t.type for t in tokens])
    if type_string in ["ED7,NZUINT,DOT,UINT", "ED7,NZUINT,DOT,NZUINT"] or \
      (config.ALLOW_ZERO_WIDTH_EDS and (type_string in \
      ["ED7,UINT,DOT,UINT", "ED7,UINT,DOT,NZUINT"])):
        ed = get_edit_descriptor_obj(tokens[0].value)
        ed.width = tokens[1].value
        ed.decimal_places = tokens[3].value
        ed.exponent = None
    elif type_string in ['ED7,NZUINT,DOT,NZUINT,ED7,NZUINT', \
            'ED7,NZUINT,DOT,NZUINT,ED7,UINT', \
            'ED7,NZUINT,DOT,NZUINT,ED7,INT', \
            'ED7,NZUINT,DOT,UINT,ED7,NZUINT', \
            'ED7,NZUINT,DOT,UINT,ED7,UINT', \
            'ED7,NZUINT,DOT,UINT,ED7,INT'] or \
        (config.ALLOW_ZERO_WIDTH_EDS and (type_string in \
            ['ED7,UINT,DOT,NZUINT,ED7,NZUINT', \
            'ED7,UINT,DOT,NZUINT,ED7,UINT', \
            'ED7,UINT,DOT,NZUINT,ED7,INT', \
            'ED7,UINT,DOT,UINT,ED7,NZUINT', \
            'ED7,UINT,DOT,UINT,ED7,UINT', \
            'ED7,UINT,DOT,UINT,ED7,INT'])):
        ed = get_edit_descriptor_obj(tokens[0].value)
        ed.width = tokens[1].value
        ed.decimal_places = tokens[3].value
        ed.exponent = tokens[5].value
    else:
        raise InvalidFormat('%s has invalid neighbouring token' % tokens[0])
    return ed

def _read_ed8(tokens):
    # Of form kX only, where k is a signed integer, may omit comma if followed
    # by Type 5 or 7 edit descriptor
    type_string = ",".join([t.type for t in tokens])
    if type_string in ["NZUINT,ED8", "UINT,ED8", "INT,ED8"]:
        ed = get_edit_descriptor_obj(tokens[1].value)
        ed.scale = tokens[0].value
    else:
        raise InvalidFormat('%s has invalid neighbouring token' % tokens[0])
    return ed

def _read_ed9(tokens):
    # Of form X only, may omit comma either side
    type_string = ",".join([t.type for t in tokens])
    if type_string == "ED9":
        ed = get_edit_descriptor_obj(tokens[0].value)
    else:
        raise InvalidFormat('%s has invalid neighbouring token' % tokens[0])
    return ed

def _read_ed10(tokens):
    # Of form X only, may omit following comma and leading comma if no repeat
    type_string = ",".join([t.type for t in tokens])
    if type_string == "ED10":
        ed = get_edit_descriptor_obj(tokens[0].value)
    else:
        raise InvalidFormat('%s has invalid neighbouring token' % tokens[0])
    return ed


# Functions that pre-process the token list

def _remove_outer_parens(tokens):
    # Finally, remove outer parens is there are any
    if (tokens[0].type == 'LEFT_PARENS') and (tokens[-1].type == 'RIGHT_PARENS'):
        tokens = tokens[1:-1]
    return tokens


# Run some tests if run as a script

#if __name__ == '__main__':
#    import doctest
#    import os
#    from _lexer import lexer
#    globs = {'lexer' : lexer, 'parser' : parser}
#    # Need to normalize whitespace since pasting into VIM converts tabs to
#    # spaces
#    doctest.testfile(os.path.join('tests', 'parser_test.txt'), \
#        globs=globs, optionflags=doctest.NORMALIZE_WHITESPACE)
