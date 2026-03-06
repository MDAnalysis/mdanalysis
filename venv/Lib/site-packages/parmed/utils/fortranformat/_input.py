import re
import pdb
from ._edit_descriptors import *
from ._misc import expand_edit_descriptors, has_next_iterator
from . import config

WIDTH_OPTIONAL_EDS = [A]
NON_WIDTH_EDS = [BN, BZ, P, SP, SS, S, X, T, TR, TL, Colon, Slash]
FORBIDDEN_EDS = [QuotedString, H]

# Some problems without pre written input vars:
#   Cannot say when reversion conditions are met
#   Cannot determine width of A edit descriptor
#   Cannot determine complex input
#   Cannot determine proper input for G edit descriptors


def input(eds, reversion_eds, records, num_vals=None):

    state = { \
        'position' : 0,
        'scale' : 0,
        'incl_plus' : False,
        'blanks_as_zeros' : config.PROC_BLANKS_AS_ZEROS,
        # TODO: Implement halt if no more record input
        'halt_if_no_vals' : False,
        'exception_on_fail' : True,
    }

    # pdb.set_trace()

    for ed in eds + reversion_eds:
        if isinstance(ed, tuple(FORBIDDEN_EDS)):
            raise InvalidFormat("%d edit descriptr not permitted on input")

    # Expand repeated edit decriptors
    eds = expand_edit_descriptors(eds)
    reversion_eds = expand_edit_descriptors(reversion_eds)
    # Assume one-to-one correspondance between edit descriptors and output
    # values if number of output values is not defined 
    num_out_eds = 0
    for ed in eds:
        if isinstance(ed, OUTPUT_EDS):
            num_out_eds += 1
    num_rev_out_eds = 0
    if num_vals is None:
        num_vals = num_out_eds
    for ed in reversion_eds:
        if isinstance(ed, OUTPUT_EDS):
            num_rev_out_eds += 1

    
    # Will loop forever is no output edit descriptors
    if (num_out_eds == 0):
        return []
    # Will loop forever if no output eds in reversion format and is more values
    # requested than in the format
    if (num_vals > num_out_eds) and (num_rev_out_eds == 0):
        raise ValueError('Not enough output edit descriptors in reversion format to output %d values' % num_vals)

    # May need to process multiple records, down to a higher function to supply
    # appropriate string for format
    if not hasattr(records, 'next'):
        records = iter(re.split('\r\n|\r|\n', records))
    record = _next(records, None)
    if record is None:
        return [] 
    
    # if a_widths is not None:
    #     a_widths = itertools.cycle(a_widths)

    vals = []
    finish_up = False
    ed_ind = -1
    while True:
        ed_ind += 1
        # Signal to stop when Colon edit descriptor or end of format or end of
        # reversion format reached. Also not to output any more data
        if len(vals) >= num_vals:
            finish_up = True
        # Select the appropriate edit descriptor
        if ed_ind < len(eds):
            ed = eds[ed_ind]
        else:
            rev_ed_ind = (ed_ind - len(eds)) % len(reversion_eds)
            # Reversion begun and has been instructed to halt
            if finish_up and (rev_ed_ind == 0):
                break
            ed = reversion_eds[rev_ed_ind]

        if isinstance(ed, QuotedString):
            raise InvalidFormat('Cannot have string literal in an input format')
        elif isinstance(ed, BN):
            state['blanks_as_zeros'] = False
        elif isinstance(ed, BZ):
            state['blanks_as_zeros'] = True
        elif isinstance(ed, P):
            state['scale'] = ed.scale
        elif isinstance(ed, SP):
            state['incl_plus'] = True
        elif isinstance(ed, SS):
            state['incl_plus'] = False
        elif isinstance(ed, S):
            state['incl_plus'] = config.PROC_INCL_PLUS
        elif isinstance(ed, (X, TR)):
            state['position'] = min(state['position'] + ed.num_chars, len(record))
        elif isinstance(ed, TL):
            state['position'] = max(state['position'] - ed.num_chars, 0)
        elif isinstance(ed, T):
            if (ed.num_chars - 1) < 0:
                state['position'] = 0
            elif ed.num_chars > len(record):
                state['position'] = len(record)
            else:
                state['position'] = ed.num_chars - 1
        elif isinstance(ed, Slash):
            # End of record
            record = _next(records, None)
            state['position'] = 0
            if record is None:
                break
        elif isinstance(ed, Colon):
            # Break if input value satisfied
            if finish_up:
                break
        elif isinstance(ed, (Z, O, B, I)):
            val, state = read_integer(ed, state, record)
            vals.append(val)
        elif isinstance(ed, A):
            val, state = read_string(ed, state, record)
            vals.append(val)
        elif isinstance(ed, L):
            val, state = read_logical(ed, state, record)
            vals.append(val)
        elif isinstance(ed, (F, E, D, EN, ES)):
            val, state = read_float(ed, state, record)
            vals.append(val)
        elif isinstance(ed, G):
            # Difficult to know what wanted since do not know type of input variable
            # Use the G_INPUT_TRIAL_EDS variable to try the variables
            # until one sticks
            # n.b. vals and state do not get written to if
            # exception id raised
            resolved = False
            g_trial_eds = iter(config.G_INPUT_TRIAL_EDS)
            while not resolved:
                ed_name = _next(g_trial_eds, '')
                if ed_name.upper() in ('F', 'E', 'D', 'EN', 'ES'):
                    trial_ed = F()
                    trial_ed.width = ed.width
                    trial_ed.decimal_places = ed.decimal_places
                    # pdb.set_trace()
                    try:
                        val, state = read_float(trial_ed, state.copy(), record)
                        vals.append(val)
                        resolved = True
                    except ValueError:
                        continue
                elif ed_name.upper() in ('Z', 'O', 'B', 'I'):
                    trial_ed = globals()[ed_name]()
                    trial_ed.width = ed.width
                    trial_ed.min_digits = ed.decimal_places
                    try:
                        val, state = read_integer(trial_ed, state.copy(), record)
                        vals.append(val)
                        resolved = True
                    except ValueError:
                        continue
                elif ed_name.upper() in ('L'):
                    trial_ed = L()
                    trial_ed.width = ed.width
                    try:
                        val, state = read_logical(trial_ed, state.copy(), record)
                        vals.append(val)
                        resolved = True
                    except ValueError:
                        continue
                elif ed_name.upper() in ('A'):
                    trial_ed = A()
                    trial_ed.width = ed.width
                    try:
                        val, state = read_string(trial_ed, state.copy(), record)
                        vals.append(val)
                        resolved = True
                    except ValueError:
                        continue
                elif ed_name in ('G'):
                    raise ValueError('G edit descriptor not permitted in config.G_INPUT_TRIAL_EDS')
                else:
                    raise ValueError('Unrecognised trial edit descriptor string in config.G_INPUT_TRIAL_EDS')
                    
    if config.RET_WRITTEN_VARS_ONLY:
        vals = [val for val in vals if val is not None]
    return vals[:num_vals]

def _interpret_blanks(substr, state):
    # Save leading blanks
    len_str = len(substr)
    if state['blanks_as_zeros']:
        # TODO: Are tabs blank characters?
        substr = substr.replace(' ', '0')
    else:
        substr = substr.replace(' ', '')
    # If were blanks but have been stripped away, replace with a zero
    if len(substr) == 0 and (len_str > 0):
        substr = '0'
    return substr

def _get_substr(w, record, state):
    start = max(state['position'], 0)
    end = start + w
    # if end > len(record):
    #     substr = ''
    #     # TODO: test if no chars transmitted, then poition does not change
    #     w = 0
    # else:
    substr = record[start:end]
    state['position'] = min(state['position'] + w, len(record))
    return substr, state


def _next(it, default=None):
    try:
        val = next(it)
    except StopIteration:
        val = default
    return val


def read_string(ed, state, record):
    if ed.width is None:
        # Will assume rest of record is fair game for the
        # unsized A edit descriptor
        ed.width = len(record) - state['position']
    substr, state = _get_substr(ed.width, record, state)
    val = substr.ljust(ed.width, config.PROC_PAD_CHAR)
    return (val, state)


def read_integer(ed, state, record):
    substr, state = _get_substr(ed.width, record, state)
    if ('-' in substr) and (not config.PROC_ALLOW_NEG_BOZ) and isinstance(ed, (Z, O, B)):
        if state['exception_on_fail']:
            raise ValueError('Negative numbers not permitted for binary, octal or hex')
        else:
            return (None, state)
    if isinstance(ed, Z):
        base = 16
    elif isinstance(ed, I):
        base = 10
    elif isinstance(ed, O):
        base = 8
    elif isinstance(ed, B):
        base = 2
    # If a negative is followed by blanks, Gfortran and ifort
    # interpret as a zero
    if re.match(r'^ *- +$', substr):
        substr = '0'
    # If a negative or negative and blanks, ifort interprets as
    # zero for an I edit descriptor
    if config.PROC_NEG_AS_ZERO and isinstance(ed, I) and re.match(r'^( *- *| +)$', substr):
        substr = '0'
    # If string is zero length (reading off end of record?),
    # interpret as zero so as to match what would be found in an
    # unwritten FORTRAN variable
    if substr == '':
        if config.RET_UNWRITTEN_VARS_NONE or config.RET_WRITTEN_VARS_ONLY:
            return (None, state)
        else:
            substr = '0'
    teststr = _interpret_blanks(substr, state)
    try:
        val = int(teststr, base)
    except ValueError:
        if state['exception_on_fail']:
            raise ValueError('%s is not a valid input for one of integer, octal, hex or binary' % substr)
        else:
            return (None, state)
    return (val, state)


def read_logical(ed, state, record):
    substr, state = _get_substr(ed.width, record, state)
    # Deal with case where there is no more input to read from
    if (substr == '') and (config.RET_UNWRITTEN_VARS_NONE or config.RET_WRITTEN_VARS_ONLY):
        return (None, state)
    # Remove preceding whitespace and take the first two letters as
    # uppercase for testing
    teststr = substr.upper().lstrip().lstrip('.')
    if len(teststr):
        teststr = teststr[0]
    else:
        # This is case where just a preceding period is read in
        raise ValueError('%s is not a valid boolean input' % substr)
    if teststr == 'T':
        val = True
    elif teststr == 'F':
        val = False
    else:
        if state['exception_on_fail']:
            raise ValueError('%s is not a valid boolean input' % substr)
        else:
            val = None
    return (val, state)


def read_float(ed, state, record):
    substr, state = _get_substr(ed.width, record, state)
    teststr = _interpret_blanks(substr, state)
    # When reading off end of record, get empty string,
    # interpret as 0
    if teststr == '':
        if config.RET_UNWRITTEN_VARS_NONE or config.RET_WRITTEN_VARS_ONLY:
            return (None, state)
        else:
            teststr = '0'
    # Python only understands 'E' as an exponential letter
    teststr = teststr.upper().replace('D', 'E')
    # Prepend an exponential letter if only a '-' or '+' denotes an exponent
    if 'E' not in teststr:
        teststr = teststr[0] + teststr[1:].replace('+', 'E+').replace('-', 'E-')
    # ifort allows '.' to be interpreted as 0
    if re.match(r'^ *\. *$', teststr):
        teststr = '0'
    # ifort allows '-' to be interpreted as 0
    if re.match(r'^ *- *$', teststr):
        teststr = '0'
    # ifort allows numbers to end with 'E', 'E+', 'E-' and 'D'
    # equivalents
    res = re.match(r'(.*)(E|E\+|E\-)$', teststr)
    if res:
        teststr = res.group(1)
    try:
        val = float(teststr)
    except ValueError:
        if state['exception_on_fail']:
            raise ValueError('%s is not a valid input as for an E, ES, EN or D edit descriptor' % substr)
        else:
            return (None, state)
    # Special cases: insert a decimal if none specified
    if ('.' not in teststr) and (ed.decimal_places is not None):
        val = val / 10 ** ed.decimal_places
    # Apply scale factor if exponent not supplied
    if 'E' not in teststr:
        val = val / 10 ** state['scale'] 
    return (val, state)
