import sys
import os

# Should all edit descriptor values be returned even if they were not
# written to?
RET_WRITTEN_VARS_ONLY = False
# Should 'None' values be returned when no record is available to read
# from or the FORTRAN 'default'?
RET_UNWRITTEN_VARS_NONE = True
# The order in which edit desciptors are tried by default when G edit
# descriptor encountered on input
G_INPUT_TRIAL_EDS = ['F', 'L', 'A']
# Contrary to specification, many compilers allow zero width edit
# descriptors
ALLOW_ZERO_WIDTH_EDS = True
# Set the characters that separate the records
RECORD_SEPARATOR = '\n'

# The maximum size for an integer
if sys.version_info[0] >= 3:
    PROC_MAXINT = sys.maxsize
else:
    PROC_MAXINT = sys.maxint
# Processor dependant default for including leading plus or not
PROC_INCL_PLUS = False 
# Option to allow signed binary, octal and hex on input (not a FORTRAN feature)
PROC_ALLOW_NEG_BOZ = False
# Prcessor dependant padding character
PROC_PAD_CHAR = ' '
# Interpret blanks or jsut a negative as a zero, as in ifort behaviour
PROC_NEG_AS_ZERO = True
# Show a sign for zero?
PROC_SIGN_ZERO = False
PROC_MIN_FIELD_WIDTH = 46
PROC_DECIMAL_CHAR = '.'
G0_NO_BLANKS = False
PROC_NO_LEADING_BLANK = False
# The default value if BN, BZ edit descriptors are not specified
PROC_BLANKS_AS_ZEROS = False

def reset():
    global RET_WRITTEN_VARS_ONLY, RET_UNWRITTEN_VARS_NONE, PROC_INCL_PLUS, \
        PROC_ALLOW_NEG_BOZ, PROC_PAD_CHAR, PROC_NEG_AS_ZERO, PROC_SIGN_ZERO, \
        PROC_MIN_FIELD_WIDTH, PROC_DECIMAL_CHAR, G0_NO_BLANKS, \
        PROC_NO_LEADING_BLANK, PROC_BLANKS_AS_ZEROS, PROC_MAXINT, G_INPUT_TRIAL_EDS, \
        ALLOW_ZERO_WIDTH_EDS
    G_INPUT_TRIAL_EDS = ['F', 'L', 'A']
    if sys.version_info[0] >= 3:
        PROC_MAXINT = sys.maxsize
    else:
        PROC_MAXINT = sys.maxint
    RET_WRITTEN_VARS_ONLY = False
    RET_UNWRITTEN_VARS_NONE = True
    PROC_INCL_PLUS = False
    PROC_ALLOW_NEG_BOZ = False
    PROC_PAD_CHAR = ' '
    PROC_NEG_AS_ZERO = True
    PROC_SIGN_ZERO = False
    PROC_MIN_FIELD_WIDTH = 46
    PROC_DECIMAL_CHAR = '.'
    G0_NO_BLANKS = False
    PROC_NO_LEADING_BLANK = False
    PROC_BLANKS_AS_ZEROS = False
    ALLOW_ZERO_WIDTH_EDS = True
