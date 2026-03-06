"""
This stores a list of arguments, tokenizing a string into a list of arguments.
"""
import logging
import warnings
from .exceptions import NoArgument, InputError

LOGGER = logging.getLogger(__name__)

class Argument(object):
    """ Is an argument """
    def __init__(self, string):
        """ Input string for the argument """
        self.string = string.strip()
        self.marked = False
        self.keep_quotes = False
        if len(string) == 0 or (len(string) == 1 and string in ("'", '"')):
            LOGGER.warning("Unclosed quotation detected! Ignoring lone quote")

    def lower(self):
        """ Make the argument lower-case """
        return self.string.lower()

    def isfloat(self):
        """ Determines if it can be a double """
        try:
            tmp = float(self.string)
            del tmp
            return True
        except ValueError:
            return False

    def isint(self):
        """ Determines if it can be an integer """
        try:
            tmp = int(self.string)
            del tmp
            return True
        except ValueError:
            return False

    def ismask(self):
        """ Determines if we can be an Amber mask or not """
        # The easiest way to do this is to just search for ":" or "@" in the
        # input string (or just the string *, which matches everything).
        return self.string == '*' or ':' in self.string or '@' in self.string

    # Define some casting behavior
    def __float__(self):
        self.marked = True
        return float(self.string)

    def __int__(self):
        self.marked = True
        return int(self.string)

    def __str__(self):
        self.marked = True
        if self.keep_quotes:
            return self.string

        # Strip leading and trailing quotes (only if they both exist)
        if self.string[0] == "'" and self.string[-1] == "'":
            return self.string[1:-1]
        elif self.string[0] == '"' and self.string[-1] == '"':
            return self.string[1:len(self.string)-1]

        return self.string

    def __eq__(self, other):
        """ String comparison """
        return self.string == str(other)

    def mark(self):
        """ Provide a way of forcibly marking an argument """
        self.marked = True

    def unmark(self):
        """ Provide a way of forcibly unmarking an argument """
        self.marked = False

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class ArgumentList(object):
    """ 
    List of arguments.  
    The arguments should be parsed in the following order to ensure it's done
    correctly:
        get_key_*
        get_next_int
        get_next_float
        get_next_mask
        get_next_string
    """
    def __init__(self, inp):
        """ Defines the string """
        if isinstance(inp, tuple) or isinstance(inp, list):
            _inp = [l for l in inp if l is not None and l != '']
            self.original = ' '.join([str(l) for l in _inp])
            self.tokenlist = [Argument(str(l)) for l in _inp]
        elif isinstance(inp, ArgumentList):
            self.tokenlist = []
            try:
                while True:
                    arg = Argument(inp.get_next_string())
                    self.tokenlist.append(arg)
            except NoArgument:
                pass
            self.original = ' '.join([str(a) for a in inp.tokenlist])
        else:
            # Treat as a string-like and pad with spaces
            self.original = f' {inp} '
            self.tokenize(self.original)

    def tokenize(self, instring):
        """ 
        Tokenizes the line, everything between quotes is a single token,
        including whitespace
        """
        import re
        tokenre = re.compile(r'''((?:\s+".*?"\s+?)|(?:\s+'.*?'\s+?)|(?:\S+))''')
        tokenlist = tokenre.findall(instring)
        # Make sure every non-whitespace character is consumed
        if len(tokenre.sub('', instring).strip()) > 0:
            unhandled_args = tokenre.sub('', instring).strip()
            warnings.warn(f"Not all of argument string were handled: [{unhandled_args}]")
        # Strip out whitespace in all arguments
        self.tokenlist = [Argument(token.strip()) for token in tokenlist]

    def __str__(self):
        """ Return the original input string """
        return self.original

    def get_next_int(self, optional=False, default=None):
        """ Returns the next unmarked integer in the token list """
        for token in self.tokenlist:
            if token.isint() and not token.marked:
                return int(token)
        if optional: return default
        raise NoArgument("Missing integer argument!")

    def get_next_float(self, optional=False, default=None):
        """ Returns the next unmarked float in the token list """
        for token in self.tokenlist:
            if token.isfloat() and not token.marked:
                return float(token)
        if optional: return default
        raise NoArgument("Missing floating point argument!")

    def get_next_mask(self, optional=False, default=None):
        """ Returns the next string as long as it matches a mask """
        for token in self.tokenlist:
            if token.marked: continue
            if token.ismask():
                return str(token)
        if optional: return default
        raise NoArgument("Missing Amber mask!")

    def get_next_string(self, optional=False, keep_quotes=False, default=None):
        """ Returns the next unmarked float in the token list """
        for token in self.tokenlist:
            if not token.marked:
                token.keep_quotes = keep_quotes
                return str(token)
        if optional: return default
        raise NoArgument("Missing string!")

    def unmarked(self):
        """ Returns all unmarked arguments as a list of strings """
        unmarked_tokens = []
        for token in self.tokenlist:
            if not token.marked:
                unmarked_tokens.append(token.string)
        return unmarked_tokens

    def _get_key_arg(self, key, argtype):
        """ Gets a key with a given argument type """
        for i, token in enumerate(self.tokenlist):
            # Skip all marked tokens. This allows us to allow users to specify
            # a token numerous times (useful if you want to enable an
            # arbitrary-length list of arguments
            if token.marked: continue
            if token == key:
                assert not self.tokenlist[i + 1].marked, "Marked token should not happen. BUG!"
                token.marked = True
                try:
                    return argtype(self.tokenlist[i+1])
                except ValueError as err:
                    raise InputError(str(err)) from err

        raise NoArgument('Catch me!')

    def get_key_int(self, key, default):
        """ Get an integer that follows a keyword """
        try:
            return self._get_key_arg(key, int)
        except NoArgument:
            return default

    def get_key_float(self, key, default):
        """ Get a float that follows a keyword """
        try:
            return self._get_key_arg(key, float)
        except NoArgument:
            return default

    def get_key_mask(self, key, default):
        """ Get a mask that follows a keyword """
        for i, token in enumerate(self.tokenlist):
            if token == key:
                if not self.tokenlist[i+1].ismask():
                    raise InputError(
                        f"Expected mask to follow {key}. Got {self.tokenlist[i+1].string} instead"
                    )
        try:
            return self._get_key_arg(key, str)
        except NoArgument:
            return default

    def get_key_string(self, key, default):
        """ Get a string that follows a keyword """
        try:
            return self._get_key_arg(key, str)
        except NoArgument:
            return default

    def has_key(self, key, mark=True):
        """ See if a particular keyword is present (optionally mark it) """
        for token in self.tokenlist:
            if not token.marked and token == key:
                if mark: token.mark()
                return True
        return False
