"""
Module for evaluating Amber Mask strings and translating them into lists in
which a selected atom is 1 and one that's not is 0.
"""
from ..exceptions import MaskError
from ..periodic_table import AtomicNum

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class AmberMask(object):
    """
    What is hopefully a fully-fledged Amber mask parser implemented in Python.

    Parameters
    ----------
    parm : Structure
        The topology structure for which to select atoms
    mask : str
        The mask string that selects a subset of atoms
    """

    #======================================================

    def __init__(self, parm, mask):
        self.parm = parm
        self.mask = mask.strip()

    #======================================================

    def __str__(self):
        return self.mask

    #======================================================

    def Selected(self, invert=False):
        """ Generator that returns the indexes of selected atoms

        Parameters
        ----------
        invert : bool, optional
            If True, all atoms *not* selected by the mask will be returned

        Returns
        -------
        generator of int
            Each iteration will yield the index of the next atom that has been
            selected by the mask. Atom indices are 0-based
        """
        for i, v in enumerate(self.Selection(invert=invert)):
            if v:
                yield i

    #======================================================

    def Selection(self, prnlev=0, invert=False):
        """
        Parses the mask and analyzes the result to return an atom
        selection array

        Parameters
        ----------
        prnlev : int, optional
            Print debug information on the processing of the Amber mask string.
            This is mainly useful if you are modifying the mask parser. Default
            value is 0 (no printout), values between 1 and 8 control the level
            of output (larger values produce more output). Default 0
        invert : bool, optional
            If True, the returned array will invert the selection of the mask
            (i.e., selected atoms will not be selected and vice-versa)

        Returns
        -------
        mask : list of int
            A list with length equal to the number of atoms in the assigned
            :class:`Structure <parmed.structure.Structure>` instance. Selected
            atoms will have a value of 1 placed in the corresponding slot in the
            return list while atoms not selected will be assigned 0.
        """
        from sys import stderr, stdout
        if prnlev > 2: stderr.write('In AmberMask.Selection(), debug active!\n')
        if prnlev > 5: stdout.write('original mask: ==%s==\n' % self.mask)

        # 0) See if we got the default "all" mask(*) and return accordingly
        if self.mask.strip() == '*':
            return [1 for atom in self.parm.atoms]

        # 1) preprocess input expression
        infix = self._tokenize(prnlev)
        if prnlev > 5: stdout.write('tokenized mask: ==%s==\n' % infix)

        # 2) construct postfix (RPN) notation
        postfix = self._torpn(infix, prnlev)
        if prnlev > 5: stdout.write('postfix mask: ==%s==\n' % postfix)

        # 3) evaluate the postfix notation
        if invert:
            return [1-i for i in self._evaluate(postfix, prnlev)]
        return self._evaluate(postfix, prnlev)

    #======================================================

    def _tokenize(self, prnlev):
        """ Tokenizes the mask string into individual selections:
            1. remove spaces
            2. isolate 'operands' into brackets [...]
            3. split expressions of the type :1-10@CA,CB into 2 parts;
               the 2 parts are joned with & operator and (for the sake
               of preserving precedence of other operators) enclosed by
               (...); i.e. :1-10@CA,CB is split into (:1-10 & @CA,CB)
            4. do basic error checking
        """
        buffer = '' # keeping track of a single operand
        infix = ''  # value that is returned at the end
        # flag == 0: means new operand or operand was completed & ended with ]
        # flag == 1: means operand with ":" read
        # flag == 2: means operand with "@" read
        # flag == 3: means '<' or '>' read, waiting for numbers
        flag = 0
        i = 0
        while i < len(self.mask):
            p = self.mask[i]
            # skip whitespace
            if p.isspace():
                i += 1
                continue
            # If p is an operator, is the last character, or is a ()...
            elif self._isOperator(p) or i == len(self.mask) - 1 or p in ['(',')']:
                # Deal with the last character being a wildcard that we have to
                # convert
                if p == '=' and i == len(self.mask) - 1: # wildcard
                    if flag > 0:
                        p = '*'
                    else:
                        raise MaskError("AmberMask: '=' not in name list syntax")
                # If this is the end of an operand, terminate the buffer, flush
                # it to infix, and reset flag to 0 and empty the buffer
                if flag > 0:
                    if i == len(self.mask) - 1 and p != ')': buffer += p
                    buffer += '])'
                    flag = 0
                    infix += buffer
                    buffer = ''
                if i != len(self.mask) - 1 or p == ')': infix += p
                # else if p is >,<
                if p in ['<','>']:
                    buffer = '([%s' % p
                    i += 1
                    try:
                        p = self.mask[i]
                    except IndexError:
                        raise MaskError('Bad distance syntax [%s]' % self.mask)
                    buffer += p
                    flag = 3
                    if self.parm.coordinates is None:
                        raise MaskError('<,> operators require coordinates')
                    if not p in [':','@']:
                        raise MaskError('Bad syntax [%s]' % self.mask)
            elif self._isOperand(p):
                if flag == 0:
                    buffer = '(['
                    flag = 1
                    if p != '*':
                        raise MaskError('Bad syntax [%s]' % self.mask)
                if p == '=': # wildcard
                    if flag > 0:
                        p = '*'
                    else:
                        raise MaskError("'=' not in name list syntax")
                buffer += p
            elif p == ':':
                if flag == 0:
                    buffer = '([:'
                    flag = 1
                else:
                    buffer += '])|([:'
                    flag = 1
            elif p == '@':
                if flag == 0:
                    buffer = '([@'
                    flag = 2
                elif flag == 1:
                    buffer += ']&[@'
                    flag = 2
                elif flag == 2:
                    buffer += '])|([@'
                    flag = 2
            else:
                raise MaskError('Unknown symbol (%s) expression' % p)
            i += 1
        # end while i < len(self.mask):
        # Check that each operand has at least 4 characters: [:1] and [@C], etc.
        i = 0
        n = 1 # number of characters in current operand
        flag = 0
        while i < len(infix):
            p = infix[i]
            if p == '[':
                n += 1
                flag = 1
            elif p == ']':
                if n < 4 and infix[i-1] != '*':
                    raise MaskError('empty token in infix')
                n = 1
            else:
                if flag == 1:
                    n += 1
            i += 1

        return infix + '\n' # terminating \n for next step

    #======================================================

    def _isOperator(self, char):
        """ Determines if a character is an operator """
        return len(char) == 1 and char in '!&|<>'

    #======================================================

    def _isOperand(self, char):
        """ Determines if a character is an operand """
        return len(char) == 1 and (char in "\\*/%-?,'.=+_" or char.isalnum())

    #======================================================

    def _torpn(self, infix, prnlev):
        """ Converts the infix to an RPN array """
        postfix = ''
        stack = ['\n']  # use a list as a stack. Then pop() works as expected
        flag = 0
        i = 0

        while i < len(infix):
            p = infix[i]
            if p == '[':
                postfix += p
                flag = 1
            elif p == ']':
                postfix += p
                flag = 0
            elif flag:
                postfix += p
            elif p == '(':
                stack.append(p)
            elif p == ')':
                pp = stack.pop()
                while pp != '(':
                    if pp == '\n':
                        raise MaskError('Unbalanced parentheses in Mask.')
                    postfix += pp
                    pp = stack.pop()
            # At this point both ()s are discarded
            elif p == '\n':
                pp = stack.pop()
                while pp != '\n':
                    if pp == '(':
                        raise MaskError('Unbalanced parentheses in Mask.')
                    postfix += pp
                    pp = stack.pop()
            elif self._isOperator(p):
                P1 = self._priority(p)
                P2 = self._priority(stack[len(stack)-1])
                if P1 > P2:
                    stack.append(p)
                else:
                    while P1 <= P2:
                        pp = stack.pop()
                        postfix += pp
                        P1 = self._priority(p)
                        P2 = self._priority(stack[len(stack)-1])
                    stack.append(p)
            else:
                raise MaskError('Unknown symbol %s' % p) # should not reach here
            i += 1
        # end while i < len(infix):
        return postfix

    #======================================================

    def _evaluate(self, postfix, prnlev):
        """ Evaluates a postfix in RPN format and returns a selection array """
        from sys import stderr
        buffer = ''
        stack = []

        pos = 0 # position in postfix
        while pos < len(postfix):
            p = postfix[pos]
            if p == '[': buffer = ''
            elif p == ']': # end of the token
                ptoken = buffer
                pmask = self._selectElemMask(ptoken)
                stack.append(pmask)
            elif self._isOperand(p) or p in [':','@']:
                buffer += p
            elif p in ['&','|']:
                pmask1 = None
                pmask2 = None
                try:
                    pmask1 = stack.pop()
                    pmask2 = stack.pop()
                    pmask = self._binop(p, pmask1, pmask2)
                except IndexError:
                    raise MaskError('Illegal binary operation')
                stack.append(pmask)
            elif p in ['<','>']:
                if pos < len(postfix)-1 and postfix[pos+1] in [':','@']:
                    buffer += p
                else:
                    try:
                        pmask1 = stack.pop() # distance criteria
                        pmask2 = stack.pop()
                        pmask = self._selectDistd(pmask1, pmask2)
                    except IndexError:
                        return [0 for a in self.parm.atoms]
                    stack.append(pmask)
            elif p == '!':
                try:
                    pmask1 = stack.pop()
                except IndexError:
                    raise MaskError('Illegal ! operation')
                pmask = self._neg(pmask1)
                stack.append(pmask)
            else:
                raise MaskError('Unknown symbol evaluating RPN: %s' % p)
            pos += 1
        # end while i < len(postfix)

        try:
            pmask = stack.pop()
        except IndexError:
            raise MaskError('Empty stack -- no available operands')

        if stack:
            raise MaskError('There may be missing operands in the mask')

        if prnlev > 7:
            stderr.write('%d atoms selected by %s' % (sum(pmask), self.mask))

        return pmask

    #======================================================

    def _neg(self, pmask1):
        """ Negates a given mask """
        return pmask1.Not()

    #======================================================

    def _selectDistd(self, pmask1, pmask2):
        """ Selects atoms based on a distance criteria """
        # pmask1 is either @<number> or :<number>, and represents the distance
        # criteria. pmask2 is the selection of atoms from which the distance is
        # evaluated.
        pmask = _mask(len(self.parm.atoms))
        # Determine if we want > or <
        if pmask1[0] == '<':
            cmp = lambda x, y: x < y
        elif pmask1[0] == '>':
            cmp = lambda x, y: x > y
        else: # Should never execute this
            raise MaskError('Unknown comparison criteria for distance mask: %s' % pmask1[0])
        pmask1 = pmask1[1:]
        if pmask1[0] not in ':@': # Should never execute this
            raise MaskError('Bad distance criteria for mask: %s' % pmask1)
        try:
            distance = float(pmask1[1:])
        except (TypeError, ValueError):
            raise MaskError('Distance must be a number: %s' % pmask1[1:])
        distance *= distance # Faster to compare square of distance
        # First select all atoms that satisfy the distance. If we ended up
        # choosing residues, then we will go back through afterwards and select
        # entire residues when one of the atoms in that residue is selected.
        idxlist = [i for i, val in enumerate(pmask2) if val == 1]
        for i, atomi in enumerate(self.parm.atoms):
            for j in idxlist:
                atomj = self.parm.atoms[j]
                dx = atomi.xx - atomj.xx
                dy = atomi.xy - atomj.xy
                dz = atomi.xz - atomj.xz
                d2 = dx*dx + dy*dy + dz*dz
                if cmp(d2, distance):
                    pmask[i] = 1
                    break
        # Now see if we have to select all atoms in residues with any selected
        # atoms
        if pmask1[0] == ':':
            for res in self.parm.residues:
                for atom in res.atoms:
                    if pmask[atom.idx] == 1:
                        for atom in res.atoms:
                            pmask[atom.idx] = 1
                        break
        return pmask

    #======================================================

    def _selectElemMask(self, ptoken):
        """ Selects an element mask """
        # some constants
        ALL = 0
        NUMLIST = 1
        NAMELIST = 2
        TYPELIST = 3
        ELEMLIST = 4
        # define the mask object and empty buffer
        pmask = _mask(len(self.parm.atoms))
        buffer = ''
        buffer_p = 0
        # This is a residue NUMber LIST
        if ptoken.startswith(':'):
            reslist = NUMLIST
            pos = 1
            while pos < len(ptoken):
                p = ptoken[pos]
                buffer += p
                buffer_p += 1
                if p == '*' and ptoken[pos-1] != '\\':
                    if buffer_p == 1 and (pos == len(ptoken) - 1 or ptoken[pos+1] == ','):
                        reslist = ALL
                    elif reslist == NUMLIST:
                        reslist = NAMELIST
                elif p.isalpha() or p in '_?*':
                    reslist = NAMELIST
                if pos == len(ptoken) - 1:
                    buffer_p = 0
                if len(buffer) != 0 and buffer_p == 0:
                    if reslist == ALL:
                        pmask.select_all()
                    elif reslist == NUMLIST:
                        self._residue_numlist(buffer, pmask)
                    elif reslist == NAMELIST:
                        self._residue_namelist(buffer, pmask)
                    reslist = NUMLIST
                pos += 1
        elif ptoken.startswith('@'):
            atomlist = NUMLIST
            pos = 1
            while pos < len(ptoken):
                p = ptoken[pos]
                buffer += p
                buffer_p += 1
                if p == '*' and ptoken[pos-1] != "\\":
                    if atomlist == NUMLIST:
                        atomlist = NAMELIST
                elif p.isalpha() or p in '?*_':
                    if atomlist == NUMLIST:
                        atomlist = NAMELIST
                elif p == '%':
                    atomlist = TYPELIST
                elif p == '/':
                    atomlist = ELEMLIST
                if pos == len(ptoken) - 1:
                    buffer_p = 0

                if len(buffer) != 0 and buffer_p == 0:
                    if atomlist == ALL:
                        pmask.select_all()
                    elif atomlist == NUMLIST:
                        self._atom_numlist(buffer, pmask)
                    elif atomlist == NAMELIST:
                        self._atom_namelist(buffer, pmask)
                    elif atomlist == TYPELIST:
                        self._atom_typelist(buffer[1:], pmask)
                    elif atomlist == ELEMLIST:
                        self._atom_elemlist(buffer[1:], pmask)
                pos += 1
        elif ptoken.strip() == '*':
            pmask.select_all()
        elif ptoken[0] in ['<','>']:
            return ptoken
        else: # Should never reach here
            raise MaskError('Mask is missing : and @')
        # end if ':' in ptoken:

        return pmask

    #======================================================

    def _atom_numlist(self, instring, mask):
        """ Fills a _mask based on atom numbers """
        buffer = ''
        pos = 0
        at1 = at2 = dash = 0
        while pos < len(instring):
            p = instring[pos]
            if p.isdigit():
                buffer += p
            if p == ',' or pos == len(instring) - 1:
                if dash == 0:
                    at1 = int(buffer)
                    self._atnum_select(at1, at1, mask)
                else:
                    at2 = int(buffer)
                    self._atnum_select(at1, at2, mask)
                    dash = 0
                buffer = ''
            elif p == '-':
                at1 = int(buffer)
                dash = 1
                buffer = ''
            if not (p.isdigit() or p in [',','-']):
                raise MaskError('Unknown symbol in atom number parsing [%s]'%p)
            pos += 1

    #======================================================

    def _atom_namelist(self, instring, mask, key='name'):
        """ Fills a _mask based on atom names/types """
        buffer = ''
        pos = 0
        while pos < len(instring):
            p = instring[pos]
            if p.isalnum() or p in "\\*?+'-_":
                buffer += p
            if p == ',' or pos == len(instring) - 1:
                if '-' in buffer and buffer[0].isdigit():
                    self._atom_numlist(buffer, mask)
                else:
                    self._atname_select(buffer, mask, key)
                buffer = ''
            if not (p.isalnum() or p in "\\,?*'+-_"):
                raise MaskError('Unrecognized symbol in atom name parsing [%s]' % p)
            pos += 1

    #======================================================

    def _atom_typelist(self, buffer, mask):
        """ Fills a _mask based on atom types """
        self._atom_namelist(buffer, mask, key='type')

    #======================================================

    def _atom_elemlist(self, buffer, mask):
        """
        Fills a _mask based on atom elements. For now it will just be Atom
        names, since elements are not stored in the prmtop anywhere.
        """
        self._atom_namelist(buffer, mask, key='element')

    #======================================================

    def _residue_numlist(self, instring, mask):
        """ Fills a _mask based on residue numbers """
        buffer = ''
        pos = 0
        at1 = at2 = dash = 0
        while pos < len(instring):
            p = instring[pos]
            if p.isdigit():
                buffer += p
            if p == ',' or pos == len(instring) - 1:
                if dash == 0:
                    at1 = int(buffer)
                    self._resnum_select(at1, at1, mask)
                else:
                    try:
                        at2 = int(buffer)
                    except ValueError:
                        raise MaskError('Bad mask: error in integer conversion')
                    self._resnum_select(at1, at2, mask)
                    dash = 0
                buffer = ''
            elif p == '-':
                at1 = int(buffer)
                dash = 1
                buffer = ''
            pos += 1

    #======================================================

    def _residue_namelist(self, instring, mask):
        """ Fills a _mask based on residue names """
        buffer = ''
        pos = 0
        while pos < len(instring):
            p = instring[pos]
            if p.isalnum() or p in ['*','?','+',"'",'-']:
                buffer += p
            if p == ',' or pos == len(instring) - 1:
                if '-' in buffer and buffer[0].isdigit():
                    self._residue_numlist(buffer, mask)
                else:
                    self._resname_select(buffer, mask)
                buffer = ''
            if not (p.isalnum() or p in ",?*'+-"):
                raise MaskError('Unknown symbol in residue name parsing [%s]' % p)
            pos += 1

    #======================================================

    def _atnum_select(self, at1, at2, mask):
        """ Fills a _mask array between atom numbers at1 and at2 """
        for i in range(at1-1, at2): mask[i] = 1

    #======================================================

    def _resnum_select(self, res1, res2, mask):
        """ Fills a _mask array between residues res1 and res2 """
        for i, atom in enumerate(self.parm.atoms):
            res = atom.residue.idx + 1
            if res >= res1 and res <= res2: mask[i] = 1

    #======================================================

    def _atname_select(self, atname, mask, key='name'):
        """ Fills a _mask array with all atom names of a given name """
        if atname.isdigit():
            atname = int(atname) - 1
            for i, atom in enumerate(self.parm.atoms):
                mask[i] = mask[i] | int(atname == i)
        elif key == 'element':
            try:
                for i, atom in enumerate(self.parm.atoms):
                    mask[i] = mask[i] | int(AtomicNum[atname] == atom.atomic_number)
            except KeyError:
                raise MaskError('Unknown element %s' % atname)
        else:
            for i, atom in enumerate(self.parm.atoms):
                mask[i] = mask[i] | int(_nameMatch(atname, getattr(atom, key)))

    #======================================================

    def _resname_select(self, resname, mask):
        """ Fills a _mask array with all residue names of a given name """
        for i, atm in enumerate(self.parm.atoms):
            if _nameMatch(resname, atm.residue.name):
                mask[i] = 1
            elif resname.isdigit():
                mask[i] = mask[i] | int(int(resname) == atm.residue.idx + 1)

    #======================================================

    def _binop(self, op, pmask1, pmask2):
        """ Does a binary operation on a pair of masks """
        if op == '&':
            return pmask1.And(pmask2)
        if op == '|':
            return pmask1.Or(pmask2)
        raise MaskError('Unknown operator [%s]' % op)

    #======================================================

    def _priority(self, op):
        if op in ['>','<']: return 6
        if op in ['!']: return 5
        if op in ['&']: return 4
        if op in ['|']: return 3
        if op in ['(']: return 2
        if op in ['\n']: return 1

        raise MaskError('Unknown operator [%s] in Mask ==%s==' % (op, self.mask))

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def _nameMatch(atnam1, atnam2):
    """
    Determines if atnam1 matches atnam2, where atnam1 can have * as a wildcard
    and spaces are ignored. atnam2 should come from the prmtop. We'll use regex
    to do this.

    We will replace * with a regex that will match any alphanumeric character
    0 or more times: * --> \\w*

    We will replace ? with a regex that will match exactly 1 alphanumeric
    character: ? --> \\w

    Then, we will substitute all instances of atnam2 in atnam1 with ''. If it's
    a complete match, then our result will be a blank string (and will evaluate
    to False for boolean conditions). If it's not blank, then it's not a full
    match and should return False
    """
    import re
    atnam1 = str(atnam1).replace(' ','')
    atnam2 = str(atnam2).replace(' ','')
    # Replace amber mask wildcards with appropriate regex wildcards and protect
    # the + (but protect backslashes)
    R = '<!PROTECT!>'
    atnam1 = atnam1.replace('\\*', R).replace('*',r'\S*').replace(R, '*')
    atnam1 = atnam1.replace('\\?', R).replace('?',r'\S').replace(R, '?')
    atnam1 = atnam1.replace('\\+', R).replace('+',r'\+').replace(R, '+')
    # Now replace just the first instance of atnam2 in atnam2 with '', and
    # return *not* that
    return atnam1 == atnam2 or not bool(re.sub(atnam1, '', atnam2, 1))

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class _mask(list):
    """ Mask array; only used by AmberMask """

    def __init__(self, natom):
        self.natom = natom
        list.__init__(self, [0 for i in range(natom)])

    def append(self, *args, **kwargs):
        raise MaskError('_mask is a fixed-length array!')

    def extend(self, *args, **kwargs):
        raise MaskError('_mask is a fixed-length array!')

    def pop(self, *args, **kwargs):
        return self[-1]

    def remove(self, *args, **kwargs):
        raise MaskError('_mask is a fixed-length array!')

    def And(self, other):
        if self.natom != other.natom:
            raise MaskError("_mask: and() requires another mask of equal size!")
        new_mask = _mask(self.natom)
        for i in range(len(self)):
            new_mask[i] = int(self[i] and other[i])
        return new_mask

    def Or(self, other):
        if self.natom != other.natom:
            raise MaskError('_mask: or() requires another mask of equal size!')
        new_mask = _mask(self.natom)
        for i in range(len(self)):
            new_mask[i] = int(self[i] or other[i])
        return new_mask

    def Not(self):
        new_mask = _mask(self.natom)
        for i in range(self.natom):
            new_mask[i] = 1 - self[i]
        return new_mask

    def select_all(self):
        for i in range(self.natom):
            self[i] = 1

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
