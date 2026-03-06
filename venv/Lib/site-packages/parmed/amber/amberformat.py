"""
This is a generalization of the readparm.AmberParm class to handle similar
Amber-style files with %FLAG/%FORMAT tags
"""
import datetime
import re
from contextlib import closing
from copy import copy
from math import ceil

from ..constants import PrmtopPointers, AMBER_ELECTROSTATIC, CHARMM_ELECTROSTATIC
from ..exceptions import AmberError
from ..formats.registry import FileFormatType
from ..utils.io import genopen
from ..utils.fortranformat import FortranRecordReader, FortranRecordWriter

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class FortranFormat:
    """
    Processes Fortran format strings according to the Fortran specification for
    such formats. This object handles reading and writing data with any valid
    Fortran format. It does this by using the `fortranformat` project
    [https://bitbucket.org/brendanarnold/py-fortranformat].

    However, while `fortranformat` is very general and adheres well to the
    standard, it is very slow. As a result, simple, common format strings have
    been optimized and processes reads and writes between 3 and 5 times faster.
    The format strings (case-insensitive) of the following form (where # can be
    replaced by any number) are optimized:
        - #E#.#
        - #D#.#
        - #F#.#
        - #(F#.#)
        - #a#
        - #I#

    Parameters
    ----------
    format_string : str
        The Fortran Format string to process
    strip_strings : bool=True
        If True, strings are stripped before being processed by stripping
        (only) trailing whitespace
    """

    strre = re.compile(r'(\d+)?a(\d+)$', re.I)
    intre = re.compile(r'(\d+)?i(\d+)$', re.I)
    floatre = re.compile(r'(\d+)?[edf](\d+)\.(\d+)$', re.I)
    floatre2 = re.compile(r'(\d+)?\([edf](\d+)\.(\d+)\)$', re.I)

    #===================================================

    def __init__(self, format_string, strip_strings=True):
        """
        Sets the format string and determines how we will read and write
        strings using this format
        """
        self.format = format_string
        self.strip_strings = strip_strings # for ease of copying

        # Define a function that processes all arguments prior to adding them to
        # the returned list. By default, do nothing, but this allows us to
        # optionally strip whitespace from strings.
        self.process_method = lambda x: x

        if FortranFormat.strre.match(format_string):
            rematch = FortranFormat.strre.match(format_string)
            # replace our write() method with write_string to force left-justify
            self.type, self.write = str, self._write_string
            nitems, itemlen = rematch.groups()
            if nitems is None:
                self.nitems = 1
            else:
                self.nitems = int(nitems)
            self.itemlen = int(itemlen)
            self.fmt = '%s'
            # See if we want to strip the strings
            if strip_strings: self.process_method = lambda x: x.strip()

        elif FortranFormat.intre.match(format_string):
            self.type = int
            rematch = FortranFormat.intre.match(format_string)
            nitems, itemlen = rematch.groups()
            if nitems is None:
                self.nitems = 1
            else:
                self.nitems = int(nitems)
            self.itemlen = int(itemlen)
            self.fmt = '%%%dd' % self.itemlen

        elif FortranFormat.floatre.match(format_string):
            self.type = float
            rematch = FortranFormat.floatre.match(format_string)
            nitems, itemlen, num_decimals = rematch.groups()
            if nitems is None:
                self.nitems = 1
            else:
                self.nitems = int(nitems)
            self.itemlen = int(itemlen)
            self.num_decimals = int(num_decimals)
            if 'F' in format_string.upper():
                self.fmt = '%%%s.%sF' % (self.itemlen, self.num_decimals)
            else:
                self.fmt = '%%%s.%sE' % (self.itemlen, self.num_decimals)

        elif FortranFormat.floatre2.match(format_string):
            self.type = float
            rematch = FortranFormat.floatre2.match(format_string)
            nitems, itemlen, num_decimals = rematch.groups()
            if nitems is None:
                self.nitems = 1
            else:
                self.nitems = int(nitems)
            self.itemlen = int(itemlen)
            self.num_decimals = int(num_decimals)
            if 'F' in format_string.upper():
                self.fmt = '%%%s.%sF' % (self.itemlen, self.num_decimals)
            else:
                self.fmt = '%%%s.%sE' % (self.itemlen, self.num_decimals)

        else:
            # We tried... now just use the fortranformat package
            self._reader = FortranRecordReader(format_string)
            self._writer = FortranRecordWriter(format_string)
            self.write = self._write_ffwriter
            self.read = self._read_ffreader

    #===================================================

    def __copy__(self):
        return type(self)(self.format, self.strip_strings)

    #===================================================

    def __str__(self):
        return self.format

    def __repr__(self):
        return "<%s: %s>" % (type(self).__name__, self.format)

    #===================================================

    def write(self, items, dest):
        """
        Writes an iterable of data (or a single item) to the passed file-like
        object

        Parameters
        ----------
        items : iterable or single float/str/int
            These are the objects to write in this format. The types of each
            item should match the type specified in this Format for that
            argument
        dest : file or file-like
            This is the file to write the data to. It must have a `write` method
            or an AttributeError will be raised

        Notes
        -----
        This method may be replaced with _write_string (for #a#-style formats)
        or _write_ffwriter in the class initializer if no optimization is
        provided for this format, but the call signatures and behavior are the
        same for each of those functions.
        """
        if hasattr(items, '__iter__') and not isinstance(items, str):
            mod = self.nitems - 1
            for i, item in enumerate(items):
                dest.write(self.fmt % item)
                if i % self.nitems == mod:
                    dest.write('\n')
            if i % self.nitems != mod:
                dest.write('\n')
        else:
            dest.write(self.fmt % items)
            dest.write('\n')

    #===================================================

    def _write_string(self, items, dest):
        """ Writes a list/tuple of strings """
        if hasattr(items, '__iter__') and not isinstance(items, str):
            mod = self.nitems - 1
            for i, item in enumerate(items):
                dest.write((self.fmt % item).ljust(self.itemlen))
                if i % self.nitems == mod:
                    dest.write('\n')
            if i % self.nitems != mod:
                dest.write('\n')
        else:
            dest.write((self.fmt % items).ljust(self.itemlen))
            dest.write('\n')

    #===================================================

    def _read_nostrip(self, line):
        """
        Reads the line and returns converted data. Special-cased for flags that
        may contain 'blank' data. ugh.
        """
        line = line.rstrip('\n')
        nitems = int(ceil(len(line) / self.itemlen))
        ret = [0 for i in range(nitems)]
        start, end = 0, self.itemlen
        for i in range(nitems):
            ret[i] = self.process_method(self.type(line[start:end]))
            start = end
            end += self.itemlen
        return ret

    #===================================================

    def read(self, line):
        """ Reads the line and returns the converted data """
        line = line.rstrip()
        nitems = int(ceil(len(line) / self.itemlen))
        ret = [0 for i in range(nitems)]
        start, end = 0, self.itemlen
        for i in range(nitems):
            ret[i] = self.process_method(self.type(line[start:end]))
            start = end
            end += self.itemlen
        return ret

    #===================================================

    def _read_ffreader(self, line):
        """ Reads the line and returns the converted data """
        return self._reader.read(line.rstrip())

    #===================================================

    def _write_ffwriter(self, items, dest):
        dest.write('%s\n' % self._writer.write(items))

    #===================================================

    def __eq__(self, other):
        return (self.format == other.format and
                self.strip_strings == other.strip_strings)

    #===================================================

    def __hash__(self):
        return hash((self.format, self.strip_strings))

    #===================================================

    def __getstate__(self):
        return dict(format=self.format, strip_strings=self.strip_strings)

    def __setstate__(self, d):
        self.__init__(d['format'], d['strip_strings'])

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AmberFormat(metaclass=FileFormatType):
    """
    A class that can parse and print files stored in the Amber topology or MDL
    format. In particular, these files have the general form:

    ```
    %VERSION VERSION_STAMP = V00001.000  DATE = XX/XX/XX  XX:XX:XX
    %FLAG <FLAG_NAME>
    %COMMENT <comments>
    %FORMAT(<Fortran_Format>)
    ... data corresponding to that Fortran Format
    %FLAG <FLAG_NAME2>
    %COMMENT <comments>
    %FORMAT(<Fortran_Format>)
    ... data corresponding to that Fortran Format
    ```

    where the `%COMMENT` sections are entirely optional

    Parameters
    ----------
    fname : str=None
        If provided, this file is parsed and the data structures will be loaded
        from the data in this file

    Attributes
    ----------
    parm_data : dict {str : list}
        A dictionary that maps FLAG names to all of the data contained in that
        section of the Amber file.
    formats : dict {str : FortranFormat}
        A dictionary that maps FLAG names to the FortranFormat instance in which
        the data is stored in that section
    parm_comments : dict {str : list}
        A dictionary that maps FLAG names to the list of COMMENT lines that were
        stored in the original file
    flag_list : list
        An ordered list of all FLAG names. This must be kept synchronized with
        `parm_data`, `formats`, and `parm_comments` such that every item in
        `flag_list` is a key to those 3 dicts and no other keys exist
    charge_flag : str='CHARGE'
        The name of the name of the FLAG that describes partial atomic charge
        data. If this flag is found, then its data are multiplied by the
        ELECTROSTATIC_CONSTANT to convert back to fractions of electrons
    version : str
        The VERSION string from the Amber file
    name : str
        The file name of the originally parsed file (set to the fname parameter)
    """
    #===================================================

    @staticmethod
    def id_format(filename):
        """
        Identifies the file type as either Amber-format file (like prmtop) or an
        old-style topology file.

        Parameters
        ----------
        filename : str
            Name of the file to check format for

        Returns
        -------
        is_fmt : bool
            True if it is an Amber-style format, False otherwise
        """
        if isinstance(filename, str):
            with closing(genopen(filename, 'r')) as f:
                lines = [f.readline() for i in range(5)]
        elif all(hasattr(filename, attr) for attr in ['readline', 'seek', 'tell']):
            cur = filename.tell()
            lines = [filename.readline() for i in range(5)]
            filename.seek(cur)

        if lines[0].startswith('%VERSION'):
            return True
        # Try old-style format
        try:
            return AmberFormat().rdparm_old(lines, check=True)
        except ValueError:
            return False

    #===================================================

    @staticmethod
    def parse(filename, *args, **kwargs):
        """
        Meant for use with the automatic file loader, this will automatically
        return a subclass of AmberFormat corresponding to what the information
        in the prmtop file contains (i.e., either an AmberParm, ChamberParm,
        AmoebaParm, or AmberFormat)
        """
        from .readparm import LoadParm
        from ._tinkerparm import BeemanRestart
        try:
            return LoadParm(filename, *args, **kwargs)
        except (IndexError, KeyError):
            parm = AmberFormat(filename, *args, **kwargs)
            if 'ATOMIC_COORDS_LIST' in parm.parm_data:
                return BeemanRestart.from_rawdata(parm)
            return parm

    #===================================================

    def __init__(self, fname=None):
        """ Constructor.  Read a file if given """
        self._ncopies = 0
        self.parm_data = {}
        self.formats = {}
        self.parm_comments = {}
        self.flag_list = []
        self.version = None
        self.charge_flag = 'CHARGE'
        self.name = fname

        if fname is not None:
            self.rdparm(fname)

    #===================================================

    def __copy__(self):
        """ Copy all of the data """
        self._ncopies += 1
        other = type(self)()
        other.flag_list = self.flag_list[:]
        other.version = self.version
        if self.name is not None:
            other.name = self.name + '_copy%d' % self._ncopies
        else:
            other.name = None
        other.charge_flag = self.charge_flag
        other.parm_data = {}
        other.parm_comments = {}
        other.formats = {}
        for flag in other.flag_list:
            other.parm_data[flag] = self.parm_data[flag][:]
            other.parm_comments[flag] = self.parm_comments[flag][:]
            other.formats[flag] = copy(self.formats[flag])
        return other

    #===================================================

    def view_as(self, cls):
        """
        Returns a view of the current object as another object.

        Parameters
        ----------
        cls : type
            Class definition of an AmberParm subclass for the current object to
            be converted into

        Returns
        -------
        instance of cls initialized from data in this object. This is NOT a deep
        copy, so modifying the original object may modify this. The copy
        function will create a deep copy of any AmberFormat-derived object
        """
        # If these are the same classes, just return the original instance,
        # since there's nothing to do. Classes are singletons, so use "is"
        if type(self) is cls:
            return self
        return cls.from_rawdata(self)

    #===================================================

    def rdparm(self, fname, slow=False):
        """ Parses the Amber format file """
        self.name = fname
        self.version = None # reset all top info each time rdparm is called
        self.formats = {}
        self.parm_data = {}
        self.parm_comments = {}
        self.flag_list = []

        # See if we have the optimized parser available
        try:
            from . import _rdparm
        except ImportError:
            return self.rdparm_slow(fname)

        # The optimized parser only works on local, uncompressed files
        # TODO: Add gzip and bzip2 support to the optimized reader
        if (hasattr(fname, 'read') or slow
            or fname.startswith('http://') or fname.startswith('https://')
            or fname.startswith('ftp://')
            or fname.endswith('.bz2') or fname.endswith('.gz')):

            return self.rdparm_slow(fname)

        # We have the optimized version and a local file
        try:
            ret = _rdparm.rdparm(fname)
        except TypeError:
            # This is raised if VERSION is not found
            with closing(genopen(fname, 'r')) as f:
                return self.rdparm_old(f.readlines())
        else:
            # Unpack returned contents
            parm_data, parm_comments, formats, unkflg, flag_list, version = ret
            # Now assign them to instance attributes and process where necessary
            self.parm_data = parm_data
            self.parm_comments = parm_comments
            for key in formats:
                self.formats[key] = FortranFormat(formats[key])
            self.flag_list = flag_list
            self.version = version
            # Now we have to process all of those sections that the optimized
            # parser couldn't figure out
            for flag in unkflg:
                rawdata = self.parm_data[flag]
                self.parm_data[flag] = []
                for line in rawdata:
                    self.parm_data[flag].extend(self.formats[flag].read(line))
            if 'CTITLE' in self.parm_data:
                CHARGE_SCALE = CHARMM_ELECTROSTATIC
            else:
                CHARGE_SCALE = AMBER_ELECTROSTATIC
            try:
                for i, chg in enumerate(self.parm_data[self.charge_flag]):
                    self.parm_data[self.charge_flag][i] = chg / CHARGE_SCALE
            except KeyError:
                pass

    #===================================================

    def rdparm_slow(self, fname):
        """
        Parses the Amber format file. This parser is written in pure Python and
        is therefore slower than the C++-optimized version
        """

        current_flag = ''
        fmtre = re.compile(r'%FORMAT *\((.+)\)')
        version = None

        if isinstance(fname, str):
            prm = genopen(fname, 'r')
            own_handle = True
        elif hasattr(fname, 'read'):
            prm = fname
            own_handle = False
        else:
            raise TypeError('%s must be a file name or file-like object' % fname)

        # Open up the file and read the data into memory
        for line in prm:
            if line[0] == '%':
                if line[0:8] == '%VERSION':
                    self.version = line.strip()
                    continue
                elif line[0:5] == '%FLAG':
                    current_flag = line[6:].strip()
                    self.formats[current_flag] = ''
                    self.parm_data[current_flag] = []
                    self.parm_comments[current_flag] = []
                    self.flag_list.append(current_flag)
                    continue
                elif line[0:8] == '%COMMENT':
                    self.parm_comments[current_flag].append(line[9:].strip())
                    continue
                elif line[0:7] == '%FORMAT':
                    fmt = FortranFormat(fmtre.match(line).groups()[0])
                    # RESIDUE_ICODE can have a lot of blank data...
                    if current_flag == 'RESIDUE_ICODE':
                        fmt.read = fmt._read_nostrip
                    self.formats[current_flag] = fmt
                    continue
            try:
                self.parm_data[current_flag].extend(fmt.read(line))
            except KeyError:
                if version is not None:
                    raise
                break # Skip out of the loop down to the old-format parser

        # convert charges to fraction-electrons
        if 'CTITLE' in self.parm_data:
            CHARGE_SCALE = CHARMM_ELECTROSTATIC
        else:
            CHARGE_SCALE = AMBER_ELECTROSTATIC
        if self.charge_flag in self.parm_data:
            for i, chg in enumerate(self.parm_data[self.charge_flag]):
                self.parm_data[self.charge_flag][i] = chg / CHARGE_SCALE
        # If we don't have a version, then read in an old-file topology
        if self.version is None:
            prm.seek(0)
            return self.rdparm_old(prm.readlines())
        if own_handle:
            prm.close()
        return

    #===================================================

    def rdparm_old(self, prmtop_lines, check=False):
        """
        This reads an old-style topology file and stores the results in the
        same data structures as a new-style topology file

        Parameters
        ----------
        prmtop_lines : list of str
            List of all lines in the prmtop file
        check : bool, optional
            If True, only the first couple sections will be read to determine if
            this is, in fact, an old-style topology file
        """
        def read_integer(line_idx, lines, num_items):
            # line_idx should be the line _before_ the first line you
            # want data from.
            i, tmp_data = 0, []
            while i < num_items:
                idx = i % 12
                if idx == 0:
                    line_idx += 1
                try:
                    tmp_data.append(int(lines[line_idx][idx*6:idx*6+6]))
                except ValueError:
                    raise ValueError(
                            'Error parsing line %d, token %d [%s]: Problem '
                            'during integer read.' % (line_idx, idx,
                            lines[line_idx][idx*6:idx*6+6])
                    )
                i += 1
            # If we had no items, we need to jump a line:
            if num_items == 0: line_idx += 1
            return tmp_data, line_idx

        def read_string(line_idx, lines, num_items):
            # line_idx should be the line _before_ the first line you
            # want data from.
            i, tmp_data = 0, []
            while i < num_items:
                idx = i % 20
                if idx == 0:
                    line_idx += 1
                tmp_data.append(lines[line_idx][idx*4:idx*4+4])
                i += 1
            # If we had no items, we need to jump a line:
            if num_items == 0: line_idx += 1
            return tmp_data, line_idx

        def read_float(line_idx, lines, num_items):
            # line_idx should be the line _before_ the first line you
            # want data from.
            i, tmp_data = 0, []
            while i < num_items:
                idx = i % 5
                if idx == 0:
                    line_idx += 1
                try:
                    tmp_data.append(float(lines[line_idx][idx*16:idx*16+16]))
                except ValueError:
                    raise ValueError(
                            'Error parsing line %d, token %d [%s]: Problem '
                            'during floating point read.' % (line_idx, idx,
                            lines[line_idx][idx*16:idx*16+16])
                    )
                i += 1
            # If we had no items, we need to jump a line:
            if num_items == 0: line_idx += 1
            return tmp_data, line_idx

        # First add a title
        self.add_flag('TITLE', '20a4', data=['| Converted old-style topology'])

        # Next, read in the pointers
        line_idx = 0
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, 30)
        # Add a final pointer of 0, which corresponds to NUMEXTRA
        tmp_data.append(0)
        self.add_flag('POINTERS', '10I8', data=tmp_data)

        # Set some of the pointers we need
        natom = self.parm_data['POINTERS'][PrmtopPointers.NATOM]
        ntypes = self.parm_data['POINTERS'][PrmtopPointers.NTYPES]
        nres = self.parm_data['POINTERS'][PrmtopPointers.NRES]
        numbnd = self.parm_data['POINTERS'][PrmtopPointers.NUMBND]
        numang = self.parm_data['POINTERS'][PrmtopPointers.NUMANG]
        nptra = self.parm_data['POINTERS'][PrmtopPointers.NPTRA]
        natyp = self.parm_data['POINTERS'][PrmtopPointers.NATYP]
        nbonh = self.parm_data['POINTERS'][PrmtopPointers.NBONH]
        nbona = self.parm_data['POINTERS'][PrmtopPointers.NBONA]
        ntheth = self.parm_data['POINTERS'][PrmtopPointers.NTHETH]
        ntheta = self.parm_data['POINTERS'][PrmtopPointers.NTHETA]
        nex = self.parm_data['POINTERS'][PrmtopPointers.NEXT]
        nphia = self.parm_data['POINTERS'][PrmtopPointers.NPHIA]
        nphb = self.parm_data['POINTERS'][PrmtopPointers.NPHB]
        nphih = self.parm_data['POINTERS'][PrmtopPointers.NPHIH]

        # This is enough to convince me that we have an old-style prmtop if we
        # have the number of integers I suspect we should
        if check:
            return len(tmp_data) == 31

        # Next read in the atom names
        tmp_data, line_idx = read_string(line_idx, prmtop_lines, natom)
        self.add_flag('ATOM_NAME', '20a4', data=tmp_data)

        # Next read the charges
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, natom)
        # Divide by the electrostatic constant
        tmp_data = [x / AMBER_ELECTROSTATIC for x in tmp_data]
        self.add_flag('CHARGE', '5E16.8', data=tmp_data)

        # Next read the masses
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, natom)
        self.add_flag('MASS', '5E16.8', data=tmp_data)

        # Next read atom type index
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, natom)
        self.add_flag('ATOM_TYPE_INDEX', '10I8', data=tmp_data)

        # Next read number excluded atoms
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, natom)
        self.add_flag('NUMBER_EXCLUDED_ATOMS', '10I8', data=tmp_data)

        # Next read nonbonded parm index
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, ntypes**2)
        self.add_flag('NONBONDED_PARM_INDEX', '10I8', data=tmp_data)

        # Next read residue label
        tmp_data, line_idx = read_string(line_idx, prmtop_lines, nres)
        self.add_flag('RESIDUE_LABEL', '20a4', data=tmp_data)

        # Next read residue pointer
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, nres)
        self.add_flag('RESIDUE_POINTER', '10I8', data=tmp_data)

        # Next read bond force constant
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, numbnd)
        self.add_flag('BOND_FORCE_CONSTANT', '5E16.8', data=tmp_data)

        # Next read bond equil value
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, numbnd)
        self.add_flag('BOND_EQUIL_VALUE', '5E16.8', data=tmp_data)

        # Next read angle force constant
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, numang)
        self.add_flag('ANGLE_FORCE_CONSTANT', '5E16.8', data=tmp_data)

        # Next read the angle equilibrium value
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, numang)
        self.add_flag('ANGLE_EQUIL_VALUE', '5E16.8', data=tmp_data)

        # Next read the dihedral force constant
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, nptra)
        self.add_flag('DIHEDRAL_FORCE_CONSTANT', '5E16.8', data=tmp_data)

        # Next read dihedral periodicity
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, nptra)
        self.add_flag('DIHEDRAL_PERIODICITY', '5E16.8', data=tmp_data)

        # Next read the dihedral phase
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, nptra)
        self.add_flag('DIHEDRAL_PHASE', '5E16.8', data=tmp_data)

        # Next read SOLTY (?)
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, natyp)
        self.add_flag('SOLTY', '5E16.8', data=tmp_data)

        # Next read lennard jones acoef and bcoef
        numvals = ntypes * (ntypes + 1) / 2
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, numvals)
        self.add_flag('LENNARD_JONES_ACOEF', '5E16.8', data=tmp_data)
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, numvals)
        self.add_flag('LENNARD_JONES_BCOEF', '5E16.8', data=tmp_data)

        # Next read bonds including hydrogen
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, nbonh*3)
        self.add_flag('BONDS_INC_HYDROGEN', '10I8', data=tmp_data)

        # Next read bonds without hydrogen
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, nbona*3)
        self.add_flag('BONDS_WITHOUT_HYDROGEN', '10I8', data=tmp_data)

        # Next read angles including hydrogen
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, ntheth*4)
        self.add_flag('ANGLES_INC_HYDROGEN', '10I8', data=tmp_data)

        # Next read angles without hydrogen
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, ntheta*4)
        self.add_flag('ANGLES_WITHOUT_HYDROGEN', '10I8', data=tmp_data)

        # Next read dihdrals including hydrogen
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, nphih*5)
        self.add_flag('DIHEDRALS_INC_HYDROGEN', '10I8', data=tmp_data)

        # Next read dihedrals without hydrogen
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, nphia*5)
        self.add_flag('DIHEDRALS_WITHOUT_HYDROGEN', '10I8', data=tmp_data)

        # Next read the excluded atoms list
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, nex)
        self.add_flag('EXCLUDED_ATOMS_LIST', '10I8', data=tmp_data)

        # Next read the hbond terms
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, nphb)
        self.add_flag('HBOND_ACOEF', '5E16.8', data=tmp_data)
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, nphb)
        self.add_flag('HBOND_BCOEF', '5E16.8', data=tmp_data)
        tmp_data, line_idx = read_float(line_idx, prmtop_lines, nphb)
        self.add_flag('HBCUT', '5E16.8', data=tmp_data)

        # Next read amber atom type
        tmp_data, line_idx = read_string(line_idx, prmtop_lines, natom)
        self.add_flag('AMBER_ATOM_TYPE', '20a4', data=tmp_data)

        # Next read tree chain classification
        tmp_data, line_idx = read_string(line_idx, prmtop_lines, natom)
        self.add_flag('TREE_CHAIN_CLASSIFICATION', '20a4', data=tmp_data)

        # Next read the join array
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, natom)
        self.add_flag('JOIN_ARRAY', '10I8', data=tmp_data)

        # Next read the irotat array
        tmp_data, line_idx = read_integer(line_idx, prmtop_lines, natom)
        self.add_flag('IROTAT', '10I8', data=tmp_data)

        # Now do PBC stuff
        if self.parm_data['POINTERS'][PrmtopPointers.IFBOX]:
            # Solvent pointers
            tmp_data, line_idx = read_integer(line_idx, prmtop_lines, 3)
            self.add_flag('SOLVENT_POINTERS', '10I8', data=tmp_data)
            nspm = tmp_data[1]

            # Atoms per molecule
            tmp_data, line_idx = read_integer(line_idx, prmtop_lines, nspm)
            self.add_flag('ATOMS_PER_MOLECULE', '10I8', data=tmp_data)

            # Box dimensions
            tmp_data, line_idx = read_float(line_idx, prmtop_lines, 4)
            self.add_flag('BOX_DIMENSIONS', '5E16.8', data=tmp_data)

        # Now do CAP stuff
        if self.parm_data['POINTERS'][PrmtopPointers.IFCAP]:
            # CAP_INFO
            tmp_data, line_idx = read_integer(line_idx, prmtop_lines, 1)
            self.add_flag('CAP_INFO', '10I8', data=tmp_data)
            tmp_data, line_idx = read_integer(line_idx, prmtop_lines, 4)
            self.add_flag('CAP_INFO2', '10I8', data=tmp_data)
        # end if self.parm_data['POINTERS'][IFCAP]

    #===================================================

    def set_version(self):
        """ Sets the version string """
        now = datetime.datetime.now()
        self.version = (
                '%%VERSION  VERSION_STAMP = V0001.000  DATE = %02d/%02d/%02d  '
                '%02d:%02d:%02d' % (now.month, now.day, now.year % 100,
                                    now.hour, now.minute, now.second)
        )

    #===================================================

    def write_parm(self, name):
        """
        Writes the current data in parm_data into a new topology file with
        the given name

        Parameters
        ----------
        name : str or file-like
            Name of the file to write the topology file to or file-like object to write
        """
        # now that we know we will write the new prmtop file, open the new file
        if isinstance(name, str):
            new_prm = genopen(name, 'w')
            own_handle = True
        else:
            new_prm = name
            own_handle = False
        try:
            # get current time to put into new prmtop file if we had a %VERSION
            self.set_version()
            # convert charges back to amber charges...
            if 'CTITLE' in self.parm_data:
                CHARGE_SCALE = CHARMM_ELECTROSTATIC
            else:
                CHARGE_SCALE = AMBER_ELECTROSTATIC

            if self.charge_flag in self.parm_data.keys():
                for i in range(len(self.parm_data[self.charge_flag])):
                    self.parm_data[self.charge_flag][i] *= CHARGE_SCALE
            # write version to top of prmtop file
            new_prm.write('%s\n' % self.version)

            # write data to prmtop file, inserting blank line if it's an empty field
            for flag in self.flag_list:
                new_prm.write('%%FLAG %s\n' % flag)
                # Insert any comments before the %FORMAT specifier
                for comment in self.parm_comments[flag]:
                    new_prm.write('%%COMMENT %s\n' % comment)
                new_prm.write('%%FORMAT(%s)\n' % self.formats[flag])
                if len(self.parm_data[flag]) == 0: # empty field...
                    new_prm.write('\n')
                    continue
                self.formats[flag].write(self.parm_data[flag], new_prm)
        finally:
            if own_handle:
                new_prm.close()

        if self.charge_flag in self.parm_data.keys():
            # Convert charges back to electron-units
            for i in range(len(self.parm_data[self.charge_flag])):
                self.parm_data[self.charge_flag][i] /= CHARGE_SCALE

    #===================================================

    def add_flag(self, flag_name, flag_format, data=None, num_items=-1,
                 comments=None, after=None):
        """
        Adds a new flag with the given flag name and Fortran format string and
        initializes the array with the values given, or as an array of 0s
        of length num_items

        Parameters
        ----------
        flag_name : str
            Name of the flag to insert. It is converted to all upper case
        flag_format : str
            Fortran format string representing how the data in this section
            should be written and read. Do not enclose in ()
        data : list=None
            Sequence with data for the new flag. If None, a list of zeros of
            length ``num_items`` (see below) is given as a holder
        num_items : int=-1
            Number of items in the section. This variable is ignored if a set of
            data are given in `data`
        comments : list of str=None
            List of comments to add to this section
        after : str=None
            If provided, the added flag will be added after the one with the
            name given to `after`. If this flag does not exist, IndexError will
            be raised

        Raises
        ------
        AmberError if flag already exists
        IndexError if the ``after`` flag does not exist
        """
        if flag_name in self.parm_data:
            raise AmberError('%s already exists' % (flag_name))
        if after is not None:
            after = after.upper()
            if not after in self.flag_list:
                raise IndexError('%s not found in topology flag list' % after)
            # If the 'after' flag is the last one, just append
            if self.flag_list[-1] == after:
                self.flag_list.append(flag_name.upper())
            else:
                # Otherwise find the index and add it after
                idx = self.flag_list.index(after) + 1
                self.flag_list.insert(idx, flag_name.upper())
        else:
            self.flag_list.append(flag_name.upper())
        self.formats[flag_name.upper()] = FortranFormat(flag_format)
        if data is not None:
            self.parm_data[flag_name.upper()] = list(data)
        else:
            if num_items < 0:
                raise AmberError("If you do not supply prmtop data, num_items "
                                 "must be non-negative!")
            self.parm_data[flag_name.upper()] = [0 for i in range(num_items)]
        if comments is not None:
            if isinstance(comments, str):
                comments = [comments]
            else:
                comments = list(comments)
            self.parm_comments[flag_name.upper()] = comments
        else:
            self.parm_comments[flag_name.upper()] = []

    #===================================================

    def delete_flag(self, flag_name):
        """ Removes a flag from the topology file """
        flag_name = flag_name.upper()
        if flag_name in self.flag_list:
            del self.flag_list[self.flag_list.index(flag_name)]
        if flag_name in self.parm_comments:
            del self.parm_comments[flag_name]
        if flag_name in self.formats:
            del self.formats[flag_name]
        if flag_name in self.parm_data:
            del self.parm_data[flag_name]

    #===================================================

    def __getstate__(self):
        return dict(parm_data=self.parm_data, flag_list=self.flag_list,
                    formats=self.formats, parm_comments=self.parm_comments,
                    charge_flag=self.charge_flag, version=self.version,
                    name=self.name)

    def __setstate__(self, d):
        self.__dict__ = d
