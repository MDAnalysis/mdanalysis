"""
Tool for parsing and writing OFF library files to and from dictionaries of
ResidueTemplate objects
"""
from collections import OrderedDict
from contextlib import closing
import numpy as np
from ..topologyobjects import Atom
from ..constants import RAD_TO_DEG
from ..exceptions import AmberWarning
from ..formats.registry import FileFormatType
from ..modeller.residue import ResidueTemplate, ResidueTemplateContainer
from ..modeller.residue import PROTEIN, NUCLEIC, SOLVENT, UNKNOWN
from .. import periodic_table as pt
from ..utils.io import genopen
import re
import warnings

class AmberOFFLibrary(metaclass=FileFormatType):
    """
    Class containing static methods responsible for parsing and writing OFF
    libraries
    """
    #===================================================

    # Useful regexes
    _headerre = re.compile(r'!!index *array *str')
    _resre = re.compile(r'\s*"(\S*?)"\s*$')
    _sec1re = re.compile(r'!entry\.(\S*?)\.unit\.atoms *table *str *name *str'
                         r' *type *int *typex *int *resx *int *flags *int'
                         r' *seq *int *elmnt *dbl *chg')
    _sec2re = re.compile(r'!entry\.(\S*?)\.unit\.atomspertinfo *table *str'
                         r' *pname *str *ptype *int *ptypex *int *pelmnt'
                         r' *dbl *pchg')
    _sec3re = re.compile(r'!entry\.(\S*?)\.unit\.boundbox *array *dbl')
    _sec4re = re.compile(r'!entry\.(\S*?)\.unit\.childsequence *single *int')
    _sec5re = re.compile(r'!entry\.(\S*?)\.unit\.connect *array *int')
    _sec6re = re.compile(r'!entry\.(\S*?)\.unit\.connectivity *table *int'
                         r' *atom1x *int *atom2x *int *flags')
    _sec7re = re.compile(r'!entry\.(\S*?)\.unit\.hierarchy *table *str'
                         r' *abovetype *int *abovex *str *belowtype *int'
                         r' *belowx')
    _sec8re = re.compile(r'!entry\.(\S*?)\.unit\.name *single *str')
    _sec9re = re.compile(r'!entry\.(\S*?)\.unit\.positions *table *dbl *x'
                         r' *dbl *y *dbl *z')
    _sec10re = re.compile(r'!entry\.(\S*?)\.unit\.residueconnect *table'
                          r' *int *c1x *int *c2x *int *c3x *int *c4x *int'
                          r' *c5x *int *c6x')
    _sec11re = re.compile(r'!entry\.(\S*?)\.unit\.residues *table *str *name'
                          r' *int *seq *int *childseq *int *startatomx *str'
                          r' *restype *int *imagingx')
    _sec12re = re.compile(r'!entry\.(\S*?)\.unit\.residuesPdbSequenceNumber'
                          r' *array *int')
    _sec13re = re.compile(r'!entry\.(\S*?)\.unit\.solventcap *array *dbl')
    _sec14re = re.compile(r'!entry\.(\S*?)\.unit\.velocities *table *dbl *x'
                          r' *dbl *y *dbl *z')

    #===================================================

    @classmethod
    def id_format(cls, filename):
        """ Sees if an open file is an OFF library file.

        Parameters
        ----------
        filename : str
            The name of the file to see if it is an OFF file format

        Returns
        -------
        is_fmt : bool
            True if it is recognized as OFF, False otherwise
        """
        with closing(genopen(filename, 'r')) as f:
            if cls._headerre.match(f.readline()):
                return True
            return False

    #===================================================

    @classmethod
    def parse(cls, filename):
        """ Parses an Amber OFF library

        Parameters
        ----------
        filename : str or file-like iterable
            The file name or file object to parse. If it is an iterable, it will
            be exhausted

        Returns
        -------
        residues : OrderedDict {str : :class:`ResidueTemplate`}
            Dictionary pairing residue names with their :class:`ResidueTemplate`
            objects

        Raises
        ------
        ValueError if the first line does not match the file format. This line
        will be consumed

        IOError if filename is the name of a file that does not exist

        RuntimeError if EOF is reached prematurely or other formatting issues
        found
        """
        if isinstance(filename, str):
            fileobj = genopen(filename, 'r')
            own_handle = True
        else:
            fileobj = filename
            own_handle = False
        # Now parse the library file
        line = fileobj.readline()
        if not cls._headerre.match(line):
            raise ValueError('Unrecognized OFF file format')
        # Build the return value
        residues = OrderedDict()
        # Pull a list of all the residues we expect to find
        line = fileobj.readline()
        rematch = cls._resre.match(line)
        while rematch and line:
            name = rematch.groups()[0]
            residues[name] = None
            line = fileobj.readline()
            rematch = cls._resre.match(line)
        if not line:
            raise RuntimeError('Unexpected EOF in Amber OFF library')
        # Now make sure we have the next expected line
        while line:
            if not line.strip():
                line = fileobj.readline()
                continue
            rematch = cls._sec1re.match(line)
            if not rematch:
                raise RuntimeError('Expected atoms table not found')
            name = rematch.groups()[0]
            residues[name] = cls._parse_residue(fileobj, name)
            line = fileobj.readline()

        if own_handle:
            fileobj.close()

        return residues

    #===================================================

    @classmethod
    def _parse_residue(cls, fileobj, name):
        """
        Parses the residue information out of the OFF file assuming the file
        is pointed at the first line of an atoms table section of the OFF file

        Parameters
        ----------
        fileobj : file-like
            Assumed to be open for read, this file is parsed until the *next*
            atom table is read
        name : str
            The name of the residue being processed right now
        """
        container = ResidueTemplateContainer(name)
        nres = 1
        templ = ResidueTemplate(name)
        line = fileobj.readline()
        while line[0] != '!':
            nam, typ, typx, resx, flags, seq, elmnt, chg = line.split()
            nam = _strip_enveloping_quotes(nam)
            typ = _strip_enveloping_quotes(typ)
            typx = int(typx)
            resx = int(resx)
            flags = int(flags)
            seq = int(seq)
            elmnt = int(elmnt)
            chg = float(chg)
            atom = Atom(atomic_number=elmnt, type=typ, name=nam, charge=chg)
            if resx == nres + 1:
                container.append(templ)
                nres += 1
                templ = ResidueTemplate(name)
            templ.add_atom(atom)
            line = fileobj.readline()
            # Skip blank lines
            while line and not line.strip():
                line = fileobj.readline()
        container.append(templ)
        if nres > 1:
            start_atoms = []
            runsum = 0
            for res in container:
                start_atoms.append(runsum)
                runsum += len(res)
        # Make sure we get the next section
        rematch = cls._sec2re.match(line)
        if not rematch:
            raise RuntimeError('Expected pertinfo table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError(
                f'Found residue {rematch.groups()[0]} while processing residue {name}'
            )
        line = fileobj.readline()
        while line[0] != '!':
            if not line:
                raise RuntimeError('Unexpected EOF in Amber OFF library')
            # Not used, just skip
            # TODO sanity check
            line = fileobj.readline()
        rematch = cls._sec3re.match(line)
        if not rematch:
            raise RuntimeError('Expected boundbox table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError(
                f'Found residue {rematch.groups()[0]} while processing residue {name}'
            )
        # Only 5 lines
        try:
            hasbox = float(fileobj.readline().strip())
            angle = float(fileobj.readline().strip())
            a = float(fileobj.readline().strip())
            b = float(fileobj.readline().strip())
            c = float(fileobj.readline().strip())
        except ValueError:
            raise RuntimeError('Error processing boundbox table entries')
        else:
            if hasbox > 0:
                if angle < 3.15:
                    # No box is this acute -- must be in radians
                    angle *= RAD_TO_DEG
                container.box = [a, b, c, angle, angle, angle]
        # Get the child sequence entry
        line = fileobj.readline()
        rematch = cls._sec4re.match(line)
        if not rematch:
            raise RuntimeError('Expected childsequence table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError(
                f"Found residue {rematch.groups()[0]} while processing residue {name}"
            )
        n = int(fileobj.readline().strip())
        if nres + 1 != n:
            warnings.warn(
                f"Unexpected childsequence ({n}); expected {nres+1} for residue {name}",
                AmberWarning,
            )
        elif not isinstance(templ, ResidueTemplate) and n != len(templ) + 1:
            raise RuntimeError("child sequence must be # of residues in the unit + 1")
        # Get the CONNECT array to set head and tail
        line = fileobj.readline()
        rematch = cls._sec5re.match(line)
        if not rematch:
            raise RuntimeError('Expected connect array not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError(
                f"Found residue {rematch.groups()[0]} while processing residue {name}"
            )
        try:
            head = int(fileobj.readline().strip())
            tail = int(fileobj.readline().strip())
        except ValueError:
            raise RuntimeError('Error processing connect table entries')
        if head > 0 and nres == 1:
            templ.head = templ[head-1]
        elif head > 0 and nres > 1:
            if head < sum((len(r) for r in container)):
                raise RuntimeError('HEAD on multi-residue unit not supported')
        if tail > 0 and nres == 1:
            templ.tail = templ[tail-1]
        elif tail > 0 and nres > 1:
            if tail < sum((len(r) for r in container)):
                warnings.warn(
                    f'TAIL on multi-residue unit not supported ({name}). Ignored...', AmberWarning
                )
        # Get the connectivity array to set bonds
        line = fileobj.readline()
        if len(templ.atoms) > 1:
            rematch = cls._sec6re.match(line)
            if not rematch:
                raise RuntimeError('Expected connectivity table not found')
            elif rematch.groups()[0] != name:
                raise RuntimeError(
                    f"Found residue {rematch.groups()[0]} while processing residue {name}"
                )
            line = fileobj.readline()
            while line[0] != '!':
                i, j, flag = line.split()
                line = fileobj.readline()
                if nres > 1:
                    # Find which residue we belong in
                    i = int(i) - 1
                    j = int(j) - 1
                    for ii, idx in enumerate(start_atoms):
                        if idx > i:
                            ii -= 1
                            break
                    start_idx = start_atoms[ii]
                    container[ii].add_bond(i-start_idx, j-start_idx)
                else:
                    templ.add_bond(int(i)-1, int(j)-1)
        # Get the hierarchy table
        rematch = cls._sec7re.match(line)
        if not rematch:
            raise RuntimeError('Expected hierarchy table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError(
                f"Found residue {rematch.groups()[0]} while processing residue {name}"
            )
        line = fileobj.readline()
        while line[0] != '!':
            # Skip this section... not used
            # TODO turn this into a sanity check
            line = fileobj.readline()
        # Get the unit name
        rematch = cls._sec8re.match(line)
        if not rematch:
            raise RuntimeError('Expected unit name string not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError(
                f"Found residue {rematch.groups()[0]} while processing residue {name}"
            )
        fileobj.readline() # Skip this... not used
        line = fileobj.readline()
        # Get the atomic positions
        rematch = cls._sec9re.match(line)
        if not rematch:
            raise RuntimeError('Expected unit positions table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError(
                f"Found residue {rematch.groups()[0]} while processing residue {name}"
            )
        for res in container:
            for atom in res:
                x, y, z = fileobj.readline().split()
                atom.xx, atom.xy, atom.xz = float(x), float(y), float(z)
        line = fileobj.readline()
        # Get the residueconnect table
        rematch = cls._sec10re.match(line)
        if not rematch:
            raise RuntimeError('Expected unit residueconnect table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError(
                f"Found residue {rematch.groups()[0]} while processing residue {name}"
            )
        for i in range(nres):
            c1,c2,c3,c4,c5,c6 = (int(x) for x in fileobj.readline().split())
            if (c1 > 0 and templ.head is not None and
                    templ.head is not templ[c1-1]):
                raise RuntimeError('HEAD atom is not connect0')
            if (c2 > 0 and templ.tail is not None and
                    templ.tail is not templ[c2-1]):
                raise RuntimeError('TAIL atom is not connect1')
            for i in (c3, c4, c5, c6):
                if i == 0: continue
                templ.connections.append(templ[i-1])
        # Get the residues table
        line = fileobj.readline()
        rematch = cls._sec11re.match(line)
        if not rematch:
            raise RuntimeError('Expected unit residues table not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError(
                f"Found residue {rematch.groups()[0]} while processing residue {name}"
            )
        for i in range(nres):
            resname, id, next, start, typ, img = fileobj.readline().split()
            resname = _strip_enveloping_quotes(resname)
            id = int(id)
            start = int(start)
            next = int(next)
            typ = _strip_enveloping_quotes(typ)
            img = int(img)
            if next - start != len(container[i]):
                warnings.warn(
                    f'residue table predicted {next - start}, not {len(container[i])} atoms for residue {name}',
                    AmberWarning,
                )
            if typ == 'p':
                container[i].type = PROTEIN
            elif typ == 'n':
                container[i].type = NUCLEIC
            elif typ == 'w':
                container[i].type = SOLVENT
            elif typ != '?':
                warnings.warn('Unknown residue type "%s"' % typ, AmberWarning)
            if nres > 1:
                container[i].name = resname
        # Get the residues sequence table
        line = fileobj.readline()
        rematch = cls._sec12re.match(line)
        if not rematch:
            raise RuntimeError('Expected residue sequence number not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError(
                f"Found residue {rematch.groups()[0]} while processing residue {name}"
            )
        for i in range(nres):
            #TODO sanity check
            fileobj.readline()
        line = fileobj.readline()
        # Get the solventcap array
        rematch = cls._sec13re.match(line)
        if not rematch:
            raise RuntimeError('Expected unit solventcap array not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError(
                f"Found residue {rematch.groups()[0]} while processing residue {name}"
            )
        # Ignore the solvent cap
        fileobj.readline()
        fileobj.readline()
        fileobj.readline()
        fileobj.readline()
        fileobj.readline()
        # Velocities
        line = fileobj.readline()
        rematch = cls._sec14re.match(line)
        if not rematch:
            raise RuntimeError('Expected unit solventcap array not found')
        elif rematch.groups()[0] != name:
            raise RuntimeError(
                f"Found residue {rematch.groups()[0]} while processing residue {name}"
            )
        for res in container:
            for atom in res:
                vx, vy, vz = (float(x) for x in fileobj.readline().split())
                atom.vx, atom.vy, atom.vz = vx, vy, vz

        if nres > 1:
            return container
        return templ

    #===================================================

    @classmethod
    def write(cls, lib, dest):
        """ Writes a dictionary of ResidueTemplate units to a file in OFF format

        Parameters
        ----------
        lib : dict {str : :class:`ResidueTemplate`}
            Items can be either :class:`ResidueTemplate` or
            :class:`ResidueTemplateContainer` instances
        dest : str or file-like
            Either a file name or a file-like object to write the file to
        """
        own_handle = False
        if not hasattr(dest, 'write'):
            dest = genopen(dest, 'w')
            own_handle = True
        # Write the residues in alphabetical order
        names = sorted(lib.keys())
        dest.write('!!index array str\n')
        for name in names:
            dest.write(' "%s"\n' % name)
        for name in names:
            cls._write_residue(dest, lib[name])

        if own_handle: dest.close()

    #===================================================

    @staticmethod
    def _write_residue(dest, res):
        """ Writes a residue to an open file handle

        Parameters
        ----------
        dest : file-like
            File object to write the residue information to
        res : :class:`ResidueTemplate` or :class:`ResidueTemplateContainer`
            The residue template (or template container) to write to the file
        """
        if isinstance(res, ResidueTemplate):
            # Put it into a template container with the same name
            tmp = ResidueTemplateContainer(res.name)
            tmp.append(res)
            res = tmp
        dest.write('!entry.%s.unit.atoms table  str name  str type  int typex  '
                   'int resx  int flags  int seq  int elmnt  dbl chg\n' %
                   res.name)
        for i, r in enumerate(res):
            for atom in r:
                dest.write(' "%s" "%s" 0 %d 131072 %d %d %.6f\n' % (atom.name,
                           atom.type, i+1, atom.idx+1, atom.atomic_number,
                           atom.charge))
        dest.write('!entry.%s.unit.atomspertinfo table  str pname  str ptype  '
                   'int ptypex  int pelmnt  dbl pchg\n' % res.name)
        for r in res:
            for atom in r:
                dest.write(' "%s" "%s" 0 -1 0.0\n' % (atom.name, atom.type))
        dest.write('!entry.%s.unit.boundbox array dbl\n' % res.name)
        if res.box is None:
            dest.write((' -1.000000\n' + ' 0.0\n' * 4))
        else:
            dest.write(' 1.000000\n')
            if res.box[3] == res.box[4] == res.box[5]:
                dest.write(' %f\n' % res.box[3])
            else:
                raise RuntimeError('Cannot write boxes with different angles')
            dest.write(' %f\n' % res.box[0])
            dest.write(' %f\n' % res.box[1])
            dest.write(' %f\n' % res.box[2])
        dest.write('!entry.%s.unit.childsequence single int\n %d\n' %
                   (res.name, len(res)+1))
        dest.write('!entry.%s.unit.connect array int\n' % res.name)
        if len(res) > 1:
            dest.write(' 0\n 0\n')
        else:
            if res[0].head is not None:
                dest.write(' %d\n' % (res[0].head.idx + 1))
            else:
                dest.write(' 0\n')
            if res[0].tail is not None:
                dest.write(' %d\n' % (res[0].tail.idx + 1))
            else:
                dest.write(' 0\n')
        if any(len(r) > 1 for r in res):
            dest.write('!entry.%s.unit.connectivity table  int atom1x  '
                       'int atom2x  int flags\n' % res.name)
            base = 1
            for r in res:
                for bond in r.bonds:
                    dest.write(' %d %d 1\n' % (bond.atom1.idx+base,
                                               bond.atom2.idx+base))
                base += len(r)
        dest.write('!entry.%s.unit.hierarchy table  str abovetype  int '
                   'abovex  str belowtype  int belowx\n' % res.name)
        c = 1
        for i, r in enumerate(res):
            dest.write(' "U" 0 "R" %d\n' % (i+1))
            for atom in r:
                dest.write(' "R" %d "A" %d\n' % (i+1, c))
                c += 1
        dest.write('!entry.%s.unit.name single str\n' % res.name)
        dest.write(' "%s"\n' % res.name)
        dest.write('!entry.%s.unit.positions table  dbl x  dbl y  dbl z\n' %
                   res.name)
        for r in res:
            for atom in r:
                dest.write(' %.6g %.6g %.6g\n' % (atom.xx, atom.xy, atom.xz))
        dest.write('!entry.%s.unit.residueconnect table  int c1x  int c2x  '
                   'int c3x  int c4x  int c5x  int c6x\n' % res.name)
        c = 1
        for r in res:
            # Make the CONECT1 and 0 default to first and last atom so that the
            # TREE gets set correctly by tleap. Not used for anything else...
            conn = [c, c+len(r)-1, 0, 0, 0, 0]
            if r.head is not None: conn[0] = r.head.idx + 1
            if r.tail is not None: conn[1] = r.tail.idx + 1
            for i, at in enumerate(r.connections[:4]):
                conn[i+2] = at.idx + 1
            dest.write(' %d %d %d %d %d %d\n' % tuple(conn))
            c += len(r)
        dest.write('!entry.%s.unit.residues table  str name  int seq  int '
                   'childseq  int startatomx  str restype  int imagingx\n' %
                   res.name)
        c = 1
        for i, r in enumerate(res):
            if r.type is PROTEIN:
                typ = 'p'
            elif r.type is NUCLEIC:
                typ = 'n'
            elif r.type is SOLVENT:
                typ='w'
            elif r.type is UNKNOWN:
                typ='?'
            else:
                warnings.warn('Unrecognized residue type %r' % r.type,
                              AmberWarning)
                typ = '?'
            dest.write(' "%s" %d %d %d "%s" %d\n' % (r.name, i+1, 1+len(r), c,
                       typ, _imaging_atom(r)+c))
            c += len(r)
        dest.write('!entry.%s.unit.residuesPdbSequenceNumber array int\n' %
                   res.name)
        for i, r in enumerate(res):
            if len(res) == 1:
                dest.write(' 0\n')
            else:
                dest.write(' %d\n' % (i+1))
        dest.write('!entry.%s.unit.solventcap array dbl\n' % res.name)
        dest.write(' -1.000000\n' + ' 0.0\n' * 4)
        dest.write('!entry.%s.unit.velocities table  dbl x  dbl y  dbl z\n' %
                   res.name)
        for r in res:
            for atom in r:
                try:
                    s = ' %g %g %g\n' % (atom.vx, atom.vy, atom.vz)
                except AttributeError:
                    dest.write(' 0.0 0.0 0.0\n')
                else:
                    dest.write(s)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Helper routines
def _strip_enveloping_quotes(inp):
    """ Strips the quotation marks enveloping a string """
    if inp[0] == inp[-1] == '"' or inp[0] == inp[-1] == "'":
        return inp[1:-1]
    return inp

def _imaging_atom(res):
    """
    Determines the imaging atom for the residue. If all atoms are hydrogen
    except 1, it is the heavy atom. Otherwise, it is the atom *closest* to the
    COM of the residue
    """
    from ..geometry import center_of_mass
    #TODO implement the docstring
    found_heavy = False
    heavy_idx = -1
    for i, atom in enumerate(res):
        if atom.atomic_number > 1:
            heavy_idx = i
            if found_heavy: break
            found_heavy = True
    else:
        if heavy_idx != -1:
            return heavy_idx
        return 0 # No heavy atoms?? No imaging atom, then.
    coords = res.coordinates.reshape((len(res), 3))
    masses = np.zeros(len(res))
    for i, atom in enumerate(res):
        if atom.mass == 0:
            masses[i] = pt.Mass[pt.Element[atom.atomic_number]]
        else:
            masses[i] = atom.mass
    com = center_of_mass(coords, masses)
    diff = coords - com
    return np.argmin((diff * diff).sum(axis=1))
