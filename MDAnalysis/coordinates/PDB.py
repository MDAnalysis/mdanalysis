# $Id: PDB.py 101 2008-05-18 13:19:06Z orbeckst $
"""
PDB structure files in MDAnalysis --- :mod:`MDAnalysis.coordinates.PDB`
========================================================================

MDAnalysis reads coordinates from PDB files and additional optional
data such as B-factors. It is also possible to substitute a PDB file
instead of PSF file in order to define the list of atoms (but no
connectivity information will be  available in this case).

The :mod:`PDB` module makes heavy use of Biopython's :mod:`Bio.PDB`:

  Hamelryck, T., Manderick, B. (2003) PDB parser and structure class 
  implemented in Python. Bioinformatics, 19, 2308-2310.

  http://biopython.org

but replaces the standard PDB file parser with one that uses the
:class:`MDAnalysis.coordinates.pdb.extensions.SloppyStructureBuilder`
to cope with very large pdb files as commonly encountered in MD
simulations.
"""
from __future__ import with_statement

try:
    # BioPython is overkill but potentially extensible (altLoc etc)
    import Bio.PDB
except ImportError:
    raise ImportError("No PDB I/O functionality. Install biopython.")

import numpy

import MDAnalysis.core
import MDAnalysis.core.util as util
import base
import base
import pdb.extensions

from MDAnalysis.topology.core import guess_atom_element

import warnings

class Timestep(base.Timestep):
	@property
	def dimensions(self):
	        """unitcell dimensions (`A, B, C, alpha, beta, gamma`)

		- `A, B, C` are the lengths of the primitive cell vectors `e1, e2, e3`
		- `alpha` = angle(`e1, e2`)
		- `beta` = angle(`e1, e3`)
		- `gamma` = angle(`e2, e3`)
		"""
		# Layout of unitcell is [A,B,C,90,90,90] with the primitive cell vectors
		return self._unitcell
    
class PDBReader(base.Reader):
    """Read a pdb file into a BioPython pdb structure.

    The coordinates are also supplied as one numpy array and wrapped
    into a Timestep object; attributes are set so that the PDBReader
    object superficially resembles the DCDReader object.
    """
    format = 'PDB'
    units = {'time': None, 'length': 'Angstrom'}
    _Timestep = Timestep

    def __init__(self, pdbfilename, convert_units=None, **kwargs):
        self.pdbfilename = pdbfilename
        self.filename = self.pdbfilename
        if convert_units is None:
            convert_units = MDAnalysis.core.flags['convert_gromacs_lengths']
        self.convert_units = convert_units  # convert length and time to base units

        pdb_id = "0UNK"
        self.pdb = pdb.extensions.get_structure(pdbfilename, pdb_id)
        pos = numpy.array([atom.coord for atom in self.pdb.get_atoms()])
        self.numatoms = pos.shape[0]
        self.numframes = 1
        self.fixed = 0          # parse B field for fixed atoms?
        self.skip = 1
        self.periodic = False
        self.delta = 0
        self.skip_timestep = 1
	#self.ts._unitcell[:] = ??? , from CRYST1?
	self.ts = self._Timestep(pos)
	del pos
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)             # in-place !
            
    def get_bfactors(self):
        """Return an array of bfactors (tempFactor) in atom order."""
        warnings.warn("get_bfactors() will be removed in MDAnalysis 0.8; "
                      "use AtomGroup.bfactors [which will become AtomGroup.bfactors()]", 
                      DeprecationWarning)
        return numpy.array([a.get_bfactor() for a in self.pdb.get_atoms()])

    def Writer(self, filename, **kwargs):
        """Returns a strict PDBWriter for *filename*.

        :Arguments:
          *filename*
              filename of the output PDB file

        :Returns: :class:`PDBWriter`

        .. Note:: This :class:`PDBWriter` 's :meth:`~PDBWriter.write` method always requires a 
                  :class:`Timestep` as an argument (it is not optional anymore when the Writer
                  is obtained through this method of :class:`PDBReader`.)
        """
        # This is messy; we cannot get a universe from the Reader, which would be
        # also needed to be fed to the PDBWriter (which is a total mess...).
        # Hence we ignore the problem and document it in the doc string... --- the
        # limitation is simply that PDBWriter.write() must always be called with an argument.
        kwargs['BioPDBstructure'] = self.pdb   # make sure that this Writer is 
        kwargs.pop('universe', None)           # always linked to this reader, don't bother with Universe
        return PDBWriter(filename, **kwargs)

    def __len__(self):
        return self.numframes
    def __iter__(self):
        def iterPDB():
            yield self.ts   # just a single frame available
            raise StopIteration
        return iterPDB()
    def __getitem__(self, frame):
        if frame != 0:
            raise IndexError('PDBReader only contains a single frame at index 0')
        return self.ts


class PDBWriter(base.Writer):
    """Write out the current time step as a pdb file.

    This is not cleanly implemented at the moment. One must supply a
    universe, even though this is nominally an optional argument. The
    class behaves slightly differently depending on if the structure
    was loaded from a PDB (then the full-fledged :mod:`Bio.PDB` writer is
    used) or if this is really only an atom selection (then a less
    sophistiocated writer is employed).

    .. Note:: The standard PDBWriter can only write the *whole system*. 
      In order to write a selection, use the :class:`PrimitivePDBWriter`, 
      which happens automatically when the
      :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.write` method of a 
      :class:`~MDAnalysis.core.AtomGroup.AtomGroup` instance is used.
    """
    format = 'PDB'
    units = {'time': None, 'length': 'Angstrom'}

    # PDBWriter is a bit more complicated than the DCDWriter in the
    # sense that a DCD frame only contains coordinate information. The
    # PDB contains atom data as well and hence it MUST access the
    # universe. In order to present a unified (and backwards
    # compatible) interface we must keep the universe argument an
    # optional keyword argument even though it really is required.
    
    def __init__(self,pdbfilename,universe=None,multi=False,**kwargs):
        """pdbwriter = PDBWriter(<pdbfilename>,universe=universe,**kwargs)
        :Arguments:
        pdbfilename     filename; if multi=True, embed a %%d formatstring
                        so that write_next_timestep() can insert the frame number
        universe        supply a universe [really REQUIRED; optional only for compatibility]
        multi           False: write a single structure to a single pdb
                        True: write all frames to multiple pdb files
        """
        import Bio.PDB.Structure
        self.universe = universe
        self.PDBstructure = kwargs.pop('BioPDBstructure', None)  # hack for PDBReader.Writer()
        if not self.PDBstructure:
            try:
                self.PDBstructure = universe.trajectory.pdb
            except AttributeError:
                pass
        self.filename = pdbfilename
        self.multi = multi
        if self.multi:
            raise NotImplementedError('Sorry, multi=True does not work yet.')
        if self.PDBstructure is not None and not isinstance(self.PDBstructure,Bio.PDB.Structure.Structure):
            raise TypeError('If defined, PDBstructure must be a Bio.PDB.Structure.Structure, eg '
                            'Universe.trajectory.pdb.')
    def write_next_timestep(self,ts=None):
        self.write(ts)
    def write(self,ts=None):
        """Write timestep as a pdb file.

        If ts=None then we try to get the current one from the universe.
        """        
        if self.PDBstructure is None:
            if self.universe is None:
                warnings.warn("PDBWriter: Not writing frame as neither Timestep nor Universe supplied.")
                return            
            # primitive PDB writing (ignores timestep argument)
            ppw = PrimitivePDBWriter(self.filename)
            ppw.write(self.universe.selectAtoms('all'))
            ppw.close()
        else:
            # full fledged PDB writer
            # Let's cheat and use universe.pdb.pdb: modify coordinates
            # and save...
            if ts is None:
                try:
                    ts = self.universe.trajectory.ts
                except AttributeError:
                    warnings.warn("PDBWriter: Not writing frame as neither universe nor timestep supplied.")
                    return
            if not hasattr(ts, '_pos'):
                raise TypeError("The PDBWriter can only process a Timestep as optional argument, not "
                                "e.g. a selection. Use the PrimitivePDBWriter instead and see the docs.")
            for a,pos in zip(self.PDBstructure.get_atoms(), ts._pos):
                a.set_coord(pos)
            io = pdb.extensions.SloppyPDBIO()
            io.set_structure(self.PDBstructure)
            io.save(self.filename)

    def close_trajectory(self):
        pass    # do nothing, keeps super classe's __del__ happy


class PrimitivePDBReader(base.Reader):
    """PDBReader that reads a PDB-formatted file, no frills.

    Records read:
     - CRYST1 for unitcell A,B,C, alpha,beta,gamm
     - ATOM. HETATM for x,y,z

    http://www.wwpdb.org/documentation/format32/sect9.html

    =============  ============  ===========  =============================================
    COLUMNS        DATA  TYPE    FIELD        DEFINITION
    =============  ============  ===========  =============================================
    1 -  6         Record name   "CRYST1"
    7 - 15         Real(9.3)     a              a (Angstroms).
    16 - 24        Real(9.3)     b              b (Angstroms).
    25 - 33        Real(9.3)     c              c (Angstroms).
    34 - 40        Real(7.2)     alpha          alpha (degrees).
    41 - 47        Real(7.2)     beta           beta (degrees).
    48 - 54        Real(7.2)     gamma          gamma (degrees).

    1 -  6         Record name   "ATOM  "
    7 - 11         Integer       serial       Atom  serial number.
    13 - 16        Atom          name         Atom name.
    17             Character     altLoc       Alternate location indicator. IGNORED
    18 - 21        Residue name  resName      Residue name.
    22             Character     chainID      Chain identifier.
    23 - 26        Integer       resSeq       Residue sequence number.
    27             AChar         iCode        Code for insertion of residues. IGNORED
    31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)     occupancy    Occupancy.
    61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    67 - 76        String        segID        (unofficial CHARMM extension ?)
    77 - 78        LString(2)    element      Element symbol, right-justified. IGNORED
    79 - 80        LString(2)    charge       Charge  on the atom. IGNORED
    =============  ============  ===========  =============================================
    """
    format = 'PDB'
    units = {'time': None, 'length': 'Angstrom'}
    _Timestep = Timestep

    def __init__(self, filename, convert_units=None, **kwargs):
        """Read coordinates from *filename*.

        *filename* can be a gzipped or bzip2ed compressed PDB file.
        """
        self.filename = filename
        if convert_units is None:
            convert_units = MDAnalysis.core.flags['convert_gromacs_lengths']
        self.convert_units = convert_units  # convert length and time to base units

        coords = []
        atoms = []
        unitcell = numpy.zeros(6, dtype=numpy.float32)
        with util.openany(filename, 'r') as pdbfile:
            for line in pdbfile:
                def _c(start,stop,typeclass=float):
                    return self._col(line, start, stop, typeclass=typeclass)
                if line[:3] == 'END':
                    break
                if line[:6] == 'CRYST1':
                    A,B,C = _c(7,15), _c(16,24), _c(25,33)
                    alpha,beta,gamma = _c(34,40), _c(41,47), _c(48,54)
                    unitcell[:] = A,B,C, alpha,beta,gamma
                if line[:6] in ('ATOM  ', 'HETATM'):
                    # directly use COLUMNS from PDB spec
                    serial = _c(7,11,int)
                    name = _c(13,16,str).strip()
                    resName = _c(18,21,str).strip()
                    chainID = _c(22,22,str)  # empty chainID is a single space ' '!
                    resSeq = _c(23,26,int)
                    x,y,z = _c(31,38), _c(39,46), _c(47,54)
                    occupancy = _c(55,60)
                    tempFactor = _c(61,66)
                    segID = _c(67,76, str).strip()
                    element = _c(77,78, str).strip()
                    coords.append((x,y,z))
                    atoms.append((serial, name, resName, chainID, resSeq, occupancy, tempFactor, segID, element))
        self.numatoms = len(coords)
        self.ts = self._Timestep(numpy.array(coords, dtype=numpy.float32))
        self.ts._unitcell[:] = unitcell
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)             # in-place !
            self.convert_pos_from_native(self.ts._unitcell[:3])    # in-place ! (only lengths)
        self.numframes = 1
        self.fixed = 0
        self.skip = 1
        self.periodic = False
        self.delta = 0
        self.skip_timestep = 1
        # hack for PrimitivePDBParser:
        self._atoms = numpy.rec.fromrecords(atoms, names="serial,name,resName,chainID,resSeq,occupancy,tempFactor,segID,element")

    def _col(self, line, start, stop, typeclass=float):
        """Pick out and convert the columns start-stop.

        Numbering starts at column 1 with *start* and includes *stop*;
        this is the convention used in FORTRAN (and also in the PDB format).

        :Returns: ``typeclass(line[start-1:stop])`` or
                  ``typeclass(0)`` if conversion fails
        """
        x = line[start-1:stop]
        try:
            return typeclass(x)
        except ValueError:
            return typeclass(0)

    def get_bfactors(self):
        """Return an array of bfactors (tempFactor) in atom order."""
        warnings.warn("get_bfactors() will be removed in MDAnalysis 0.8; "
                      "use AtomGroup.bfactors [which will become AtomGroup.bfactors()]", 
                      DeprecationWarning)
        return self._atoms.tempFactor

    def get_occupancy(self):
        """Return an array of occupancies in atom order."""
        return self._atoms.occupancy

    def Writer(self, filename, **kwargs):
        """Returns a permissive (simple) PDBWriter for *filename*.

        :Arguments:
          *filename*
              filename of the output PDB file

        :Returns: :class:`PrimitivePDBWriter`

        """
        return PrimitivePDBWriter(filename, **kwargs)

    def __len__(self):
        return self.numframes
    def __iter__(self):
        yield self.ts  # Just a single frame
        raise StopIteration
    def __getitem__(self, frame):
        if frame != 0:
            raise IndexError('PrimitivePDBReader can only read a single frame at index 0')
        return self.ts
    def _read_next_timestep(self):
        raise IndexError("PrimitivePDBReader can only read a single frame")

class PrimitivePDBWriter(base.Writer):
    """PDB writer that implements a subset of the PDB 3.2 standard.
    http://www.wwpdb.org/documentation/format32/v3.2.html
    """
    #          1         2         3         4         5         6         7         8
    # 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.
    # ATOM__seria nameAres CressI   xxxxxxxxyyyyyyyyzzzzzzzzOCCUPAtempft          elCH
    # ATOM  %5d   %-4s %-3s %4d %1s %8.3f   %8.3f   %8.3f   %6.2f %6.2f           %2s
    #                 %1s  %1s                                                      %2d
    #            =        =      ===                                    ==========
    # ATOM  %5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2d
    # ATOM  %(serial)5d %(name)-4s%(altLoc)1s%(resName)-3s %(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f          %(element)2s%(charge)2d

    # Strict PDB format:
    #fmt = {'ATOM':   "ATOM  %(serial)5d %(name)-4s%(altLoc)1s%(resName)-3s %(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f          %(element)2s%(charge)2d\n",
    # PDB format as used by NAMD/CHARMM: 4-letter resnames and segID, altLoc ignored
    fmt = {'ATOM':   "ATOM  %(serial)5d %(name)-4s %(resName)-4s%(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f      %(segID)-4s%(element)2s%(charge)2d\n",
           'REMARK': "REMARK     %s\n",
           'TITLE':  "TITLE    %s\n",
           'CRYST1': "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
           }
    format = 'PDB'
    units = {'time': None, 'length': 'Angstrom'}
    pdb_coor_limits = {"min":-999.994, "max":9999.994}

    def __init__(self,filename,**kwargs):
        self.filename = util.filename(filename,ext='pdb')
        self.pdb = open(self.filename,'w')

    def close(self):
        self.pdb.close()

    def write(self,selection,frame=None):
        """Write selection at current trajectory frame to file.

        write(selection,frame=FRAME)

        :Arguments:
          *selection*
            a :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
          *frame*
            optionally move to frame *FRAME*

        .. Note:: The first letter of the :attr:`~MDAnalysis.core.AtomGroup.Atom.segid` 
                  is used as the PDB chainID.
        """
        u = selection.universe
        if frame is not None:            
            u.trajectory[frame]  # advance to frame
        else:
            try:
                frame = u.trajectory.ts.frame
            except AttributeError:
                frame = 1   # should catch cases when we are analyzing a single PDB (?)
        
        self.TITLE("FRAME "+str(frame)+" FROM "+str(u.trajectory.filename))
        self.CRYST1(self.convert_dimensions_to_unitcell(u.trajectory.ts))
        atoms = selection.atoms    # make sure to use atoms (Issue 46)
        coor = atoms.coordinates() # can write from selection == Universe (Issue 49)
        
        for i, atom in enumerate(atoms):
            if not self.has_valid_coordinates(coor[i]):
                raise ValueError("PDB files can have maximum coordinate values of 9999.994/-999.994. Problem with particle/atom %d" % (i+1))
                break
            self.ATOM(serial=i+1, name=atom.name.strip(), resName=atom.resname.strip(), resSeq=atom.resid,
                      chainID=atom.segid.strip(), segID=atom.segid.strip(),
                      x=coor[i,0], y=coor[i,1], z=coor[i,2])
            # get bfactor, too, and add to output?
            # 'element' is auto-guessed from atom.name in ATOM()
        self.close()

    def has_valid_coordinates(self, coor_list):
        """
        @Input: numpy array of 3 (x, y, z) coordinates for a particle/atom
        @Output: boolean True or False, True 
        True, if all values are within 9999.994/-999.994, as required for PDBs.
        """
        for coor in coor_list:
            # this expression could be performance-poor
            if not self.pdb_coor_limits["min"] < coor < self.pdb_coor_limits["max"]:
                return False        
        return True

    def TITLE(self,*title):
        """Write TITLE record.
        http://www.wwpdb.org/documentation/format32/sect2.html
        """        
        line = " ".join(title)    # should do continuation automatically
        self.pdb.write(self.fmt['TITLE'] % line)

    def REMARK(self,*remark):
        """Write generic REMARK record (without number).
        http://www.wwpdb.org/documentation/format32/remarks1.html
        http://www.wwpdb.org/documentation/format32/remarks2.html
        """
        line = " ".join(remark)
        self.pdb.write(self.fmt['REMARK'] % line)

    def CRYST1(self,dimensions, spacegroup='P 1', zvalue=1):
        """Write CRYST1 record.
        http://www.wwpdb.org/documentation/format32/sect8.html
        """
        self.pdb.write(self.fmt['CRYST1'] % (tuple(dimensions)+(spacegroup, zvalue)))

    def ATOM(self,serial=None,name=None,altLoc=None,resName=None,chainID=None,
             resSeq=None,iCode=None,x=None,y=None,z=None,occupancy=1.0,tempFactor=0.0,
             segID=None,element=None,charge=0):
        """Write ATOM record. 
        http://www.wwpdb.org/documentation/format32/sect9.html
        Only some keword args are optional (altLoc, iCode, chainID), for some defaults are set.

        All inputs are cut to the maximum allowed length. For integer
        numbers the highest-value digits are chopped (so that the
        serial and reSeq wrap); for strings the trailing characters
        are chopped.

        Note: Floats are not checked and can potentially screw up the format.
        """
        for arg in ('serial','name','resName','resSeq','x','y','z',
                    'occupancy','tempFactor','charge'):
            if locals()[arg] is None:
                raise ValueError('parameter '+arg+' must be defined.')
        serial = int(str(serial)[-5:])  # check for overflow here?
        name = name[:4]
        if len(name) < 4:
            name = " "+name   # customary to start in column 14
        altLoc = altLoc or " "
        altLoc= altLoc[:1]
        resName = resName[:4]
        chainID = chainID or ""   # or should we provide a chainID such as 'A'?
        chainID = chainID.strip()[-1:] # take the last character
        resSeq = int(str(resSeq)[-4:]) # check for overflow here?
        iCode = iCode or ""
        iCode = iCode[:1]
        element = element or guess_atom_element(name)  # element == 0|False|None will be guessed
        element = str(element)[:2]            # make sure that is a string for user input
        segID = segID or chainID
        segID = segID[:4]
        self.pdb.write(self.fmt['ATOM'] % vars())        
        
    def __del__(self):
        self.close()
