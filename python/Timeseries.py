
import AtomGroup

class TimeseriesCollection(object):
    def __init__(self):
        self.timeseries = []
    def __len__(self):
        return len(self.timeseries)
    def __getitem__(self, item):
        if (type(item) is int) or (type(item) is slice):
            return self.timeseries[item]
        else: raise IndexError
    def addTimeseries(self, ts):
        # if ts is not Timeseries do something
        self.timeseries.append(ts)
    def getAtomList(self):
        self._atomlist = []
        for ts in self.timeseries:
            self._atomlist += ts.getAtomList()
        return self._atomlist
    def getFormat(self):
        return ''.join([ts.getFormatCode() for ts in self.timeseries])
    def getBounds(self):
        if not hasattr(self, '_atomlist'):
            self.getAtomList()
        if not len(self._atomlist) == 0:
            lowerb = min(self._atomlist)
            upperb = max(self._atomlist)
        else: lowerb = upperb = 0
        return lowerb, upperb
    def getDataSize(self):
        return sum([ts.getDataSize() for ts in self.timeseries])
    def getAtomCounts(self):
        atomCounts = []
        for ts in self.timeseries:
            atomCounts += ts.getAtomCounts()
        return atomCounts
    def getAuxData(self):
        auxData = []
        for ts in self.timeseries:
            auxData += ts.getAuxData()
        return auxData
    def clear(self):
        self.timeseries = []
    def compute(self, trj, start=0, stop=-1, skip=1):
        if (start < 0): start += len(trj)
        if (stop < 0): stop += len(trj)
        if (stop <= start): raise Exception("Stop frame is lower than start frame")
        if ((start < 0) or (start >= len(trj)) or
           (stop < 0) or (stop >= len(trj))):
               raise IndexError
        atomlist = self.getAtomList()
        format = self.getFormat()
        lowerb, upperb = self.getBounds()
        sizedata = self.getDataSize()
        atomcounts = self.getAtomCounts()
        auxdata = self.getAuxData()
        self.data = trj._read_timecorrel(atomlist, atomcounts, format, auxdata, sizedata, lowerb, upperb, start, stop, skip)
        # Now remap the timeseries data to each timeseries
        #typestr = "|f8"
        start = 0
        for t in self.timeseries:
            finish = t.getDataSize()
            datasize = len(t.getFormatCode())
            subarray = self.data[start:start+finish]
            if datasize != 1: subarray.shape = (datasize, subarray.shape[0]/datasize, -1)
            t.__data__ = subarray
            # XXX The __array_interface__ only works with numpy
            #t.__array_interface__ = {"data": (, False), "shape": subarray.shape, "typestr": typestr,"version": 3, "stride": None}
            start += finish


class Timeseries(object):
    def __init__(self, code, atoms, dsize):
        if isinstance(atoms, AtomGroup.AtomGroup):
            self.atoms = atoms._atoms
        elif isinstance(atoms, list):
            self.atoms = atoms
        elif isinstance(atoms, AtomGroup.Atom):
            self.atoms = [atoms]
        else: raise Exception("Invalid atoms passed to %s timeseries"%self.__class__.__name__)
        self.code = code
        self.numatoms = len(self.atoms)
        self.dsize = dsize
    def __array__(self):
        return self.__data__
    def __getitem__(self, item):
        if (type(item) is int) or (type(item) is slice):
            return self.__data__[item]
        else: raise IndexError
    def __getattr__(self, attr):
        if attr == "shape":
            return self.__data__.shape
        else: return super(Timeseries, self).__getattribute__(attr)
    def getAtomList(self):
        return [a.number for a in self.atoms]
    def getFormatCode(self):
        return self.code
    def getDataSize(self):
        return self.dsize
    def getNumAtoms(self):
        return self.numatoms
    def getAtomCounts(self):
        return [self.numatoms]
    def getAuxData(self):
        return [0.]*self.numatoms

class Atom(Timeseries):
    def __init__(self, code, atoms):
        if code not in ('x', 'y', 'z', 'v'):
            raise Exception("Bad code")
        if code == 'v': size = 3
        else: size = 1
        if isinstance(atoms, AtomGroup.AtomGroup):
            numatoms = len(atoms._atoms)
        elif isinstance(atoms, list):
            numatoms = len(atoms)
        elif isinstance(atoms, AtomGroup.Atom):
            numatoms = 1
        else: raise Exception("Invalid atoms passed to %s timeseries"%self.__class__.__name__)
        Timeseries.__init__(self, code*numatoms, atoms, size*numatoms)
    def getAtomCounts(self):
        return [1,]*self.numatoms

class Angle(Timeseries):
    def __init__(self, atoms):
        if not len(atoms) == 3:
            raise Exception("Angle timeseries requires a 3 atom selection")
        Timeseries.__init__(self, 'a', atoms, 1)

class Dihedral(Timeseries):
    def __init__(self, atoms):
        if not len(atoms) == 4:
            raise Exception("Dihedral timeseries requires a 4 atom selection")
        Timeseries.__init__(self, 'h', atoms, 1)

class Distance(Timeseries):
    def __init__(self, code, atoms):
        if code not in ('d', 'r'):
            raise Exception("Bad code")
        if code == 'd': size = 3
        else: size = 1
        if not len(atoms) == 2:
            raise Exception("Distance timeseries requires a 2 atom selection")
        Timeseries.__init__(self, code, atoms, size)

class CenterOfGeometry(Timeseries):
    def __init__(self, atoms):
        Timeseries.__init__(self, 'm', atoms, 3)
    def getAuxData(self):
        return [1.]*self.numatoms

class CenterOfMass(Timeseries):
    def __init__(self, atoms):
        Timeseries.__init__(self, 'm', atoms, 3)
    def getAuxData(self):
        return [a.mass for a in self.atoms]
