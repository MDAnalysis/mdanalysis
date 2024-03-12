from __future__ import absolute_import
import numpy as np

import struct
import collections

from ..lib import util
from . import base
from .DDC import DDC

class DDCMDReader(base.SingleFrameReaderBase):
    """
    DDCMD Format:

    """

    format = 'DDCMD'
    units = {'time': None, 'length': 'Angstrom'}

    def read_BinaryDDCMD(self, headerDict):
        # index for fields
        requiredFields = ['id', 'rx', 'ry', 'rz']
        for field in requiredFields:
            if field not in headerDict['fieldName']:
                print ("Missing field name: ", field)
                exit(-1)

        fileNames = []
        fileBaseName = self.filename.split("#")[0]
        for i in range(headerDict['nfiles']):
            numbering = str(i).zfill(6)
            newfileName = fileBaseName + "#" + numbering
            fileNames.append(newfileName)

        fieldBits, fieldFmts=DDC._getBitFormat(headerDict)
        fieldNames = headerDict['fieldName']

        atomCount = 0

        fieldData={}
        ddcmd_dict={}

        for idx, filename in enumerate(fileNames):
            with open(filename, "rb") as f:
                if idx == 0:  # the first file has header
                    offset=headerDict['offset']
                    f.seek(offset)
                    #print("offset=", offset)

                readSuccess=True
                while readSuccess:
                    # loop through the field
                    for idx, bit in enumerate(fieldBits):
                        data = f.read(bit)
                        if data:
                            fmt=fieldFmts[idx]
                            fieldData[fieldNames[idx]] = struct.unpack(fmt, data)[0]
                        else:
                            readSuccess = False
                            break
                    if readSuccess:
                        # get data from dict
                        gid = fieldData['id']
                        rx = fieldData['rx']
                        ry = fieldData['ry']
                        rz = fieldData['rz']
                        ddcmd_dict[gid] = [rx, ry, rz]
                        atomCount = atomCount + 1

        nrecord=headerDict['nrecord']
        if atomCount != nrecord:
            print("Warning number of atoms read {} is not equal to ddcMD nrecord {}".format(atomCount, nrecord))

        return ddcmd_dict


    def read_DDCMD(self, headerDict):
        # index for fields
        requiredFields = ['id', 'rx', 'ry', 'rz']
        for field in requiredFields:
            if field not in headerDict['fieldName']:
                print ("Missing field name: ", field)
                exit(-1)

        fields=headerDict['fieldName']
        gidIndex = fields.index("id")
        rxIndex = fields.index("rx")
        ryIndex = fields.index("ry")
        rzIndex = fields.index("rz")
        fieldSize=headerDict['nfield']

        hex=DDC._isHex(headerDict)

        fileNames = []
        fileBaseName = self.filename.split("#")[0]
        for i in range(headerDict['nfiles']):
            numbering = str(i).zfill(6)
            newfileName = fileBaseName + "#" + numbering
            fileNames.append(newfileName)
        print(fileNames)
        atomCount = 0
        ddcmd_dict = {}
        for idx, filename in enumerate(fileNames):
            with open(filename, "r") as f:
                if idx == 0:  # the first file has header
                    offset=headerDict['offset']
                    f.seek(offset)
                    #print("offset=", offset)

                for line in f:
                    strs = line.split()
                    if len(strs) >= fieldSize:
                        if hex:
                            gid = int(strs[gidIndex], 16)
                        else:
                            gid = int(strs[gidIndex])
                        rx=float(strs[rxIndex])
                        ry=float(strs[ryIndex])
                        rz=float(strs[rzIndex])
                        ddcmd_dict[gid] = [rx, ry, rz]
                        atomCount = atomCount + 1
        print("Atom count=", atomCount)
        nrecord=headerDict['nrecord']
        if atomCount != nrecord:
            print("Warning number of atoms read {} is not equal to ddcMD nrecord {}".format(atomCount, nrecord))

        return ddcmd_dict

    def _read_first_frame(self):

        with open(self.filename, "rb") as f:
            headerDict=DDC.parseHeader(f)

        # If fieldUnits is provided, double-check if rx unit is angstom or nm
        # Otherwise use default angstom
        if 'fieldUnits' in headerDict:
            DDCMDReader.units['length']=DDC._checkUnit(headerDict)

        ddcmd_dict=None
        if headerDict['datatype'] == 'FIXRECORDBINARY':
            ddcmd_dict = self.read_BinaryDDCMD(headerDict)
        else:
            ddcmd_dict = self.read_DDCMD(headerDict)

        ddcmd_od = collections.OrderedDict(sorted(ddcmd_dict.items()))

        coords = []

        for k, v in ddcmd_od.items():
            coords.append(v)


        self.n_atoms = len(coords)
        self.ts = self._Timestep.from_coordinates(
            coords,
            **self._ts_kwargs)

        unitcell = np.zeros(6, dtype=np.float32)
        unitcell[:] = headerDict['lx'], headerDict['ly'], headerDict['lz'], 90, 90, 90
        self.ts._unitcell[:] = unitcell
        self.ts.frame = 0  # 0-based frame number
        if self.convert_units:
            # in-place !
            self.convert_pos_from_native(self.ts._pos)
            self.convert_pos_from_native(self.ts._unitcell[:3])