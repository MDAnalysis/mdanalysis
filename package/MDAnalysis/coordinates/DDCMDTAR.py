from __future__ import absolute_import

from six import StringIO
from six import BytesIO
import tarfile
import re
import struct

import logging
import collections
import numpy as np

import itertools

from ..lib import util
from . import base

from .DDC import DDC

logger = logging.getLogger("MDAnalysis.coordinates.DDCMDTAR")


class DDCMDTARReader(base.ReaderBase):
    """
    DDCMD TAR Format:

    """

    format = 'DDCMDTAR'
    units = {'time': 'ps', 'length': 'Angstrom'}


    def __init__(self, filename, **kwargs):
        """Read coordinates from *filename*.

        *filename* can be a gzipped or bzip2ed compressed PDB file.

        If the pdb file contains multiple MODEL records then it is
        read as a trajectory where the MODEL numbers correspond to
        frame numbers.
        """
        super(DDCMDTARReader, self).__init__(filename, **kwargs)

        try:
            self.n_atoms = kwargs['n_atoms']
        except KeyError:
            # hackish, but should work and keeps things DRY
            # regular MDA usage via Universe doesn't follow this route
            with tarfile.open(self.filename, "r|") as tf:
                with open(self.filename, "rb") as fdata:
                    for snapshot in tf:
                        fdata.seek(snapshot.offset_data)
                        data = fdata.read(snapshot.size)
                        self.sf = StringIO(data)
                        self.headerDict = DDC.parseHeader()
                        self.n_atoms = self.headerDict['nrecord']
                        break

        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        self.taridx=None
        if 'taridx' in kwargs:
            self.taridx=kwargs['taridx']

        self.onesnapshot = False
        if 'snapshot' in kwargs:
            self.onesnapshot=True
            self.onesnapshotname=kwargs['snapshot']
            self.n_frames = 1
        else:
            #self.n_frames = len(self.tarframes)
            self.n_frames = self._get_nframes()
            #self.n_frames = self._get_n_frames()

        #self.tarframes, self.tarfilenames = self._get_tar_frames()
        self.targenerator = self._get_tar_frames()

        print('A total of '+str(self.n_frames)+' frames are found')
        #self._read_frame(0)

    #def _get_snapshot(self, snapshot):
    #    for targen in self.targenerator:
    #        if targen[1] == snapshot:
    #            yield targen

    #def _get_n_frames(self):
    #    with tarfile.open(self.filename, "r|") as tf:
    #        count = sum(1 for member in tf if member.isreg())
    #    return count

    def _get_nframes(self):
        if self.onesnapshot:
            return 1

        count=0
        with tarfile.open(self.filename, "r|") as tf:
            for frame in tf:
                count=count+1
        return count

    def _get_tar_frames(self):
        #tar_data_list=[]
        #tar_filename_list=[]
        if self.taridx == None:
            # if no taridx file is provided
            with tarfile.open(self.filename, "r|") as tf:
                with open(self.filename, "rb") as fdata:
                    if self.onesnapshot:
                        count=0
                        for frame in tf:
                            if count == 1:
                                break;
                            if frame.name == self.onesnapshotname:
                                fdata.seek(frame.offset_data)
                                data = fdata.read(frame.size)
                                count=count+1
                                yield data
                    else:
                        for frame in tf:
                            #print(frame.name + " is " + str(frame.size) + " bytes in size")
                            fdata.seek(frame.offset_data)
                            data = fdata.read(frame.size)

                            yield data
                            #tar_data_list.append(data)
                            #tar_filename_list.append(frame.name)
        else:
            # Read the list file, to efficiently get the set of each
            members = {}
            m = 0
            with open(self.taridx, "r") as idx:
                for line in idx:
                    items = line.rstrip("\r\n").split(',')
                    name = items[0]
                    pos = int(items[1])
                    sz = int(items[2])
                    if sz > 0:  # This (name) is not a directory...
                        members[name] = (pos, sz)
                        m = m + 1

            with open(self.filename, "rb") as tf:
                if self.onesnapshot:
                    if self.onesnapshotname in members:
                        (pos, sz) = members[self.onesnapshotname]
                        tf.seek(pos)  # Seek to member data
                        data = tf.read(int(sz))  # Read member data
                        yield data
                    else:
                        raise Exception("Snapshot "+ self.onesnapshotname + "does not exist")

                for name in sorted(members.keys()):
                    (pos, sz) = members[name]

                    tf.seek(pos)  # Seek to member data
                    data = tf.read(int(sz))  # Read member data
                    yield data


        #return tar_data_list, tar_filename_list


    def _reopen(self):
        # Pretend the current TS is -1 (in 0 based) so "next" is the
        # 0th frame
        #self.tarframes = self._get_tar_frames()
        self.ts.frame = -1
        self.targenerator = self._get_tar_frames()

    def _read_next_timestep(self, ts=None):
        if ts is None:
            ts = self.ts
        else:
            # TODO: cleanup _read_frame() to use a "free" Timestep
            raise NotImplementedError("PDBReader cannot assign to a timestep")
        # frame is 1-based. Normally would add 1 to frame before calling
        # self._read_frame to retrieve the subsequent ts. But self._read_frame
        # assumes it is being passed a 0-based frame, and adjusts.
        frame = self.frame + 1
        targen = next(self.targenerator)
        return self._read_frame_tdata(frame, targen)


    def _read_frame(self, frame):
        try:
            self.targenerator = self._get_tar_frames()
            tar_data = next(itertools.islice(self.targenerator, frame, frame+1))

            #tar_data = self.tarframes[frame]
        except IndexError:  # out of range of known frames
            raise IOError

        return self._read_frame_tdata(frame, tar_data)


    def _read_frame_tdata(self, frame, tar_data):

        self.sf = BytesIO(tar_data)

        headerDict = DDC.parseHeader(self.sf)

        ddcmd_dict = None
        if headerDict['datatype'] == 'FIXRECORDBINARY':
            ddcmd_dict = self.read_BinaryDDCMD(headerDict)
        else:
            "The datatype should be FIXRECORDBINARY"

        ddcmd_od = collections.OrderedDict(sorted(ddcmd_dict.items()))

        coords = []

        for k, v in ddcmd_od.items():
            coords.append(v)

        if len(coords) != self.n_atoms:
            raise ValueError("Read an incorrect number of atoms\n"
                             "Expected {expected} got {actual}"
                             "".format(expected=self.n_atoms, actual=len(coords)))

        self.ts = self._Timestep.from_coordinates(
            coords,
            **self._ts_kwargs)

        unitcell = np.zeros(6, dtype=np.float32)
        unitcell[:] = headerDict['lx'], headerDict['ly'], headerDict['lz'], 90, 90, 90
        self.ts.time = headerDict['time']
        self.ts._unitcell[:] = unitcell
        self.ts.frame = 0  # 0-based frame number
        if self.convert_units:
            # in-place !
            self.convert_pos_from_native(self.ts._pos)
            self.convert_pos_from_native(self.ts._unitcell[:3])

        self.ts.frame = frame

        return self.ts

    def close(self):
        #self._pdbfile.close()
        pass

################################################################################
# ddcMD parser
################################################################################

    def read_BinaryDDCMD(self, headerDict):
        # index for fields
        #requiredFields = ['id', 'species', 'rx', 'ry', 'rz']
        requiredFields = ['id', 'rx', 'ry', 'rz']
        for field in requiredFields:
            if field not in headerDict['fieldName']:
                print ("Missing field name: ", field)
                exit(-1)

        fieldBits, fieldFmts=DDC._getBitFormat(headerDict)
        fieldNames = headerDict['fieldName']

        atomCount = 0

        fieldData={}
        ddcmd_dict={}

        offset=headerDict['offset']
        self.sf.seek(0)
        self.sf.seek(offset)

        readSuccess=True
        while readSuccess:
            # loop through the field
            for idx, bit in enumerate(fieldBits):
                data = self.sf.read(bit)
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


