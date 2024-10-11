import re
import sys

class DDC:

    # check if gid is a hex
    @staticmethod
    def _isHex(headerDict):

        # the hex print out is assume only indicated in the fieldFormat
        if 'fieldFormat' not in headerDict:
            return False

        fields=headerDict['fieldName']
        idx = fields.index("id")
        fieldFormats=headerDict['fieldFormat']
        gidFmt=fieldFormats[idx]

        if gidFmt[-1]=='x':
            return True

        return False

    @staticmethod
    def _getBitFormat(headerDict):

        # figure out bits for each field
        fieldTypes = headerDict['fieldTypes']
        fieldFmts=[]
        fieldBits=[]
        for type in fieldTypes:
            bit=int(type[1])
            fieldBits.append(bit)
            if type[0] == 'u' and bit == 4:
                fieldFmts.append('<L')
            elif type[0] == 'u' and bit == 8:
                fieldFmts.append('<Q')
            elif type[0] == 'f' and bit == 4:
                fieldFmts.append('<f')
            elif type[0] == 'f' and bit == 8:
                fieldFmts.append('<d')
            else:
                raise Exception("Add format support for "+type)

        if len(fieldTypes) != len(fieldFmts):
            raise Exception("Number of field types doesn't equal to that of field formats")

        if len(fieldTypes) != len(fieldBits):
            raise Exception("Number of field types doesn't equal to that of field bits")

        return fieldBits, fieldFmts

    @staticmethod
    def _checkUnit(headerDict):
        """
        ddcMD has Ang Angstrom, Bohr, a0, meter, mm, um, nm units
        MDanalysis supports Angstrom, nm, pm, fm
        Exception will raise if unit is not Angstrom or nm
        :param headerDict:
        :return:
        """
        fieldNames=headerDict['fieldName']
        rxIdx=fieldNames.index('rx')
        fieldUnits=headerDict['fieldUnits']
        rxUnit=fieldUnits[rxIdx]
        if rxUnit=='Ang' or rxUnit=='Angstrom':
            # use default unit
            rxUnit = 'Angstrom'
        elif rxUnit=='nm':
            pass
        else:
            raise Exception("rx unit is not supported "+ rxUnit)

        return rxUnit

#    @staticmethod
#    def _getHeaderOffset(filehandle):
#        # find out the offset for the first file
#        filehandle.seek(0)
#        offset=0
#        endFlg=False
#
#        for line in filehandle:
#            if endFlg and not line.isspace():
#                break
#            if line[0] == '}':
#                endFlg=True
#
#            offset = offset + len(line)
#
#        return offset

    @staticmethod
    def _parseObj(object):
        replaceObj=object.replace("\n", " ")
        splitObjs=re.split('{|}', replaceObj)

        if len(splitObjs)<3:
            print("Object doesn't enclose with a pair of {}")

        contentList=splitObjs[1].split(";")

        contentDict={}

        for content in contentList:
            pair=content.split("=")
            if len(pair) == 2:
                contentDict[pair[0].strip()]=pair[1].strip()
            #else:
            #   print ("Wrong pair in object:", content)

        return contentDict

    @staticmethod
    def parseHeader(filehandle):

        header = ""
        offset = 0
        endFlg=False

        #try:
        #    for line in filehandle:
        #        if endFlg and not line.isspace():
        #            break
        #        if line[0] == '}':
        #            endFlg = True
        #        header = header + line
        #        offset = offset + len(line)
        #except:
        #    print("read bytes error")

        while True:
            line = filehandle.readline()
            if sys.version_info[0] >= 3:
                try:
                    line=str(line, 'utf-8')
                except:
                    break

            if endFlg and not line.isspace():
                break
            if line[0] == '}':
                endFlg = True

            header = header + line
            offset = offset + len(line)

            if not line:
                break

        headerObj=DDC._parseObj(header)

        requiredKeys=['datatype', 'h', 'field_names', 'field_types', 'nfields', 'nfiles', 'nrecord', 'time']
        for key in requiredKeys:
            if key not in headerObj:
                print ("Missing required key in header :", key)
                exit(-1)


        headerDict={}
        headerDict['offset'] = offset
        headerDict['datatype'] = headerObj['datatype']

        boxSplit=headerObj['h'].split()
        if len(boxSplit)<9:
            print ("Item of box values should be 9 but it is ", len(boxSplit))
            exit(-1)
        headerDict['lx'] = float(boxSplit[0])
        headerDict['ly'] = float(boxSplit[4])
        headerDict['lz'] = float(boxSplit[8])

        headerDict['nfield'] = int(headerObj['nfields'])
        headerDict['nfiles'] = int(headerObj['nfiles'])
        headerDict['nrecord'] = int(headerObj['nrecord'])
        headerDict['fieldName'] = headerObj['field_names'].split()
        headerDict['fieldTypes'] = headerObj['field_types'].split()
        headerDict['time'] = float(headerObj['time'].split()[0])*0.001

        if 'species' in headerObj:
            headerDict['species'] = headerObj['species'].split()

        if 'field_format' in headerObj:
            headerDict['fieldFormat'] = headerObj['field_format'].split()

        if 'field_units' in headerObj:
            headerDict['fieldUnits'] = headerObj['field_units'].split()

        return headerDict

