"""
This module will create an amber mdin file for either sander or
pmemd (or others). The program specification loads the appropriate
dictionaries with default values, etc. It can read and write mdins.
"""
from io import TextIOBase
from re import findall

# This module will create and read a sander/pmemd input
from .cntrl import cntrl
from .ewald import ewald
from .pb import pb
from .qmmm import qmmm
from .rism import rism
from .gbnsr6 import gbnsr6
from ...exceptions import InputError


def addOn(line, string, file):
    if len(line.strip()) == 0:
        return line + string
    elif len(line) + len(string) > 40:
        file.write(line + '\n')
        return ' ' + string
    else:
        return line + string

def splitStr(string):
    # split string at commas not within single or double quotes
    pattern = r"""('[^']*'|"[^"]*"|[^'",\s]+)"""
    matches = findall(pattern, string)

    return matches


class Mdin:

    def __init__(self, program = 'sander', verbosity = 1):
        # define instance data
        self.program = program   # which program we're creating the input file for
        self.cntrl_obj = cntrl() # object with cntrl namelist vars in a dictionary
        self.ewald_obj = ewald() # object with ewald namelist vars in a dictionary
        self.pb_obj = pb()       # object with pb namelist vars in a dictionary
        self.qmmm_obj = qmmm()   # object with qmmm namelist vars in a dictionary
        self.rism_obj = rism()   # object with rism namelist vars in a dictionary
        self.gbnsr6_obj = gbnsr6()   # object with gbnsr6 namelist vars in a dictionary
        self.verbosity = 0       # verbosity level: 0 -- print nothing
                                 #                  1 -- print errors
                                 #                  2 -- 1 + warnings
                                 #                  3 -- 2 + notes
        self.cards = []          # array that has all of the input cards that come
                                 # after namelists
        self.cntrl_nml = {}      # dictionary with cntrl namelist vars
        self.cntrl_nml_defaults = {} # dictionary with default cntrl namelist vars
        self.ewald_nml = {}          # dictionary with ewald namelist vars
        self.ewald_nml_defaults = {} # dictionary with default ewald namelist vars
        self.pb_nml = {}             # dictionary with pb namelist vars
        self.pb_nml_defaults = {}    # dictionary with default pb namelist vars
        self.qmmm_nml = {}           # dictionary with qmmm namelist vars
        self.qmmm_nml_defaults = {}  # dictionary with default qmmm namelist vars
        self.rism_nml = {}           # dictionary with rism namelist vars
        self.rism_nml_defaults = {}  # dictionary with default rism namelist vars
        self.gbnsr6_nml = {}           # dictionary with gbnsr6 namelist vars
        self.gbnsr6_nml_defaults = {}  # dictionary with default gbnsr6 namelist vars
        self.valid_namelists = {}    # array with valid namelists for each program
        self.title = 'mdin prepared by ParmEd'   # title for the mdin file


        if self.program == "sander":
            self.cntrl_nml = self.cntrl_obj.sander
            self.ewald_nml = self.ewald_obj.sander
            self.pb_nml = self.pb_obj.sander
            self.qmmm_nml = self.qmmm_obj.sander
            self.rism_nml = self.rism_obj.sander
            self.valid_namelists = {'cntrl', 'ewald', 'qmmm', 'pb', 'rism'}
        elif self.program == "sander.APBS":
            self.cntrl_nml = self.cntrl_obj.sander
            self.pb_nml = self.pb_obj.sanderAPBS
            self.valid_namelists = {'cntrl', 'apbs'}
        elif self.program == "pmemd":
            self.cntrl_nml = self.cntrl_obj.pmemd
            self.ewald_nml = self.ewald_obj.pmemd
            self.valid_namelists = {'cntrl', 'ewald'}
        elif self.program == "gbnsr6":
            self.cntrl_nml = self.cntrl_obj.gbnsr6
            self.gbnsr6_nml = self.gbnsr6_obj.gbnsr6
            self.valid_namelists = ['cntrl', 'gb']
        else:
            raise InputError(f"Unrecognized program [{self.program}]")

        self.cntrl_nml_defaults = self.cntrl_nml.copy()
        self.ewald_nml_defaults = self.ewald_nml.copy()
        self.pb_nml_defaults = self.pb_nml.copy()
        self.qmmm_nml_defaults = self.qmmm_nml.copy()
        self.rism_nml_defaults = self.rism_nml.copy()
        self.gbnsr6_nml_defaults = self.gbnsr6_nml.copy()

    def write(self, filename: str = 'mdin'):

        def write_nml(nml, defaults, file: TextIOBase, header, print_header=False) -> None:
            # automatic indent of single space
            line = ' '
            header_printed = False
            # add any variable that is different from the default to the mdin file
            for var, val in nml.items():
                if val == defaults[var]:
                    continue
                if not header_printed:
                    file.write(f"{header}\n")
                    header_printed = True
                if isinstance(val, (list, tuple)):
                    if isinstance(val[0], str):
                        temp_out = ','.join(map(lambda x: f"'{x}'", val))
                        line = addOn(line, f"{var}={temp_out}, ", file)
                    else:
                        line = addOn(line, f"{var}={','.join(map(str, val))}, ", file)
                elif isinstance(val, str):
                    line = addOn(line, f"{var}='{val}', ", file)
                else:
                    line = addOn(line, f"{var}={val}, ", file)
            # flush any remaining items that haven't yet been printed to the mdin file
            if line.strip():
                file.write(f"{line}\n")

            # gbnsr6 requieres always print the cntrl header, so we check if was already printed
            if not header_printed and print_header:
                file.write(f"{header}\n")
                header_printed = True

            # End namelist
            if header_printed:
                file.write('/\n')

        # open the file for writing and write the header and &cntrl namelist
        file = open(filename,'w')
        file.write(self.title + '\n')
        write_nml(self.cntrl_nml, self.cntrl_nml_defaults, file, "&cntrl", True)
        write_nml(self.ewald_nml, self.ewald_nml_defaults, file, "&ewald")
        pb_nml_name = "&apbs" if self.program == "sander.APBS" else "&pb"
        write_nml(self.pb_nml, self.pb_nml_defaults, file, pb_nml_name)
        write_nml(self.gbnsr6_nml, self.gbnsr6_nml_defaults, file, "&gb")
        if self.cntrl_nml.get('irism'):
            write_nml(self.rism_nml, self.rism_nml_defaults, file, "&rism", True)
        if self.cntrl_nml.get("ifqnt") == 1:
            write_nml(self.qmmm_nml, self.qmmm_nml_defaults, file, "&qmmm")

        # Write the cards to the input file
        file.write("\n".join([card.strip() for card in self.cards]))
        if len(self.cards) != 0:
            file.write('\nEND\n')

        file.close()

    def read(self, filename = 'mdin'):
        lines = open(filename, 'r').readlines()

        # split up input file into separate fields by comma
        blocks = []  # namelists in the order they appear
        block_fields = [] # array of arrays that correspond to entries in
                          # namelists found in "blocks" above
        inblock = False
        lead_comment = True
        for i in range(len(lines)):
            if not inblock and not lines[i].strip().startswith('&') and lead_comment:
                continue
            elif not inblock and not lines[i].strip().startswith('&') and not lead_comment:
                final_ended = True
                for j in range(i,len(lines)):
                    if lines[j].strip().startswith('&'):
                        final_ended = False
                if final_ended and len(lines[i].strip()) != 0:
                    self.cards.append(lines[i])
            elif not inblock and lines[i].strip().startswith('&'):
                lead_comment = False
                inblock = True
                block = lines[i].strip()[1:].lower()
                blocks.append(block)    # add the name of the namelist to "blocks"
                block_fields.append([]) # add empty array to be filled with entries
                                        # for given namelist
                if not block in self.valid_namelists:
                    raise InputError(f'Invalid namelist ({lines[i].strip()}) in input file ({filename}) for {self.program}')
            elif inblock and (lines[i].strip() == '/' or lines[i].strip() == '&end'):
                inblock = False
            elif inblock and lines[i].strip().startswith('&'):
                raise InputError(f"Invalid input file ({filename}). Namelist not terminated")
            elif inblock:
                items = splitStr(lines[i].strip())
                j = 0
                while j < len(items):
                    items[j] = items[j].strip()
                    if len(items[j]) == 0:
                        items.pop(j)
                    else:
                        j += 1
                block_fields[len(block_fields)-1].extend(items)

        # take out the last END in the cards if it's there
        if len(self.cards) != 0 and self.cards[len(self.cards)-1].strip().upper() == 'END':
            self.cards.pop()

        # combine any multi-element fields: e.g. rstwt=1,2,3,
        begin_field = -1
        for i in range(len(block_fields)):
            for j in range(len(block_fields[i])):
                if block_fields[i][j].startswith("'") or block_fields[i][j].startswith("\"") or not '=' in block_fields[i][j]:
                    if begin_field == -1:
                        raise InputError(f'Invalid input file ({filename}).')
                    else:
                        if block_fields[i][j].startswith("'") or block_fields[i][j].startswith("\""):
                            block_fields[i][begin_field] += block_fields[i][j]
                        
                            # prevents treating '=' in quoted field values (e.g. wildcards for restraintmask)
                            # from delimiting new fields
                            if '=' in block_fields[i][j]:
                                block_fields[i][j] = block_fields[i][j].replace('=','')
                        
                        else:
                            block_fields[i][begin_field] += ',' + block_fields[i][j]
                else:
                    begin_field = j

        # now parse through the options and add them to the dictionaries
        for i in range(len(block_fields)):
            for j in range(len(block_fields[i])):
                if not '=' in block_fields[i][j]:
                    continue
                else:
                    var = block_fields[i][j].split('=', 1)
                    self.change(blocks[i], var[0].strip(), var[1].strip())

    def change(self, namelist, variable, value):
        """ Change the value of a variable without adding a new key-pair """

        namelist_map = dict(
            cntrl=(self.cntrl_nml, self.cntrl_nml_defaults),
            ewald=(self.ewald_nml, self.ewald_nml_defaults),
            pb=(self.pb_nml, self.pb_nml_defaults),
            apbs=(self.pb_nml, self.pb_nml_defaults),
            qmmm=(self.qmmm_nml, self.qmmm_nml_defaults),
            rism=(self.rism_nml, self.rism_nml_defaults),
            gb=(self.gbnsr6_nml, self.gbnsr6_nml_defaults)
        )

        variable = variable.lower()
        if isinstance(value, str) and len(value) > 0:
            if (value[0] == value[-1] == "'") or (value[0] == value[-1] == '"'):
                value = value[1:-1]

        if namelist in namelist_map:

            nml, nml_defaults = namelist_map[namelist]
            if variable in nml:
                mytype = type(nml_defaults[variable])
                nml[variable] = mytype(value)
            else:
                raise InputError(f"Unknown variable {variable} in &{namelist}")
        else:
            raise InputError(f'Unknown namelist ({namelist})!')

    def check(self):
        return True

    def SHAKE(self):
        self.change('cntrl','ntf', 2)
        self.change('cntrl','ntc', 2)
        self.change('cntrl','dt', 0.002)

    def constPressure(self, press=1.0, taup=1.0):
        self.change('cntrl','ntb', 2)
        self.change('cntrl','ntp', 1)
        self.change('cntrl','pres0', press)
        self.change('cntrl','taup', taup)

    def constVolume(self):
        self.change('cntrl','ntb', 1)
        self.change('cntrl','ntp', 0)

    def constTemp(self, ntt=3, temp=300.0, gamma_ln=2.0, ig=-1, tautp=1.0):
        self.change('cntrl','ntt', ntt)
        self.change('cntrl','temp0', temp)
        self.change('cntrl','tempi', temp)
        self.change('cntrl','gamma_ln', gamma_ln if ntt==3 else 0)
        self.change('cntrl','ig', ig)
        self.change('cntrl','tautp', tautp)

    def constpH(self, solvph=7.0, igb=2, ntcnstph=10):
        self.change('cntrl','icnstph', 1)
        self.change('cntrl','solvph', solvph)
        self.change('cntrl','ntcnstph', ntcnstph)
        self.change('cntrl','igb', igb)
        self.change('cntrl','ntb', 0)
        self.change('cntrl','saltcon', 0.1)

    def restrainHeavyAtoms(self, restraint_wt=0.0):
        self.change('cntrl','ntr', 1)
        self.change('cntrl','restraint_wt', restraint_wt)
        self.change('cntrl','restraintmask', '!@H=')

    def restrainBackbone(self, restraint_wt=0.0):
        self.change('cntrl','ntr', 1)
        self.change('cntrl','restraint_wt', restraint_wt)
        self.change('cntrl','restraintmask', '@N,CA,C')

    def genBorn(self, igb=5, rgbmax=25.0):
        self.change('cntrl','igb', igb)
        self.change('cntrl','ntb', 0)
        self.change('cntrl','ntp', 0)
        self.change('cntrl','rgbmax', rgbmax)

    def time(self, time=1000.0, dt=None): # time in ps
        if dt is None:
            if self.cntrl_nml['ntc'] == 2 and self.cntrl_nml['ntf'] == 2:
                dt = 0.002
            else:
                dt = 0.001
        time = int(time / dt)

        self.change('cntrl','dt', dt)
        self.change('cntrl','nstlim', time)
        self.change('cntrl','imin', 0)

    def heat(self, tempi=0.0, temp0=300.0, ntt=3, tautp=5.0, ig=-1, gamma_ln=5.0):
        self.constVolume()
        self.change('cntrl','tempi', tempi)
        self.change('cntrl','temp0', temp0)
        self.change('cntrl','ntt', ntt)
        self.change('cntrl','tautp', tautp)
        self.change('cntrl','ig', ig)
        self.change('cntrl','gamma_ln', gamma_ln if ntt==3 else 0)

    def restart(self,ntx=5):
        self.change('cntrl','irest',1)
        self.change('cntrl','ntx',ntx)

    def TI(self, clambda=0.0):
        self.change('cntrl','clambda', clambda)
        self.change('cntrl','icfe',1)

    def softcore_TI(self, scalpha=0.5, scmask='', crgmask='', logdvdl=0):
        self.change('cntrl','icfe',1)
        self.change('cntrl','ifsc',1)
        self.change('cntrl','scalpha',scalpha)
        self.change('cntrl','scmask',scmask)
        self.change('cntrl','crgmask',crgmask)
        self.change('cntrl','logdvdl',logdvdl)

    def minimization(self, imin=1, maxcyc=1, ncyc=10, ntmin=1):
        self.change('cntrl','imin', imin)
        self.change('cntrl','maxcyc', maxcyc)
        self.change('cntrl','ncyc', ncyc)
        self.change('cntrl','ntmin', ntmin)

    def rism(self, imin=5, irism=1):
        self.change('cntrl', 'imin', imin)
        self.change('cntrl', 'irism', irism)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    def AddCard(self, title='Residues in card', cardString='RES 1'):
        self.cards.append('%s\n%s\nEND' % (title, cardString))
