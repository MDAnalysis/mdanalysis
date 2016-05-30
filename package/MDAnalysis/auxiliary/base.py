import numpy as np

class AuxReader(object):
    """ Base class for auxiliary readers 
    To read in from file or e.g. array
    """
      
    # TODO deal with changing/different units

    def __init__(self, auxdata, auxnames, repres_ts_as='closest', **kwargs):
        # TODO default auxnames...
        # TODO allow to name each column + pick which to keep
        # TODO pass in dt, initial_time for when time not in data
        # TODO cutoff for representive closest 
        self.names = auxnames
        self.repres_ts_as = repres_ts_as

        self.step = 0

    def __iter__(self):
        self.reopen() # ...if array...?
        return self

    def next(self):
        """ Move to next step of data """
        return self.read_next_step()

    def __next__(self):
        """ Move to next step of data """
        return self.next()

    # def read_next_step():  
        # DEFINE FOR EACH

    # TODO - add __enter__ and __exit__ methods



class AuxFileReader(AuxReader):
    """ Base class for auxiliary readers dealing with files """
    def __init__(self, filename, auxname, **kwargs):
        super(AuxFileReader, self).__init__(filename, auxname, **kwargs)

        self.auxfilename = filename
        self.auxfile = open(filename)

    def rewind(self):
        """ position back at start and read first step """
        self.reopen()
        self.read_next_step()

    def close(self):
        """ close file if open """
        if self.auxfile == None:
            return
        self.auxfile.close()
        self.auxfile = None

    def reopen(self):
        """ reposition to just before first step """
        self.auxfile.close()  
        self.auxfile = open(self.auxfilename) 
        self.step = 0

    def go_to_step(self, i):
        """ Move to and read i-th step """
        # switch to seek?
        self.reopen()
        while self.step != i:
            value = self.read_next_step()
        return value

