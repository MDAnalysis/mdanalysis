import numpy as np
from . import base

class XVGReader(base.AuxFileReader):
    """ Read data from .xvg file
    
    Assumes data is time-ordered and first column is time

    Currently reading from file; can probably load a once and read through array
    instead
    """
 
    def __init__(self, auxname, filename, **kwargs):
        super(XVGReader, self).__init__(auxname, time_first_col=True, **kwargs)
        ## should generalise in case time not first column...
        
    def _read_next_step(self):
        """ Read next recorded timepoint in xvg file """
        line = self.auxfile.readline()
        if line:
            # xvg has both comments '#' and grace instructions '@'
            while line.lstrip()[0] in ['#', '@']:
                line = self.auxfile.readline()
            # remove end of line comments
            line_no_comment = line.split('#')[0]
            self._time = float(line_no_comment.split()[0])
            self.step_data = [float(i) for i in line_no_comment.split()[1:]]
            # TODO check number of columns is as expected...
            self.step = self.step + 1
            return self.step_data
        else:
            self.go_to_first_step()
            raise StopIteration
 
    def go_to_ts(self, ts):
        """ Move to and read auxilairy steps corresponding to *ts* """
        self.go_to_first_step()
        while not self.step_in_ts(ts):
            self._read_next_step()
        return self.read_next_ts(ts)

    def go_to_step(self, i):
        """ Move to and read i-th step """
        ## probably not the best way to do this - seek?
        self.go_to_first_step()
        while self.step != i:
            value = self._read_next_step()
        return value
        

