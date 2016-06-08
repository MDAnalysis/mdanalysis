import numpy as np
from . import base

class XVGReader(base.AuxFileReader):
    """ Read data from .xvg file
    
    Assumes data is time-ordered and first column is time
    """
    # TODO - swtich to reading file all at once
 
    def __init__(self, auxname, filename, **kwargs):
        super(XVGReader, self).__init__(auxname, filename, time_col=0, **kwargs)

        
    def _read_next_step(self):
        """ Read next recorded step in xvg file """
        line = self.auxfile.readline()
        if line:
            # xvg has both comments '#' and grace instructions '@'
            while line.lstrip()[0] in ['#', '@']:
                line = self.auxfile.readline()
            # remove end of line comments
            line_no_comment = line.split('#')[0]
            self.step = self.step + 1
            self._data = [float(i) for i in line_no_comment.split()]
            if self.n_cols and len(self._data) != self.n_cols:
                raise ValueError('Step {0} has {1} columns instead of '
                                 '{2}'.format(self.step, len(self._data),
                                              self.n_cols))
            return self._data
        else:
            self.go_to_first_step()
            raise StopIteration
 
    def go_to_ts(self, ts):
        """ Move to and read auxilairy steps corresponding to *ts* """
        # only restart if we're currently beyond *ts*
        if self.step_to_frame(self.step, ts) >= ts.frame:
            self._restart()
        while self.step_to_frame(self.step+1, ts) != ts.frame:
            if self.step == self.n_steps-1:
                return ts
            self._read_next_step()    
 
        return self.read_ts(ts)

    def go_to_step(self, i):
        """ Move to and read i-th auxiliary step """
        ## could seek instead?
        if i >= self.n_steps:
            raise ValueError("{0} is out of range range of auxiliary"
                             "(num. steps {1}!".format(i, self.n_steps))
        if i < 0:
            raise ValueError("Step numbering begins from 0")

        self.go_to_first_step()
        while self.step != i:
            value = self._read_next_step()
        return value

    def count_n_steps(self):
        """ Read through all steps to count total number of steps.
        (Also create times list while reading through) """
        self._restart()
        times = []
        count = 0
        for step in self:
            count = count + 1
            times.append(self.time)
        self._times = times
        return count

    def read_all_times(self):
        self.count_n_steps()

