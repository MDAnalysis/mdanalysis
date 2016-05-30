import numpy as np
from . import base

class XVGReader(base.AuxFileReader):
    """ Read data from .xvg file
    
    Assumes data is time-ordered and first column is time
    """

    def __init__(self, filename, auxname, **kwargs):
        super(XVGReader, self).__init__(filename, auxname, **kwargs)

        self.read_next_step() # TODO read to ts of trajectory aux is being added to
        
    def read_next_step(self):
        """ Read next recorded timepoint in xvg file """
        line = self.auxfile.readline()

        if line:
            # xvg has both comments '#' and grace instructions '@'
            while line.lstrip()[0] in ['#', '@']:
                line = self.auxfile.readline()
            # TODO what about empty lines in header?
            # remove end of line comments
            line_no_comment = line.split('#')[0]
            self.time = float(line_no_comment.split()[0])
            self.step_data = [float(i) for i in line_no_comment.split()[1:]]
            # TODO check number of columns is as expected...
            self.step = self.step + 1
            return self.step_data
        else:
            raise StopIteration
            # can't call next without reseting, but if reopen, messes with 
            # reading timesteps when last ts > last aux



    ## [structure of following a bit iffy atm]
    ## [can move mostly to base]
    def read_next_ts(self, ts):
        """ Read and record data from steps closest to *ts*
        Calculate representative value for ts """
        # TODO make sure starts at 'begining' of ts! --> go_to_ts?
        # TODO fix when reaches EOF during a ts (add except for StopIteration)
        #      and in this case should not reopen when StopIteration...

        self.reset_ts(ts)
        while self.step_in_ts(ts):
            self.add_step_to_ts(ts)
        self.ts_rep = self.calc_rep()
        #ts.aux.__dict__['self.name'] = self.ts_rep   # add aux to ts!
        return self.ts_rep
        ## currently means that after reading in a timestep, ts_data and
        ## ts_diffs correspond to that ts but the current step/step_data of 
        ## the auxreader is the first step of the next ts... change?

    def reset_ts(self, ts):
        self.ts_data = np.array([])
        self.ts_diffs = []

    def step_in_ts(self, ts):
        # or determine from 'dt's; would assume regular spacing...
        if (self.time-ts.time) <= ts.dt/2 and (self.time-ts.time) > -ts.dt/2.:
            return True
        else:
            return False

    def add_step_to_ts(self, ts):
        if len(self.ts_data) == 0:
            self.ts_data = [self.step_data]
        else:
            self.ts_data = np.append(self.ts_data, [self.step_data], axis=0)
        self.ts_diffs.append(abs(self.time - ts.time))
        self.read_next_step()
       
    def calc_rep(self):
        if len(self.ts_data) == 0:
            value = [] #flag as missing
        elif self.repres_ts_as == 'closest':
            value = self.ts_data[np.argmin(self.ts_diffs),:]
            ## TODO - add cutoff
        elif self.repres_ts_as == 'average':
            value = np.mean(self.ts_data, axis=0)
        return value

