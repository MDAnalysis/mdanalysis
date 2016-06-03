import numpy as np

class AuxReader(object):
    """ Base class for auxiliary readers
    Does iteration and aligning between trajectory ts/auxiliary step
    Actual reading/parsing of data to be handled by individual readers.
    """
      
    # TODO deal with changing/different units

    def __init__(self, auxnames, represent_ts_as='closest', 
                 cutoff=None, dt=None, time_first_col=False, initial_time=None):
        # TODO default auxnames; allow to name each column + pick which to keep

        self.names = auxnames
        self.represent_ts_as = represent_ts_as
        self.cutoff = cutoff ## UNITS?
        self.time_first_col = time_first_col  ## what if time in another column?

        self.ts_data = None
        self.ts_rep = None

        self.step = 0
        self._read_next_step()

        # set initial time; default to 0
        if time_first_col:
            init_t = self.time
        elif initial_time:
            init_t = intial_time
        else:
            init_t = 0 # TODO - warn using default?
        self.initial_time = init_t

        # set dt (assuming constant!); default to 1 ps.
        if time_first_col:
            self._read_next_step()
            dt = self.time - self.initial_time
            self.go_to_first_step()
        elif dt:
            self.dt = dt
        else:
            init_t = 1 # TODO - warn using default; (check units!)
        self.dt = dt
        

    def next(self):
        """ Move to next step of data """
        return self._read_next_step()

    def __next__(self):
        """ Move to next step of data """
        return self.next()

    def __iter__(self):
        self._restart()
        return self

    def _restart(self):
        """ Reset back to start; calling next should read in first step """
        # Overwrite when reading from file to seek to start of file
        self.step = 0
                
    def go_to_first_step(self):
        """ Return to and read first step """
        ## May need to overwrite, depending e.g. how header dealt with
        ## Could also use go_to_step(1) ...
        self._restart()
        self._read_next_step()
                
    def _read_next_step(self):  
        # Define in each auxiliary reader
        raise NotImplementedError(
            "BUG: Override _read_next_timestep() in auxilairy reader!")

    def read_next_ts(self, ts):
        """ Read and record data from steps closest to *ts*, then 
        calculate representative value for ts """
        # Make sure auxiliary and trajectory are still aligned
        if not self.first_in_ts(ts):
            return self.go_to_ts(ts)
        self.reset_ts()
        while self.step_in_ts(ts):
            self.add_step_to_ts(ts.time)
            try:
                self._read_next_step()
            except StopIteration:
                break
        self.ts_rep = self.calc_representative()
        ts.aux.__dict__[self.names] = self.ts_rep
        return ts
        ## currently means that after reading in a timestep, ts_data and
        ## ts_diffs correspond to that ts but the current step/step_data of 
        ## the auxreader is the first step 'belonging' of the next ts...

    def first_in_ts(self, ts):
        """ Check if current step is first step 'belonging' to *ts*
        Assumes auxiliary *dt* is constant! """
        if (self.time-(ts.time-ts.dt/2)) < self.dt:
            return True
        else:
            return False

    def step_in_ts(self, ts):
        """ Check if current step 'belongs' to *ts* """
        if (self.time-ts.time) <= ts.dt/2. and (self.time-ts.time) > -ts.dt/2.:
            return True
        else:
            return False
           
    def go_to_ts(self, ts):
        """ Move to and read auxilairy steps corresponding to *ts* """
        # Need to define in each auxiliary reader
        raise NotImplementedError(
            "BUG: Override go_to_ts() in auxilairy reader!")

    def reset_ts(self):
        self.ts_data = np.array([])
        self.ts_diffs = []

    def add_step_to_ts(self, ts_time):
        if len(self.ts_data) == 0:
            self.ts_data = [self.step_data]
        else:
            self.ts_data = np.append(self.ts_data, [self.step_data], axis=0)
        self.ts_diffs.append(abs(self.time - ts_time))

    def calc_representative(self):
        if self.cutoff:
            ## starting to get nasty, should probably change...
            cutoff_data = np.array([self.ts_data[i] 
                                    for i,d in enumerate(self.ts_diffs)
                                    if d < self.cutoff])
            cutoff_diffs = [d for d in self.ts_diffs if d < self.cutoff]
        else:
            cutoff_data = self.ts_data
            cutoff_diffs = self.ts_diffs
        if len(cutoff_data) == 0:
            value = [] # TODO - flag as missing
        elif self.represent_ts_as == 'closest':
            value = cutoff_data[np.argmin(cutoff_diffs),:]
        elif self.represent_ts_as == 'average':
            value = np.mean(cutoff_data, axis=0)
        # TODO - interpolation?
        return value
    
    def close():
        # Overwrite when reading from file to close open file
        pass    

    @property
    def time(self):
        """ Time in ps of current auxiliary step
        As read from the auxiliary data, if present; otherwise calcuated
        as (step - 1) * dt + initial_time
        """
        try:
            return self._time
        except AttributeError:
            return (self.step - 1) * self.dt + self.initial_time

    # TODO - add __enter__ and __exit__ methods when reading from file

    
class AuxFileReader(AuxReader):
    """ Base class for auxiliary readers that read from file 
    Extends AuxReader with methods particular to reading from file"""
    
    def __init__(self, auxname, filename, **kwargs):
        super(AuxFileReader, self).__init__(auxname, **kwargs)

        self.auxfilename = filename
        self.auxfile = open(filename)

    def close(self):
        """ close file if open """
        if self.auxfile == None:
            return
        self.auxfile.close()
        self.auxfile = None

    def _restart(self):
        """ reposition to just before first step """
        auxfile.seek(0)
        self.step = 0
        

