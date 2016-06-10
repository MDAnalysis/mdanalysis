import numpy as np
import math

class AuxReader(object):
    """ Base class for auxiliary readers.

    Handles iteration and aligning between trajectory timesteps and auxiliary 
    step.

    Reading/parsing of data should be handled by individual format-specific 
    readers.

    When reading in a trajectory, auxiliary data steps are 'assigned' the 
    closest trajectory timestep.

    Attributes
    ==========
      step: 
        number of the current auxiliary step, starting at 0
      dt:   
        change in time between auxiliary steps (ps)
      initial_time: 
        time of first auxilairy step (ps)
      n_cols:
        number of columns of data for each auxiliary step
      time_col:
        index of column in auxiliary data storing time (None if not present)
      time:  
        time of current auxiliary step, as read from data (if present) or 
        calculated using dt and initial_time
      times:
        list of time of each step in auxiliary data
      constant_dt: 
        if true, will use dt/initial_time to calculate time even when time
        stored in data
      data_cols:
        indicies of columns containing data of interest
      step_data:
        values corresponding to data_cols for current step
      
      name: 
        name for auxilairy data; when added to a trajectory the reader will be
        stored as trajectory._auxs[name] and representative timestep data as 
        trajectory.ts.aux 
      represent_ts_as: 
        method to calculate representative value of auxiliary data for a 
        trajectory timestep. Currently available:
            'closest': value from step closest to the trajectory timestep
            'average': average of values from auxiliary steps assigned to 
                       the trajectory timestep
      cutoff:
        auxiliary steps further from the trajectory timestep than *cutoff* 
        will be ignored when calculating representative values
      ts_data:
        list of 'step_data' from each auxiliary step assigned to the current
        trajectory timestep
      ts_rep:
        represenatative value of auxiliary data for current trajectory timestep

    """
      
    def __init__(self, auxname, represent_ts_as='closest', cutoff=None, 
                 dt=None, initial_time=None, time_col=None, data_cols=None, 
                 constant_dt=True, **kwargs):

        self.name = auxname
        self.represent_ts_as = represent_ts_as
        self.cutoff = cutoff

        # set initially to avoid error on first read
        self.n_cols = None

        self.ts_data = None
        self.ts_rep = None

        self._initial_time = initial_time
        self._dt = dt
        self.constant_dt=True
        self.time_col = time_col
        self.data_cols = data_cols

        self.step = -1
        self._read_next_step()
        self.n_cols = len(self._data)

        if time_col >= self.n_cols:
            raise ValueError("Index {0} for time column out of range (num. "
                             "cols is {1})".format(time_col, self.n_cols))
        if data_cols:
            for col in data_cols:
                if col >= self.n_cols:
                    raise ValueError("Index {0} for data column out of range (num."
                                     " cols is {1})".format(col, self.n_cols))

        # get dt and initial time from auxiliary data if stores the step time
        if self.time_col is not None and self.constant_dt:
            self._initial_time = self.time
            self._read_next_step()
            self._dt = self.time - self._initial_time
            self.go_to_first_step()


    def next(self):
        """ Move to next step of auxiliary data """
        return self._read_next_step()

    def __next__(self):
        """ Move to next step of auxiliary data """
        return self.next()

    def __iter__(self):
        self._restart()
        return self

    def _restart(self):
        """ Reset back to start; calling next should read in first step """
        # Overwrite when reading from file to seek to start of file
        self.step = -1
                
    def go_to_first_step(self):
        """ Return to and read first step """
        ## May need to overwrite, depending e.g. how header dealt with
        ## Could also use go_to_step(0) ...
        self._restart()
        self._read_next_step()
                
    def _read_next_step(self):  
        # Define in each auxiliary reader
        raise NotImplementedError(
            "BUG: Override _read_next_timestep() in auxilairy reader!")

    def read_ts(self, ts):
        """ Read the auxiliary steps 'assigned' to *ts* (the steps that are
        within *ts.dt*/2 of of the trajectory timestep/frame - ie. closer to *ts*
        than either the preceeding or following frame). 

        Calculate a 'representative value' for the timestep from the 
        data in each of these auxiliary steps, and add to *ts*.
        """
        # Make sure our auxiliary step starts at the right point (just before
        # the frame being read): the current step should be assigned to a 
        # previous frame, and the next step to either the frame being read of a 
        # following frame. Move to right position if not.
        if not (self.step_to_frame(self.step, ts) < ts.frame
                and self.step_to_frame(self.step+1, ts) >= ts.frame):
            return self.go_to_ts(ts)

        self.reset_ts() # clear previous ts data
        while self.step_to_frame(self.step+1, ts) == ts.frame:
            self._read_next_step()
            self.add_step_to_ts(ts.time)
        self.ts_rep = self.calc_representative()
        ts.aux.__dict__[self.name] = self.ts_rep
        return ts

        # the auxiliary reader should end up positioned at the last step asigned
        # to the trajectory frame (or, if the frame includes no auxiliary steps,
        # (as when auxiliary data is less frequent), the most recent auxiliary 
        # step before the frame)


    def step_to_frame(self, step, ts):
        """ Calculate the frame number that auxiliary step number *step* is 
        assigned to, given dt and offset from *ts*
        (An auxiliary step is assigned to the frame it is closest in time
        to - ) """
        if step >= self.n_steps or step < 0:
            ## make sure step is in the valid range. Raise error?
            return None 
        offset = ts.data.get('time_offset', 0)
        return math.floor((self.times[step]-offset+ts.dt/2.)/ts.dt)

    def go_to_ts(self, ts):
        """ Move to and read auxiliary steps corresponding to *ts* """
        # Need to define in each auxiliary reader
        raise NotImplementedError(
            "BUG: Override go_to_ts() in auxilairy reader!")

    def reset_ts(self):
        self.ts_data = np.array([])
        self.ts_diffs = []

    def add_step_to_ts(self, ts_time):
        """ Add data from the current step to *ts_data* """
        if len(self.ts_data) == 0:
            self.ts_data = np.array([self.step_data])
        else:
            self.ts_data = np.append(self.ts_data, [self.step_data], axis=0)
        self.ts_diffs.append(abs(self.time - ts_time))

    def calc_representative(self):
        """ Calculate a represenatative value from the calues stored in 
        *ts_data* """
        if self.cutoff:
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
        return value
    
    def close(self):
        # Overwrite when reading from file to close open file
        pass    

    @property
    def time(self):
        """ Time in ps of current auxiliary step.
        As read from the appropriate column of the auxiliary data, if present; 
        otherwise calcuated as step * dt + initial_time
        """
        if self.time_col is not None:
            return self._data[self.time_col]
        else:
            return self.step * self.dt + self.initial_time

    @property
    def step_data(self):
        """ Auxiliary values of interest for the current step.
        As taken from the appropriate columns of the full data list *_data*. """
        if self.data_cols:
            return [self._data[i] for i in self.data_cols]
        else:
            return [self._data[i] for i in range(self.n_cols) 
                    if i != self.time_col]

    @property
    def n_steps(self):
        """ Total number of steps in the auxiliary data. """
        try:
            return self._n_steps
        except AttributeError:
            self._n_steps = self.count_n_steps()
            return self._n_steps

    @property
    def times(self):
        """ List of times of each step in the auxiliary data. 
        Calculated using dt and initial_time if constant_dt is true; otherwise
        as read from each step in turn. """
          
        try:
            return self._times
        except AttributeError:
            if self.constant_dt:
                self._times = [i*self.dt+self.initial_time 
                               for i in range(self.n_steps)]
            else:
                self._times = self.read_all_times()
            return self._times

    @property
    def dt(self):
        """ Change in time between steps. Defaults to 1 ps if not provided/
        read from auxiliary data """
        if self._dt:
            return self._dt
        else:
            return 1  ## default to 1ps; WARN?

    @property
    def initial_time(self):
        """ Time corresponding to first auxiliary step. Defaults to 0 ps if 
        not provided/read from auxilairy data """
        if self._initial_time:
            return self._initial_time
        else:
            return 0 ## default to 0; WARN?      


    # TODO - add __enter__ and __exit__ methods when reading from file

    
class AuxFileReader(AuxReader):
    """ Base class for auxiliary readers that read from file 
    Extends AuxReader with methods particular to reading from file"""
    
    def __init__(self, auxname, filename, **kwargs):
        self.auxfilename = filename
        self.auxfile = open(filename)
        
        super(AuxFileReader, self).__init__(auxname, **kwargs)

    def close(self):
        """ close file if open """
        if self.auxfile == None:
            return
        self.auxfile.close()
        self.auxfile = None

    def _restart(self):
        """ reposition to just before first step """
        self.auxfile.seek(0)
        self.step = -1
        
    def _reopen(self):
        self.auxfile.close()
        self.auxfile = open(self.auxfilename)
        self.step = -1


