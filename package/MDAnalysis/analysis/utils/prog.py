""" Add docstring later
"""

import time
import datetime
import numpy as np
import sys
import math

class ProgressbarMulticore(object):
    """ Add class docstring later
    """
    def __init__(self,list_of_totals,name="Analysis", bar_length=40):
        self.threads_compute_units = np.array(list_of_totals)
        self.total_configs = np.sum(self.threads_compute_units)
        self.threads = len(self.threads_compute_units)
        self.working_threads = self.threads

        self.has_started = [False] * self.threads
        self.counter = 0
        self.last_update = np.zeros(self.threads) # times of last update
        self.cumulative_time = 0
        
        self.eta = 0
        self.elaps = 0
        self.speed = 0

        self.remaining_cfgs = np.amax(self.threads_compute_units)
        self.remaining_changed = False

        # cosmectics
        self.name = name
        self.bar_length = bar_length
        self.dots = 0
        self.eta_started = False
        self.cfg_len = len(str(self.total_configs))
        
    def _update_timings(self,core_id):
        istant = time.time()

        if self.has_started[core_id]:
            self.counter += 1
            self.threads_compute_units[core_id] -= 1
            self.cumulative_time += istant-self.last_update[core_id]
            self.last_update[core_id] = istant            
        else:
            self.has_started[core_id] = True
            self.last_update[core_id] = istant

        remainings = np.amax(self.threads_compute_units)
        
        if remainings != self.remaining_cfgs:
            self.remaining_changed = True
            self.remaining_cfgs = remainings

    def _compute_eta(self):
        if self.remaining_changed:
            self.speed = self.cumulative_time/self.counter
            self.eta = self.speed*self.remaining_cfgs
            self.remaining_changed = False

    def update(self,core_id):
        self._update_timings(core_id)
        self._compute_eta()

    def _print_bar(self):
        percentage=self.counter*100./self.total_configs
        bars = int(percentage/100.*self.bar_length)
        empty = self.bar_length-bars

        #timing = " ("+str(round(self.speed,1))+" s/cfg)"
        
        eta = ""

        left_cfgs = " "+str(self.total_configs-self.counter).rjust(self.cfg_len)+"/"+str(self.total_configs).rjust(self.cfg_len)
        
        if (self.eta < 1 and self.eta_started is False):
            eta=str(self.dots*'.')+str((3-self.dots)*' ')
            self.dots += 1
            timing = ""
            if self.dots > 3:
                self.dots = 0
        else:
            self.eta_started = True
            eta=str(datetime.timedelta(seconds=int(self.eta)))

        print self.name+" ["+str(bars*"=")+str(empty*" ")+"] "+str(round(percentage,1)).rjust(4)+"% Elapsed: "+str(datetime.timedelta(seconds=self.elaps))+" ETA: "+eta+left_cfgs+"\r",
        sys.stdout.flush()

    def timer(self,diff):
        if self.eta > diff:
            self.eta -= diff
            self._print_bar()
        else:
            self._print_bar()

        if any(self.has_started):
            self.elaps += diff
