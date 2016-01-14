import multiprocessing as mp
import copy
from operator import itemgetter

# For fancy progressbar
try:
    from blessings import Terminal
    from progressive.bar import Bar
    from progressive.tree import ProgressTree, Value, BarDescriptor
except:
    pass

class ParallelProcessor():
    
    def __init__(self, list_of_jobs, trajectory, start=None, stop=None, step=None, threads=None,
                       progressbar=True):
        self._trajectory = trajectory

        start, stop, step = trajectory.check_slice_indices(start, stop, step)

        self.start   = start
        self.stop    = stop
        self.step    = step
        self.nframes = len(xrange(start, stop, step))

        self.list_of_jobs = list_of_jobs

        if threads == None:
            self.threads = mp.cpu_count()
        elif threads > mp.cpu_count(): 
            self.threads = mp.cpu_count()
        else:
            self.threads = threads
        
        self.progressbar = progressbar


    def slices(self):
        # This function returns a list containing the start and
        # last configuration to be analyzed from each thread

        threads    = self.threads # Get number of threads from initialization
        step       = self.step
        configs    = (self.stop-self.start)/step

        # Check whether the number of threads is higher than
        # the number of trajectories to be analyzed
        while (configs/threads == 0):
            threads -= 1
            
        self.threads = threads # Update the number of threads

        print "Total configurations: "+str(configs)
        print "Analysis running on ", threads, " threads.\n"

        # Number of cfgs for each thread, and remainder to be added
        thread_cfg   = configs/threads
        reminder   = configs%threads

        slices     = []
        beg = self.start

        # Compute the start and last configurations
        for thread in range(0,threads):
            if thread < reminder:
                end = beg+step*thread_cfg
            else:
                end = beg+step*(thread_cfg-1)
            
            slices.append([beg, end+1])
            beg = end+step

        # Print on screen the configurations assigned to each thread
        for thread in range(threads):
            confs = 1+(slices[thread][1]-1-slices[thread][0])/step
            digits = len(str(self.stop))
            line = "Thread "+str(thread+1).rjust(len(str(threads)))+": " \
                          +str(slices[thread][0]).rjust(digits)+"/"  \
                          +str(slices[thread][1]-1).rjust(digits)    \
                          +" | Configurations: "+str(confs).rjust(1+len(str(thread_cfg)))
            print line

        print 
        return slices

    def job_copies(self):
        # This function creates n=threads copies of the analysis
        # objects submitted from the user
        threads    = self.threads
        job_list   = []

        for thread in range(threads):
            job_list.append([ copy.deepcopy(job) for job in self.list_of_jobs])
            
        return job_list

    def compute(self,job_list,start,stop,step, out_queue, order,progress,error):
        # Run the single_frame method for each analysis object for all
        # the trajectories in the batch
        count = 0 # just to track IO errors

        try:
            for i, ts in enumerate(self._trajectory[start:stop:step]):
                count = i
                for job in job_list:
                    job._single_frame()
                progress.put(order) # Updates the progress bar
            
            out_queue.put((job_list,order)) # Returns the results
        except:
            error.put([order, start+count*step])

    def prepare(self,pool,slices):
        # This function runs the setup_frame method from AnalysisBase
        # and the prepare method from the analysis object. Each job
        # object is initialized using start and stop from the slice function
        for thread, elem in zip(pool,slices):
            for job in thread:
                job._setup_frames(self._trajectory,start=elem[0],stop=elem[1],step=self.step)
                job._prepare()

    def conclude(self,list_of_jobs):
        # Run conclude for each job object
        for job in list_of_jobs:
            job._conclude()
    
    def parallel_run(self):
        # Returns indices for splitting the trajectory
        slices = self.slices()

        # Create copies of the original object to be
        # dispatched to multiprocess
        pool    = self.job_copies()
        threads = self.threads

        self.prepare(pool,slices)

        # Queues for the communication between parent and child processes
        out_queue = mp.Manager().Queue()
        progress  = mp.Manager().Queue()
        error     = mp.Manager().Queue()

        # Prepare multiprocess objects
        processes = [ mp.Process(target=self.compute, args=(batch, elem[0],elem[1],self.step,out_queue,order,progress,error)) 
                      for order, (batch, elem) in enumerate(zip(pool,slices)) ]


        # Run processes
        for p in processes: p.start()

        thread_configs = [ 1+(elem[1]-elem[0]-1)/self.step for elem in slices ]

        if self.progressbar:
            try:
                # Initialize progress bar
                readings  = [ Value(0) for i in range(threads) ]
                total_cfg = Value(0)
                bd_defaults = [ dict(type=Bar, kwargs=dict(max_value=cfg)) for cfg in thread_configs ]
                bd_defaults.append( dict(type=Bar, kwargs=dict(max_value=sum(thread_configs))) ) # add total cfgs for total counter
                
                data_tuple = [ ("Core "+str(thread+1).rjust(len(str(threads))),
                                BarDescriptor(value=readings[thread], **bd_defaults[thread])) for thread in range(threads) ]
                
                data = dict( (key,value) for (key,value) in data_tuple )
                
                data.update({'Total': BarDescriptor(value=total_cfg, **bd_defaults[-1])})
                
                # Create blessings.Terminal instance
                t = Terminal()
                # Initialize a ProgressTree instance
                n = ProgressTree(term=t)
                # Make room for the progress bar
                n.make_room(data)
                
                # Updates progress bar
                while any([ p.is_alive() for p in processes ]):
                    while not progress.empty():
                        core = progress.get()
                        n.cursor.restore()
                        readings[core].value += 1
                        total_cfg.value += 1
                        n.draw(data, BarDescriptor(bd_defaults))
                        
                print # lets leave some free space
            except:
                data = {}
                for thread in range(threads):
                    data[thread]=0
                while any([ p.is_alive() for p in processes ]):
                    while not progress.empty():
                        core = progress.get()
                        data[core] += 1
                        
                        for thread in range(threads):
                            print "Core "+str(thread)+": "+str(data[thread])+"/"+str(thread_configs[thread])
                            
                        print "\r",

                        print "\033["+str(threads+1)+"A"
                    
                print threads*"\n"
                

        # Exit the completed processes
        for p in processes: p.join()

        unread_configurations = []

        # Collect errors from queues
        while not error.empty():
            unread_configurations.append(error.get())

        if len(unread_configurations) != 0:
            unread_counter = 0
            for error in unread_configurations:
                unread_counter += (slices[error[0]][1]-error[1])/self.step
            print "Sorry, there were", unread_counter, "unread configurations."
            for error in sorted(unread_configurations, key=itemgetter(0)):
                print "Core", error[0]+1, "did not read from configuration", error[1]+1,\
                      "to", slices[error[0]][1]
                
            print "\n"


        results = []

        # Collects results from the queue
        while not out_queue.empty():
            results.append(out_queue.get())

        jobs_num = len(self.list_of_jobs)
        
        result_list = []

        # Sum the job objects from each thread
        for job in range(jobs_num):
            for order, thread in enumerate(sorted(results, key=itemgetter(1))):
                if order == 0:
                    result_list.append(thread[0][job])
                else:
                    result_list[job] += thread[0][job]

        # Run the conclude function for each job
        self.conclude(result_list)
