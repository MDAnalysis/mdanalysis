""" Add a docstring later

"""
import MDAnalysis as mda
import multiprocessing as mp
import copy
import time
from operator import itemgetter
import utils.prog as prog

class ParallelProcessor(object):
    """ Add class docstring later
    """
    def __init__(self, jobs_list, universe, start=None, stop=None,
                 step=None, threads=None):
        self._universe = universe
        self.topology = universe.filename
        self.trajname = universe._trajectory.filename

        start, stop, step = universe.trajectory.check_slice_indices(start,
                                                                     stop, step)

        self.start = start
        self.stop = stop
        self.step = step
        self.nframes = len(xrange(start, stop, step))

        self.jobs_list = jobs_list

        if threads is None:
            self.threads = mp.cpu_count()
        elif threads > mp.cpu_count():
            self.threads = mp.cpu_count()
        else:
            self.threads = threads

        self.slices = self.compute_slices()

    def compute_slices(self):
        """
        This function returns a list containing the start and
         last configuration to be analyzed from each thread
        """
        threads = self.threads # Get number of threads from initialization
        step = self.step
        configs = 1+(self.stop-self.start-1)/step

        self.nframes = configs

        # Check whether the number of threads is higher than
        # the number of trajectories to be analyzed
        while configs/threads == 0:
            threads -= 1

        self.threads = threads # Update the number of threads

        print "Total configurations: "+str(configs)
        print "Analysis running on ", threads, " threads.\n"

        # Number of cfgs for each thread, and remainder to be added
        thread_cfg = configs/threads
        reminder = configs%threads

        slices = []
        beg = self.start

        # Compute the start and last configurations
        for thread in range(0, threads):
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
                          +" | Configurations: "\
                          +str(confs).rjust(1+len(str(thread_cfg)))
            print line

        return slices


    def compute(self, out_queue, order, progress):
        """
        Run the single_frame method for each analysis object for all
        the trajectories in the batch
        """
        start = self.slices[order][0]
        stop = self.slices[order][1]
        step = self.step

        jobs_list = []
        universe = mda.Universe(self.topology, self.trajname)
        traj = universe.trajectory

        for job in self.jobs_list:
            jobs_list.append(copy.deepcopy(job))

        for job in jobs_list:
            # Initialize job objects. start, stop and step are
            # given so that self.nframes is computed correctly
            job._setup_frames(universe=universe, start=start,
                              stop=stop, step=self.step)
            job._prepare()

        progress.put(order)
        for timestep in traj[start:stop:step]:
            for job in jobs_list:
                job._single_frame(timestep)
            progress.put(order) # Updates the progress bar

        out_queue.put((jobs_list, order)) # Returns the results

    def conclude(self, jobs_list):
        """
        Run conclude for each job object
        """
        for job in jobs_list:
            job._conclude()

    def parallel_run(self):
        """
        Create copies of the original object to be
        dispatched to multiprocess
        """
        threads = self.threads

        # Queues for the communication between parent and child processes
        out_queue = mp.Manager().Queue()
        progress = mp.Manager().Queue()

        # Prepare multiprocess objects
        processes = [mp.Process(target=self.compute,
                                 args=(out_queue, order, progress))
                      for order in range(threads)]


        # Run processes
        for process in processes:
            process.start()

        thread_configs = [1+(elem[1]-elem[0]-1)/self.step
                          for elem in self.slices]

        pb = prog.ProgressbarMulticore(thread_configs,bar_length=50)

        while any([process.is_alive() for process in processes]):
            time.sleep(1)
            pb.timer(1)
            while not progress.empty():
                core = progress.get()
                pb.update(core)
        print

        # Exit the completed processes
        for process in processes:
            process.join()

        results = []

        # Collects results from the queue
        while not out_queue.empty():
            results.append(out_queue.get())

        jobs_num = len(self.jobs_list)

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
