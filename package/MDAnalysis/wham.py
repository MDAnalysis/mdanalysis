# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

import numpy as np
import warnings
import os

import datreant.core as dtr


## TODO - NAMING. Currently named in line with docs for Grossfield wham, but
## some of these aren't very clear/nice so likely to change...
# TODO - should probably split this up
def wham(bundle, spring='spring', loc_win_min='loc_win_min', 
         temperature='temperature', correl_time=None, 
         timeseries_data='timeseries_data', timeseries_type='coord', 
         energy='energy',
         calc_temperature=None, hist_min=None, hist_max=None,
         periodicity='', num_bins=200, tol=1e-6, numpad=0, 
         run_bootstrap=True, num_MC_trials=200,
         start_time=None, end_time=None, 
         energy_units='kcal', keep_files=False):
    """ Wrapper for the Grossfield implementation of WHAM.

    [link documentation]

    Each simulation must have the appropriate metadata (spring constant,
    restrained value, temperature) and auxiliary data (timeseries data as
    either reaction coordinate value, difference from restrained value, or
    value of the restraining force). If using MC bootstrap error analysis,
    can also specify a correlation time for each window. If simulations were
    performed at different temperatures, a potential energy auxiliary must also
    be provided. If any names differ from the default values they must be 
    specified.

    Various wham paramaters, etc etc [TBA]


    Parameters
    ----------
    bundle: datreant Bundle
        Bundle of the simulations WHAM is to be performed for.
    spring: str, optional
        Name of metadata storing the spring constant that each simulation uses,
        assuming biasing potential has form 1/2 k(x-x0)^2. [[<--format]]
        [[Some simulation packages don't use the 1/2 - so have to pass 2*the 
        restraint const instead?]].
        Must match the units used for energy + the reaction coordinate. [examples?].
    loc_win_min : str, optional
        Name of metadata field storing the reaction coordinate value that is
        the minimum of the biasing potential in each window, ie x0 above.
    temperature : str, optional
        Name of metadata field storing each window's temperature (in Kelvin)
    correl_time : str, optional
        Name of metadata field storing decorrelation time for each window, in 
        time step units. Only used for boostrap error; will be set to 1 for 
        each window if not provided.
    timeseries_data : str, optional
        Name of the auxiliary containing the force/reaction coordinate value 
        throughout simulation.
    timeseries_type : str, optional
        What value is recorded in timeseries_data; for available options see
        ``calc_reaction_cood``.
    energy : str, optional
        Name of the auxiliary containing potential energy. [Currently must be
        at exactly the same same steps as timeseries above (ie same dt and 
        initial time)]. [I assume units must match energy_units?]. Only
        required if simulations performed at different temperatures.
    calc_temperature : float, optional
        Temperature at which to perform wham calculation (in Kelvin). If not 
        specified, assume simulation are performed at the same temperature and 
        set to that.
    hist_min, hist_max : float, optional
        Min/max values of reaction coordinate to use in calculation. If None 
        (default), will set to the lowest/highest value of the reaction 
        coordinate.
    perodicity : str, optional
        Periodicity of system. Default ('') indicates a nonperiodic reaction
        coordiante. [[Currently passed straight on to wham, should change]]
    num_bins : int, optional
        Number of bins to use in histogram (= number of points in final profile).
    tol : float, optional
        Reference value to assess convergence (will stop iteration when the 
        biggest change in F-value is less than this).
    numpad : int, optional
        ['padding' values, for periodic PMFs; for nonperiodic use 0 (default)]
    run_boostrap : bool, optional
        Whether to run Monte Carlo bootstrap error analysis.
    num_MC_trials : int, optional
        Number of 'fake' data sets to create, if running bootstrap error analysis
    start_time, end_time : float, optional
        Start/end time (in ps) to use when calculating profile; data outside
        of this time range will be ignored. 
    units : str, optional
        Free energy units to use [[kcal, kJ, ...]]
    keep_files : bool, optional
        Whether to keep files generated by wham
    """
    #### CHECK INPUT
    if not isinstance(bundle, dtr.Bundle):
        TypeError('{} is not a bundle'.format(bundle))

    # TODO - temp catch timeseries_type is valid; will make this nicer
    calc_reaction_coord(1, timeseries_type, 1 ,1)

    if hist_min is not None and hist_max is not None:
        if float(hist_min) >= float(hist_max):   ## check floatable?
            raise ValueError('hist_min {} is greater than hist_max {}'.format(
                                                            hist_min, hist_max))

    # TODO - PERIODICITY. Passed to wham as 'P', 'Ppi', 'P[val]' (or ''); 
    # make input nicer here?

    if start_time is not None and end_time is not None:
        if start_time >= end_time:
            raise ValueError('start time {} is greater than end time {}'.format(
                                                          start_time, end_time))

    # TODO - UNITS. Check option is valid. Then change somehow? Would adjusting
    # temperature work?

    # required metadata/auxiliaries...
    metadata = [spring, loc_win_min, temperature]
    auxiliaries = [timeseries_data]
    if correl_time:
        metadata.append(correl_time)

    check_bundle_metadata(bundle, metadata)
    # check if all simulations are at the same temperature...
    temps = bundle.categories[temperature]
    if all(t == temps[0] for t in temps):
        # all same temperature. This'll be the temperature we calculate at
        multi_temp = False
        if calc_temperature is None:
            calc_temperature = temps[0]
        elif float(calc_temperature) != temps[0]:
            ## would there be a situation where we want them to be different?
            raise ValueError('Simulations all have temperature {} but calc_temperature '
                             '{} does not match'.format(temps[0], calc_temperature))
    else:
        # different temperatures. Will need a run temperature and a potential 
        # energy for every step
        multi_temp = True
        if calc_temperature is None:
            raise TypeError('Must provide a calc_temperature if simulations are '
                            'performed at different temperatures')
        auxiliaries.append(energy) 

    check_bundle_auxiliaries(bundle, auxiliaries)

    # TODO - check values for other options are valid (num_bins, numpad, 
    # num_MC_trials should be int >0; tol (and all temps) should be float >0


    #### FILES
    # TODO - if keeping files, allow to specifiy file names/directory
    # will the default file names always be valid? Particularly if sticking
    # Sim name in it... (could just use index)
    timeseriesfile_root = 'timeseries_{}.dat'
    metadatafile = 'metadatafile.dat'
    outfile = 'outfile.dat' # called freefile in Grossfield wham docs

    #### WRITE INPUT FILES
    ## TODO - getting the aux values is currently a bit nasty, will have
    ## a fiddle with auxiliary stuff to hopefully make this nicer...
    with open(metadatafile, 'w') as meta_file:
        global_min_val = None
        global_max_val = None
        passed_sims=[] #keep track of which simulations we actually feed through
                       # to wham (so we don't try run with none...)
        for sim in bundle:
            timeseries_file = timeseriesfile_root.format(sim.name)
            k = sim.categories[spring]
            x0 = sim.categories[loc_win_min]

            # figure out the time range. Assuming energy will be recorded at same points.
            data = sim.universe.trajectory._auxs[timeseries_data]
            step = 0
            if start_time is None:
                start_step = 0
            else:
                while data.step_to_time(step) < start_time:
                    step = step+1
                    if step == len(data):
                        break
                start_step = step            
            if end_time is None:
                end_step = len(data)
            else:
                while data.step_to_time(step) < end_time:
                    step = step+1
                    if step == len(data):
                        break
                end_step = step
            if start_step == len(data) and end_time is not None:
                warnings.warn('Simulation {} will be skipped (no data before '
                              'end_time ({} ps)'.format(sim.name, end_time))
            elif end_step == 0 and start_time is not None:
                warnings.warn('Simulation {} will be skipped (no data after '
                              'start_time ({} ps)'.format(sim.name, start_time))
            else:
                if end_step ==  len(data):
                    end_step = end_step-1 # TODO temp because errors with len(aux) - need to fix
                # write the timeseries file...
                max_val = None
                min_val = None
                # assuming for now that we have energy/'timeseries data' at
                # the exact same set of points.
                with open(timeseries_file, 'w') as data_file:
                    if multi_temp:
                        ener = sim.universe.trajectory._auxs[energy]
                        ener_data = [i.data[0] for i in ener[start_step:end_step]]
                    for i, as_data in enumerate(data[start_step:end_step]):
                        x = calc_reaction_coord(as_data.data[0], 
                                                      timeseries_type, k, x0)
                        min_val = update_min(x, min_val)
                        max_val = update_max(x, max_val)
                        line = [as_data.time, x]
                        if multi_temp:
                            line.append(ener_data[i])
                        data_file.write(list_to_string(line)+'\n')

                if hist_max is not None and min_val > float(hist_max):
                    warnings.warn('Simulation {} will be skipped (minimum value {} '
                                  'is above upper cutoff {}. Consider lowering '
                              'hist_max.'.format(sim.name, min_val, hist_max))
                elif hist_min is not None and max_val < float(hist_min):
                    warnings.warn('Simulation {} will be skipped (maximum value {} '
                                  'is below lower cutoff {}. Consider increasing '
                                  'hist_min.'.format(sim.name, max_val, hist_min))
                else:
                    # write the line in the metadata file
                    if correl_time:
                        correl = sim.categories[correl_time]
                    else:
                        # only used if doing bootstrap but need as a placeholder if
                        # specifying temperatures; this seems to be the default
                        # value used in wham so should be fine here too
                        correl = 1
                    metafile_info = [timeseries_file, x0, k, correl]
                    if multi_temp:
                        metafile_info.append(sim.categories[temperature])
                    meta_file.write(list_to_string(metafile_info)+'\n')
                    passed_sims.append(sim.name)

                global_min_val = update_min(min_val, global_min_val)
                global_max_val = update_max(max_val, global_max_val)

        if len(passed_sims) == 0:
            raise ValueError('Aborting (all simulations skipped). Try increasing '
                             'time or reaction coordinate range.')

    if hist_min is None:
        hist_min = global_min_val
    if hist_max is None:
        hist_max = global_max_val


    #### RUN!
    wham_command = '~/wham/wham/wham' # TODO - might need to specify path to runfile?
    wham_args = [periodicity, hist_min, hist_max, num_bins, tol, calc_temperature, 
                 numpad, metadatafile, outfile]
    randSeed = 1 # TODO - how to deal with random seed - should make an argument?
    if run_bootstrap:
        wham_args=wham_args+[num_MC_trials, randSeed]
    os.system(wham_command+' '+list_to_string(wham_args))
    ## TODO - switch to subprocess; [+ catch any errors etc]


    #### PARSE OUTPUT FILE
    outfiledata = np.genfromtxt(outfile)
    if not run_bootstrap:
        profile = outfiledata[:,:2]
    if run_bootstrap:
        profile = outfiledata[:,:3]

    #### CLEANUP FILES
    if not keep_files:
        # TODO - best way to remove files? (all in a temp directory?)
        os.remove(metadatafile)
        os.remove(outfile)
        for sim in passed_sims:
            os.remove(timeseriesfile_root.format(sim)) # or use a wildcard

    ####
    return profile
    # TODO - in the outfile we also get the probability + it's error in [:,3] 
    # and [:,4]; and the 'F-values' for each simulation (that we didn't skip);
    # option to get prob instead of PMF? option to return Fvalues as well (tuple?)




def calc_reaction_coord(value, value_type, k, x0):
    """ Calculate value of reaction coordinate corresponding to *value*.

    Calculate the reaction coordinate at a particular time point from *value*,
    depending on *value_type*. Currently allowed types are:
        - coord: reaction coordinate value
        - force: value of the restraining force; harmonic potential is assumed
                 so x = -F/k + x0
        - delta: difference in reaction coord value and minimum of restraining
                 potential; x = delta_x + x0
    [...]
    """
    if value_type == 'coord':
        x = value
    elif value_type == 'force':
        # F = -k delta_x; does the missing 1/2 factor matter here?
        x = -value/k + x0
    elif value_type == 'delta':
        # delta_x = x-x0
        x = value + x0
    else:
        raise ValueError('{} is not a valid timeseries data type'.format(value_type))
    return x

def check_bundle_metadata(bundle, expected):
    """ check each simulation in bundle has the expected metadata """
    common_metadata = bundle.categories.keys()
    for meta in expected:
        if meta not in common_metadata:
            raise ValueError("Not all simulations contain metadata {}."
                             "(Common metadata: {})".format(meta, common_metadata))
        ## TODO - also check if all the values are of the expected type?


def check_bundle_auxiliaries(bundle, expected):
    """ check each simulation in bundle has expected auxiliary """
    # when auxiliaries added to mdsynthesis, this might become more direct
    for aux in expected:
        for sim in bundle:
            if aux not in sim.universe.trajectory.aux_list:
                raise ValueError("Simulation {} does not contain auxiliary data"
                                 " {}".format(sim.name, aux))
            ## TODO - also check it's got the right length/type?

def list_to_string(lst):
    return ' '.join([str(i) for i in lst])

def update_min(new, curr):
    return new if curr is None else new if new < curr else curr

def update_max(new, curr):
    return new if curr is None else new if new > curr else curr

