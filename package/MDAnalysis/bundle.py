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

"""
Bundle module
=============

This module provides classes for grouping related simulations as datreant
Bundles and storing additional data describing each under common 
names. [[ link datreant ]]

'Additoinal data' may be timeseries (e.g. a distance measured throughout 
the simulation), stored as auxiliary data; or single values (e.g. temperature),
stored in datreant Categories.

Simulations may be added to bundles as part of a list of simulations using
:func:`make_bundle` (optionally adding meta/auxiliary data) or individually (at
a specified position) using :func:`add_to_bundle`. 

#[[ In each case, simulations are added to the Bundle as (persistantly stored) 
#MDSynthesis Sims [[link]]; but may also be passed in as a MDAnalysis Universe, 
#trajectory filename or (topology filename, trajectory filename) tuple; in which
#case a Sim is first created using a provided/default name. ]]
#[[ providing top/univ_kwargs separately; examples; etc... avoid doubling up with 
#func documentation... ]]

The resulting Bundle may the be used as outlined in [MDSynthesis/datreant 
documentation] - e.g. individual simulations may be accessed by index or name::
    bund = make_bundle([Univ1, Univ2], names=['A', 'B'])
    bund[0]                 # returns the Sim for Univ1
    bund['B'].universe      # returns Univ2
Metadata are stored in Sim categories, and so may be accessed/assigned 
individually as e.g. ``bund[0].categories[name]`` or collectively as 
``bund.categories[name]``. 

Auxiliary data is not yet integrated with Sims; they may still be 
added/accessed directly through each universe, e.g. 
``bundle[0].universe.trajectory.add_auxiliary()`` and 
``bundle[0].universe.trajectory.ts.aux``, but currently are not persistently
stored.

.. autofunction:: make_bundle
.. autofunction:: add_to_bundle

"""
import numpy as np

from mdsynthesis import Sim
import datreant.core as dtr

from .core.AtomGroup import Universe

def make_bundle(simulations, topology=None, univ_kwargs=None, names=None, 
                data=None, auxs=None, aux_kwargs=None):
    """ Create a MDSynthesis Bundle from a set of simulations.

    Simulations are passed in as a list; each may of any of the allowed
    simulation types (see [[ ref? or list here anyway? ]]); e.g.::

        make_bundle([Sim1, Sim2, Sim3, ...])
        make_bundle([(top_file1, traj_file1), (top_file2, traj_file2), ...])
        make_bundle([Universe1, Sim2, traj_file3, (top_file4, traj_file4)...])

    If providing just a trajectory filename, a topology file may be specified 
    with the *topology* kwarg, either as a single common topology or a 
    list of files corresonding to each simulation. e.g.:: 

        make_bundle([traj_file1, traj_file2, ... ], topology=common_top_file)
        make_bundle([traj_file1, traj_file2, ...], 
                    topology=[top_file1, top_file2, ...])

    If no topology file is provided, the trajectory file is assumed to double
    as a topology.

    When providing filenames, if additional keyword arguments are required to
    load the corresponding Universe, these can be provided as *univ_kwargs*. 
    Only one set of *univ_kwargs* may be provided, assumed to be common for all;
    otherwise, simulations should be loaded separately before bundling or 
    added separately using :func:`add_to_bundle`.

    If not passing in Sims, these are created, using names passed in with the
    *names* kwarg (as a list). If names are not provided, the index is 
    used instead. 

    Metadata and auxiliary data may be added on initialisation using the *data*  
    and *auxs* kwargs, respectively, passing in the (ordered) list of 
    values/auxiliary filenames etc for each simulation keyworded by metadata/auxiliary
    name; for metadata where all simulations have the same value, a single value
    may be gien in place of a list. e.g.::

        make_bundle([Sim1, Sim2], auxs={'pull_f': ['sim1.xvg', 'sim2.xvg']},
                    data = {'dist':[1,2], 'temp':300})

    If additional kwargs are requried for adding auxiliaries, these may 
    be supplied with *aux_kwargs*::

        make_bundle([Sim1, Sim2], auxs={'pull_f': ['sim1.xvg', 'sim2.xvg']},
                     aux_kwargs={'pull_f': {'dt': 400}})

    Parameters
    ----------
    simulation : lst
        list of simulation to add (ordered if required) [[ list possible types? ]]
    topology : (str, lst), optional
        filename for single topology corresponding to all simulations for which
        only a trajectory filename is provided, or ordered list of topology
        filenames corresponding to each simulation
    univ_kwargs : dict, optional
        dictionary of additional arguments to pass to Universe for each 
        simulation passed as filenames
    names : lst of str
        ordered list of names for the simulation
    data : dict
        dictionary of metadata names and corresponding values to add
    auxs : dixt
        dictionary of auxiliary names and corresponding data to add
    aux_kwargs : dict of dict
        dictionary of additional arguments to use when loading each set of
        auxiliary data
    """
    if names is None:
        # if names not provided, use the index for name
        names = map(str, range(len(simulations)))
    elif len(names) != len(simulations):
        raise ValueError("Number of simulations ({}) and number of names ({}) do "
                         "not match".format(len(simulations), len(names)))

    if not isinstance(topology, (list, np.ndarray)):
        # assume all use the same topology
        tops = [topology]*len(simulations)
    elif len(topology) != len(simulations):
        raise ValueError("Number of simulations ({}) and number of topologies ({}) "
                         "do not match".format(len(simulations), len(topology)))

    bundle = dtr.Bundle()
    for i, sim in enumerate(simulations):
        bundle = add_to_bundle(bundle=bundle, simulation=sim, topology=tops[i],
                               name=names[i])
       

    ## Add any auxiliaries
    auxargs = aux_kwargs if aux_kwargs is not None else {}
    if auxs is not None:
        for auxname, auxvals in auxs.items():
            args = aux_kwargs.get(auxname, {})
            if len(auxvals) != len(simulations):
                raise ValueError("Number of simulations ({}) and number of values"
                                 " for auxiliary {} ({}) do not match".format(
                                       len(simulations), auxname, len(auxvals)))
            for i, sim in enumerate(bundle):
                sim.universe.trajectory.add_auxiliary(auxname=auxname, 
                                                      auxdata=auxvals[i], **args) 

    ## Add any metadata
    data = data if data is not None else {}
    for key, values in data.items():
        if (isinstance(values, (list, np.ndarray)) and 
           len(values) != len(simulations)):
            raise ValueError("Number of simulations ({}) and number of values for"
                             " metadata {} ({}) do not match".format(
                                            len(simulations), key, len(values)))
        bundle.categories[key] = values

    return bundle


def add_to_bundle(bundle=None, simulation=None, name=None, index=None, 
                  topology=None, univ_kwargs=None):
    """ Insert *simulation* into *bundle* at the position *index* and return new
    bundle.

    If *bundle* is None, a new bundle is created containing just *simulation*
    and returned. 

    *simulation* may be [[ ... see above, or list here as well? ]]

    Parameters
    ----------
    bundle : datreant Bundle object, optional
        bundle to insert simulation into
    simulation :
        simulation to add; [[ list possible types? ]]
    name : str
        name for the simulation
    index : int, optional
        position in which to insert *simulation*. If *None* (default), the 
        simulation will be added in the last position.
    topology : str, optional
        filename for corresponding topology if simulation is provided as a 
        trajectory filename
    univ_kwargs : dict, optional
        dictionary of additional arguments to pass to Universe.


    Returns
    -------
    bundle : datreant Bundle object
        updated bundle
    """
    new_sim = _make_sim(simulation, name=name, topology=topology, 
                        univ_kwargs=univ_kwargs)
    if bundle is None:
        # make a new bundle and return
        # TODO - should we warn this is what we're doing?
        return dtr.Bundle(simulation)
    elif isinstance(bundle, dtr.Bundle):
        if index is None:
            bundle.add(new_sim)
        else:
            sims = [sim for sim in bundle]
            sims.insert(index, new_sim)
            bundle = dtr.Bundle(sims)
        return bundle
    else:
         raise TypeError("Invalid bundle {}".format(bundle))


def _make_sim(simulation, name=None, topology=None, univ_kwargs=None):
    """ Return a MDSynthesis Sim for *simulation* with name *name*. """
    univ_kwargs = univ_kwargs or {}
    if isinstance(simulation, Sim):
        # TODO - warn that name will be ignored?
        return simulation
    # TODO - allow to just provide name of existing saved Sim

    # check we have a valid *name* and make a new Sim
    if not isinstance(name, str):
        raise TypeError('Invalid Sim name {}'.format(name))
    sim = Sim(name, new=True)
   
    if isinstance(simulation, str):
        # assume passing in trajectory/topology file
        if isinstance(topology, str):
            sim.universe = Universe(topology, simulation, **univ_kwargs)
        elif topology is None:
            # try load with the trajectory as topology
            # TODO - warn that we're assuming traj acts as top?
            sim.universe = Universe(simulation, **univ_kwargs)
        else:
            raise TypeError('Invalid topology {}'.format(topology))
    elif (isinstance(simulation, tuple) and 
         all(isinstance(x, str) for x in simulation)):
         # assume we're passing in a tuple with top, traj filenames
         sim.universe = Universe(simulation[0], simulation[1], **univ_kwargs)
         # TODO - figure out if the order is wrong?
    elif isinstance(simulation, Universe):
         sim.universe = simulation
    else:
        raise TypeError('Invalid simulation {}; must be MDSynthesis Sim, '
                        'Universe, trajectory filename or tuple of (topology, '
                        'trajectory file names'.format(simulation))

    return sim
    
