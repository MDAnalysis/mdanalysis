from mdsynthesis import Sim
import datreant.core as dtr

from .core.AtomGroup import Universe

def make_bundle(*args, **kwargs):
    """ Create a MDSynthesis Bundle from a set of simulations.

    Simulations may be passed in as a list of MDSynthesis Sims, Universes or
    trajectory files; in the latter case, a single topology file or corresponding
    list of topology files must also be provided as a second argument, so the 
    following are all valid::
        make_bundle([Sim1, Sim2, ...])
        make_bundle([Universe1, Universe2, ...])
        make_bundle([traj_file1, traj_file2, ...], common_topology_file)
        make_bundle([traj_file1, traj_file2, ...], [top_file1, top_file2, ...])

    For the trajectory case, if additional keyword arguments are required, these 
    can be provided with the kwarg *univ_kwargs*. [currently assumes we're using 
    the same kwargs for all].

    If not passing in Sims, these are created, using names passed in 
    with the *names* kwarg (as a list). If names are not provided, the index is 
    used instead. 

    Metadata and auxiliary data may be added on initialisation using the *data*  
    and *auxs* kwargs, respectively, passing in the (ordered) list of 
    values/auxiliary filenames etc for each simulation keyworded by metadata/auxiliary
    name. If additional kwargs are requried for adding auxiliaries, these may 
    be supplied with *aux_kwargs*. e.g.::
        make_bundle([Sim1, Sim2], auxs={'pull_f': ['sim1.xvg', 'sim2.xvg']},
                    data = {'dist':[1,2]}, aux_kwargs={'pull_f': {'dt': 400}})


    The returned Bundle may the be used as outlines in [MDSynthesis/datreant 
    documentation] - e.g. individual simulations may be accessed by index or name::
        bund = make_bundle([Univ1, Univ2], names=['A', 'B'])
        bund[0]                 # returns the Sim for Univ1
        bund['B'].universe      # returns Univ2
    Metadata are stored in Sim categories, and so may be acessed/assigned 
    individually as e.g. ``bund[0].categories[name]`` or collectively as 
    ``bund.categories[name]``. 

    Auxiliary data is not yet integrated with Sims; they may still be 
    added/accessed directly through each universe, e.g. 
    ``bundle[0].universe.trajectory.add_auxiliary()`` and 
    ``bundle[0].universe.trajectory.ts.aux``, but currently are not persistently
    stored.

    """
    # TODO - check length of various lists match
    # TODO - use Group
    # TODO - careful about overwriting existing sims/category values

    names = kwargs.get('names', map(str, range(len(args[0]))))
    if len(args) == 2:
        # assume passing in traj/top
        trajs = args[0]
        tops = args[1]
        uargs = kwargs.get('univ_kwargs', {})
        if isinstance(tops, str):
            # assume we've passed in a single common topology. 
            tops = [tops]*len(trajs)
        bundle = dtr.Bundle([Sim(name, new=True) for name in names])
        for i, sim in enumerate(bundle):
            sim.universe = Universe(tops[i], trajs[i], **uargs)
    elif isinstance(args[0][0], Universe):
        # assume we've passed in list of Universes
        bundle = dtr.Bundle([Sim(name, new=True) for name in names])
        for i, sim in enumerate(bundle):
            sim.universe = args[0][i]
    elif isinstance(args[0][0], Sim):
        # assume we've passed in a list of Sims
        bundle = dtr.Bundle(args[0])

    ## Add any auxiliaries
    auxs = kwargs.get('auxs', {})
    all_aux_args = kwargs.get('aux_kwargs', {})
    for auxname, auxargs in auxs.items():
        aux_args = all_aux_args.get(auxname, {})
        for i, sim in enumerate(bundle):
            sim.universe.trajectory.add_auxiliary(auxname=auxname, 
                                                  auxdata=auxargs[i], 
                                                  **aux_args) 

    ## Add any metadata
    data = kwargs.get('data', {})
    for i, sim in enumerate(bundle):
        sim.categories.add({key: value[i] for key, value in data.items()})

    return bundle
