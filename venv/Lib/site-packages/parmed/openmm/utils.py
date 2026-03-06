"""
This package contains some useful functionality for common tasks in OpenMM
"""
from .. import unit as u
from ..utils.decorators import needs_openmm

def energy_decomposition(structure, context, nrg=u.kilocalories_per_mole):
    """
    This computes the energy of every force group in the given structure and
    computes the energy for each force group for the given Context.  Note, the
    context must have positions already assigned.

    Parameters
    ----------
    structure : Structure
        This should be the Structure object from which the System object in the
        Context was created
    context : mm.Context
        The OpenMM context set up for computing forces and energies
    nrg : energy unit, optional
        The unit to convert all energies into. Default is kcal/mol

    Returns
    -------
    dict {str:float}
        A dictionary mapping the name of the force group (taken from the
        attribute names of the format XXX_FORCE_GROUP in the structure object)
        with the energy of that group in
    """
    all_names = dict()
    force_group_names = dict()
    energy_components = dict()

    for attr in dir(structure):
        if attr.endswith('_FORCE_GROUP'):
            val = getattr(structure, attr)
            all_names[val] = attr.replace('_FORCE_GROUP', '').lower()

    for force in context.getSystem().getForces():
        gp = force.getForceGroup()
        force_group_names[gp] = all_names[gp]

    for grp, name in force_group_names.items():
        state = context.getState(getEnergy=True, groups=1<<grp)
        energy_components[name] = state.getPotentialEnergy().value_in_unit(nrg)

    e = context.getState(getEnergy=True).getPotentialEnergy()
    energy_components['total'] = e.value_in_unit(nrg)

    return energy_components

@needs_openmm
def energy_decomposition_system(structure, system, platform=None,
                                nrg=u.kilocalories_per_mole):
    """
    This function computes the energy contribution for all of the different
    force groups.

    Parameters
    ----------
    structure : Structure
        The Structure with the coordinates for this system
    system : mm.System
        The OpenMM System object to get decomposed energies from
    platform : str
        The platform to use. Options are None (default), 'CPU', 'Reference',
        'CUDA', and 'OpenCL'. None will pick default OpenMM platform for this
        installation and computer
    nrg : energy unit, optional
        The unit to convert all energies into. Default is kcal/mol

    Returns
    -------
    energies : list of tuple
        Each entry is a tuple with the name of the force followed by its
        contribution
    """
    import openmm as mm
    # First get all of the old force groups so we can restore them
    old_groups = [f.getForceGroup() for f in system.getForces()]
    old_recip_group = []
    def _ene(context, grp):
        st = context.getState(getEnergy=True, groups=1<<grp)
        return (type(system.getForce(grp)).__name__,
                st.getPotentialEnergy().value_in_unit(nrg))

    try:
        for i, f in enumerate(system.getForces()):
            if isinstance(f, mm.NonbondedForce):
                old_recip_group.append(f.getReciprocalSpaceForceGroup())
                f.setReciprocalSpaceForceGroup(i)
            f.setForceGroup(i)
        if platform is None:
            con = mm.Context(system, mm.VerletIntegrator(0.001))
        else:
            con = mm.Context(
                system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName(platform)
            )
        con.setPositions(structure.positions)
        if structure.box is not None:
            con.setPeriodicBoxVectors(*structure.box_vectors)

        return list(map(lambda x: _ene(con, x), range(system.getNumForces())))
    finally:
        idx = 0
        for grp, force in zip(old_groups, system.getForces()):
            if isinstance(force, mm.NonbondedForce):
                force.setReciprocalSpaceForceGroup(old_recip_group[idx])
                idx += 1
            force.setForceGroup(grp)
