"""
This module contains functionality for running simulations with OpenMM using
parameters defined in an input file for sander and pmemd.
"""
import logging
import os
from math import sqrt

from ...amber import AmberMdcrd, AmberMask, NetCDFTraj, Rst7
from ...amber.mdin import Mdin
from ...openmm import (StateDataReporter, NetCDFReporter, MdcrdReporter, RestartReporter,
                       ProgressReporter, EnergyMinimizerReporter)
from ...utils.timer import Timer
from ... import unit as u
from ..exceptions import SimulationError, SimulationWarning, UnhandledArgumentWarning
import sys
import warnings
from ...vec3 import Vec3
from ...geometry import reduce_box_vectors
try:
    from openmm.app import (
        forcefield as ff, OBC1, OBC2, GBn, HCT, GBn2, Simulation, DCDReporter, amberprmtopfile
    )
    import openmm as mm
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False

LOGGER = logging.getLogger(__name__)

_SCRIPT_HEADER = """\
#!/usr/bin/env python3

import os, sys

# Import the OpenMM modules that are necessary
from openmm.app import *
from openmm import *

# Import the Amber/OpenMM modules
from parmed.amber import AmberParm, Rst7, AmberMdcrd, AmberMask, NetCDFTraj
from parmed.openmm import (StateDataReporter, NetCDFReporter, MdcrdReporter,
        RestartReporter, ProgressReporter, EnergyMinimizerReporter)
from parmed import unit as u
from parmed.geometry import reduce_box_vectors

# Load the Amber topology file
parm = AmberParm('%s', '%s')
"""

def positional_restraints(mask, weights, refc, scriptfile=None):
    """
    Creates a positional restraint force from CustomExternalForce using the
    coordinates given in the refc restart file for the atoms specified in 'mask'
    """
    parm = mask.parm # store the parm object
    try:
        refc = refc.coordinates.reshape((-1, len(parm.atoms), 3))
    except ValueError:
        raise SimulationError('Invalid shape of coordinate array')

    # Make the (very simple) force term
    if scriptfile is not None:
        wt = weights.value_in_unit(u.kilocalories_per_mole/u.angstrom**2)
        scriptfile.write("# Create the positional restraint force\n")
        scriptfile.write("refc = refc.reshape((-1, len(parm.atoms), 3))\n")
        scriptfile.write("frc = CustomExternalForce('k*((x-x0)*(x-x0)+"
                         "(y-y0)*(y-y0)+(z-z0)*(z-z0))')\n")
        scriptfile.write("frc.addGlobalParameter('k', %s*u.kilocalorie_per_mole"
                         "/u.angstroms**2)\n" % wt)
        scriptfile.write("frc.addPerParticleParameter('x0')\n")
        scriptfile.write("frc.addPerParticleParameter('y0')\n")
        scriptfile.write("frc.addPerParticleParameter('z0')\n")
        scriptfile.write("for idx in mask.Selected():\n")
        scriptfile.write("    frc.addParticle(idx, refc[0, idx, :] / 10)\n")

    frc = mm.CustomExternalForce('k*((x-x0)*(x-x0)+(y-y0)*(y-y0)+'
                                    '(z-z0)*(z-z0))')
    frc.addGlobalParameter('k', weights)
    frc.addPerParticleParameter('x0')
    frc.addPerParticleParameter('y0')
    frc.addPerParticleParameter('z0')
    for idx in mask.Selected():
        frc.addParticle(idx, refc[0, idx, :] / 10) # angstroms to nanometers

    return frc

def simulate(parm, args):
    """ Runs a simulation using OpenMM """
    global HAS_OPENMM
    if not HAS_OPENMM:
        raise SimulationError('Could not import OpenMM!')

    # Set up timers
    timer = Timer()
    timer.add_timer('system', 'OpenMM System creation')
    # Handle Amber command-line arguments
    overwrite = args.has_key('-O')
    inputfile = args.get_key_string('-i', 'mdin')
    inpcrd = args.get_key_string('-c', None)
    outputfile = args.get_key_string('-o', 'mdout')
    trajectory = args.get_key_string('-x', 'mdcrd')
    mdvel = args.get_key_string('-v', 'mdvel')
    mdfrc = args.get_key_string('-frc', 'mdfrc')
    restart = args.get_key_string('-r', 'restrt')
    refc = args.get_key_string('-ref', 'refc')
    mdinfo = args.get_key_string('-inf', 'mdinfo')
    inptrajname = args.get_key_string('-y', 'inptraj')
    args.get_key_string('-cpin', None)
    args.get_key_string('-cprestrt', None)
    args.get_key_string('-cpout', None)
    args.get_key_string('-remlog', None)
    rem = args.get_key_int('-rem', 0)
    # Some OpenMM-specific options
    plat = args.get_key_string('-platform', None)
    prec = args.get_key_string('-precision', 'mixed')
    dcd = args.has_key('dcd')
    scriptfile = args.get_key_string('script', None)
    runsim = not args.has_key('norun')
    if args.has_key('progress'):
        LOGGER.warning('progress keyword is deprecated - set log level to INFO for progress')
    # Get any unmarked arguments
    unmarked_cmds = args.unmarked()
    if len(unmarked_cmds) > 0:
        warnings.warn("Un-handled arguments: " + ' '.join(unmarked_cmds), UnhandledArgumentWarning)

    # Open up the script file and write the header if requested
    if scriptfile is not None:
        scriptfile = open(scriptfile, 'w')
        if inpcrd is not None:
            scriptfile.write(_SCRIPT_HEADER % (parm, inpcrd))
        else:
            scriptfile.write(_SCRIPT_HEADER % (parm, parm.crdname))

    # Parse the input file
    mdin = Mdin()
    mdin.read(inputfile)

    # We run MD if imin == 0
    runmd = mdin.cntrl_nml['imin'] == 0

    # Process out all of the illegal options that are not supported here
    if rem:
        raise SimulationError('Replica exchange simulations are not supported.')
    if mdin.cntrl_nml['icnstph']:
        raise SimulationError('Constant pH simulations are not supported.')
    if mdin.cntrl_nml['ifqnt']:
        raise SimulationError('QM/MM simulations are not supported.')
    if mdin.cntrl_nml['ievb']:
        raise SimulationError('Empirical Valence Bond simulations not supported.')
    if mdin.cntrl_nml['ipb']:
        raise SimulationError('PB simulations not supported')
    if mdin.cntrl_nml['irism']:
        raise SimulationError('3D-RISM simulations are not supported.')
    if mdin.cntrl_nml['ifcr']:
        raise SimulationError('Charge relocation is not supported.')
    if mdin.cntrl_nml['icfe']:
        raise SimulationError('Thermodynamic integration is not supported.')
    if mdin.cntrl_nml['idecomp']:
        raise SimulationError('Residue decomposition is not supported.')
    if mdin.cntrl_nml['isgld']:
        raise SimulationError('Self-guided langevin dynamics is not supported.')
    if mdin.cntrl_nml['itgtmd']:
        raise SimulationError('Targeted MD is not supported.')
    if mdin.cntrl_nml['ineb']:
        raise SimulationError('Nudged elastic band is not supported.')
    if mdin.cntrl_nml['iamd']:
        raise SimulationError('Accelerated MD is not supported.')
    if mdin.cntrl_nml['ipimd']:
        raise SimulationError('Path-integral MD is not supported.')
    if mdin.cntrl_nml['ilscivr']:
        raise SimulationError('Linearized semi-classical simulations are not supported.')
    if parm.ptr('ifbox') > 0 and mdin.cntrl_nml['igb'] > 0:
        raise SimulationError('Cannot treat periodic system with GB')
    if mdin.cntrl_nml['nmropt']:
        raise SimulationError('NMR restraints are not supported.')
    if runmd and mdin.cntrl_nml['gamma_ln'] == 0 and mdin.cntrl_nml['ntt'] == 3:
        raise SimulationError('Langevin integrator needs a non-zero friction coefficient')
    if mdin.cntrl_nml['ioutfm'] not in (0, 1):
        raise SimulationError('ioutfm must be 0 or 1')
    # Catch some illegal trajectory printing options
    if mdin.cntrl_nml['ioutfm'] == 0 and dcd:
        warnings.warn('Forcefully setting ioutfm=1 to print DCD-formatted trajectories.',
                      SimulationWarning)
        mdin.change('cntrl', 'ioutfm', 1)
    if mdin.cntrl_nml['ntwv'] < 0 or mdin.cntrl_nml['ntwf'] < 0:
        if mdin.cntrl_nml['ioutfm'] != 1 or dcd:
            raise SimulationError('Only NetCDF trajectories can have velocities'
                                  ' and/or forces printed to the same file.')
    if (mdin.cntrl_nml['ntwv'] != 0 or mdin.cntrl_nml['ntwf'] != 0) and dcd:
        raise SimulationError('cannot print velocities or forces to DCD files.')
    if abs(mdin.cntrl_nml['ntwr']) < 500:
        warnings.warn('ntwr is set quite low (%d). Consider much larger values '
                      'to improve computational performance.' %
                      mdin.cntrl_nml['ntwr'], SimulationWarning)

    # Make sure no files will be overwritten unless overwrite is enabled
    if not overwrite:
        ERROR_MESSAGE = '%s exists. Use -O to overwrite.'
        if os.path.exists(outputfile):
            raise SimulationError(ERROR_MESSAGE % outputfile)
        if os.path.exists(mdinfo):
            raise SimulationError(ERROR_MESSAGE % mdinfo)
        if os.path.exists(restart) and mdin.cntrl_nml['imin'] != 5:
            raise SimulationError(ERROR_MESSAGE % restart)
        if os.path.exists(trajectory) and runmd and mdin.cntrl_nml['ntwx'] != 0:
            raise SimulationError(ERROR_MESSAGE % trajectory)
        if os.path.exists(mdvel) and runmd and mdin.cntrl_nml['ntwv'] > 0:
            raise SimulationError(ERROR_MESSAGE % mdvel)
        if os.path.exists(mdfrc) and runmd and mdin.cntrl_nml['ntwf'] > 0:
            raise SimulationError(ERROR_MESSAGE % mdfrc)

    # Set some dependent defaults
    if mdin.cntrl_nml['ntb'] == -1:
        if mdin.cntrl_nml['igb'] == 0:
            if mdin.cntrl_nml['ntp'] > 0:
                mdin.change('cntrl', 'ntb', 2)
            else:
                mdin.change('cntrl', 'ntb', 1)
        else:
            mdin.change('cntrl', 'ntb', 0)
    if not mdin.cntrl_nml['cut']:
        if mdin.cntrl_nml['ntb'] > 0:
            mdin.cntrl_nml['cut'] = 8.0
        else:
            mdin.cntrl_nml['cut'] = 1000.0

    # Determine our cutoff and electrostatic method
    gbmeth, kappa = None, 0.0
    if mdin.cntrl_nml['ntb'] == 0:
        # Interpret cutoffs greater than 500 Angstroms as infinite
        if mdin.cntrl_nml['cut'] >= 500:
            nbmeth = ff.NoCutoff
        else:
            nbmeth = ff.CutoffNonPeriodic
        if mdin.cntrl_nml['igb'] == 1:
            gbmeth = HCT
        elif mdin.cntrl_nml['igb'] == 2:
            gbmeth = OBC1
        elif mdin.cntrl_nml['igb'] == 5:
            gbmeth = OBC2
        elif mdin.cntrl_nml['igb'] == 7:
            gbmeth = GBn
        elif mdin.cntrl_nml['igb'] == 8:
            gbmeth = GBn2
        # Inverse Debye length assuming dielectric of 78.5 and T=298.15K
        # (multiplied by 0.73 to account for lack of ion-ion exclusions)
        kappa = 0.73 * sqrt(mdin.cntrl_nml['saltcon'] * 0.10806)
    else:
        # Periodic
        if mdin.ewald_nml['ew_type'] == 0:
            nbmeth = ff.PME
        elif mdin.ewald_nml['ew_type'] == 1:
            nbmeth = ff.Ewald
        if mdin.ewald_nml['eedmeth'] == 4:
            if mdin.cntrl_nml['cut'] > 500:
                nbmeth = ff.NoCutoff
            else:
                nbmeth = ff.CutoffPeriodic
        elif mdin.ewald_nml['eedmeth'] != 1:
            warnings.warn('eedmeth must be 1 or 4. Other values are ignored.',
                          SimulationWarning)
        if mdin.ewald_nml['vdwmeth'] not in (0, 1):
            raise SimulationError('vdwmeth must be 0 or 1')

    # Determine our constraints
    if mdin.cntrl_nml['ntc'] == 1:
        constraints = None
        rw = False
    elif mdin.cntrl_nml['ntc'] == 2:
        constraints = ff.HBonds
        rw = True
    elif mdin.cntrl_nml['ntc'] == 3:
        constraints = ff.AllBonds
        rw = True
    else:
        raise SimulationError('ntc must be 1, 2, or 3')

    # Determine if these constraints are flexible (their energies are computed)
    if mdin.cntrl_nml['ntf'] == 1:
        flexconst = True
    elif mdin.cntrl_nml['ntf'] == 2:
        flexconst = False
    elif mdin.cntrl_nml['ntf'] > 0 and mdin.cntrl_nml['ntf'] < 8:
        flexconst = False
        warnings.warn('ntf > 2 just implies no constraints are flexible.')
    else:
        raise SimulationError('ntf must be between 0 and 8')

    # Create the OpenMM System object from our current options
    timer.start_timer('system')
    LOGGER.info('Setting up the OpenMM system...')
    system = parm.createSystem(nonbondedMethod=nbmeth,
                     nonbondedCutoff=mdin.cntrl_nml['cut']*u.angstrom,
                     constraints=constraints, rigidWater=rw,
                     implicitSolvent=gbmeth,
                     implicitSolventKappa=kappa*(1.0/u.angstrom),
                     soluteDielectric=mdin.cntrl_nml['intdiel'],
                     solventDielectric=mdin.cntrl_nml['extdiel'],
                     removeCMMotion=mdin.cntrl_nml['nscm']>0,
                     ewaldErrorTolerance=mdin.ewald_nml['rsum_tol'],
                     flexibleConstraints=flexconst,
                     verbose=False
    )
    # Log it if requested
    if scriptfile is not None:
        scriptfile.write("system = parm.createSystem(nonbondedMethod=%s,\n"
                "                nonbondedCutoff=%s*u.angstrom,\n"
                "                constraints=%s, rigidWater=%s,\n"
                "                implicitSolvent=%s,\n"
                "                implicitSolventKappa=%s*(1.0/u.angstrom),\n"
                "                soluteDielectric=%s,\n"
                "                solventDielectric=%s,\n"
                "                removeCMMotion=%s,\n"
                "                ewaldErrorTolerance=%s,\n"
                "                flexibleConstraints=%s,\n"
                "                verbose=False,\n)\n\n" %
                (nbmeth, mdin.cntrl_nml['cut'], constraints, rw, gbmeth, kappa,
                mdin.cntrl_nml['intdiel'], mdin.cntrl_nml['extdiel'],
                mdin.cntrl_nml['nscm']>0, mdin.ewald_nml['rsum_tol'], flexconst)
        )

    # See if we need to turn off the long-range dispersion correction
    if mdin.cntrl_nml['ntb'] > 0 and mdin.ewald_nml['vdwmeth'] == 0:
        # Disable the long-range vdW correction
        for i, frc in enumerate(system.getForces()):
            if (isinstance(frc, mm.NonbondedForce)):
                frc.setUseDispersionCorrection(False)
                if scriptfile is not None:
                    scriptfile.write('# Disable long-range vdW correction\n')
                    scriptfile.write('system.getForces()[%d].setUseDispersionCorrection(False)\n'%i)
            if isinstance(frc, mm.CustomNonbondedForce):
                frc.setUseLongRangeCorrection(False)
                if scriptfile is not None:
                    scriptfile.write('# Disable long-range vdW correction\n')
                    scriptfile.write('system.getForces()[%d].setUseLongRangeCorrection(False)\n'%i)
        if scriptfile is not None:
            scriptfile.write('\n')

    timer.stop_timer('system')

    # See if we need to add positional restraints
    if mdin.cntrl_nml['ntr'] != 0:
        timer.add_timer('restraints', 'Setting up positional restraints')
        timer.start_timer('restraints')
        LOGGER.debug('Setting up the restraints...')
        if mdin.cntrl_nml['restraint_wt'] <= 0:
            raise SimulationError('Restraint weight must be > 0 for restrained dynamics')
        if not mdin.cntrl_nml['restraintmask'].strip():
            raise SimulationError('Restrained atom selection must be made via restraintmask')
        if scriptfile is not None:
            scriptfile.write("# Set up the restrained atom selection mask\n")
            scriptfile.write("mask = AmberMask(parm, '%s')\n" % mdin.cntrl_nml['restraintmask'])
            scriptfile.write('refc = Rst7.open("%s")\n' % refc)
        mask = AmberMask(parm, mdin.cntrl_nml['restraintmask'])
        system.addForce(
            positional_restraints(mask, mdin.cntrl_nml['restraint_wt'] *
                                  u.kilocalorie_per_mole/u.angstrom/u.angstrom,
                                  Rst7.open(refc), scriptfile=scriptfile)
        )
        if scriptfile is not None:
            scriptfile.write('system.addForce(frc)\n\n')
        timer.stop_timer('restraints')

    # See if we need to add a barostat
    if (mdin.cntrl_nml['ntb'] == 2 or mdin.cntrl_nml['ntp'] > 0) and runmd:
        if mdin.cntrl_nml['ntt'] == 0:
            raise SimulationError('constant pressure requires a thermostat')
        if mdin.cntrl_nml['barostat'] == 1:
            raise SimulationError('Berendsen barostat is not implemented in OpenMM. Use the MC '
                                  'barostat instead (barostat=2)')
        if mdin.cntrl_nml['ntp'] == 1:
            # Isotropic scaling
            barostat = mm.MonteCarloBarostat(mdin.cntrl_nml['pres0']*u.bar,
                                             mdin.cntrl_nml['temp0']*u.kelvin,
                                             mdin.cntrl_nml['mcbarint'],
                                            )
            if scriptfile is not None:
                scriptfile.write('barostat = MonteCarloBarostat(%s*u.bar,\n                    '
                                 '%s*u.kelvin, %s\n)\n' % (mdin.cntrl_nml['pres0'],
                                 mdin.cntrl_nml['temp0'], mdin.cntrl_nml['mcbarint'])
                                )
        elif mdin.cntrl_nml['ntp'] == 2:
            # Anisotropic scaling
            barostat = mm.MonteCarloAnisotropicBarostat(mdin.cntrl_nml['pres0']*u.bar,
                                                        mdin.cntrl_nml['temp0']*u.kelvin, True,
                                                        True, True, mdin.cntrl_nml['mcbarint'],
                                                       )
            if scriptfile is not None:
                scriptfile.write('barostat = MonteCarloAnisotropicBarostat(%s*u.bar,\n  '
                                 '               %s*u.kelvin, True, True, True, %s\n)\n' %
                                 (mdin.cntrl_nml['pres0'], mdin.cntrl_nml['temp0'],
                                  mdin.cntrl_nml['mcbarint'])
                                )
        else:
            raise SimulationError('ntp must be 1 or 2 for constant pressure')
        if scriptfile is not None:
            scriptfile.write('system.addForce(barostat)\n\n')
        system.addForce(barostat)

    # Set the integrator (we need this to set constraint parameters even in
    # minimization if someone requested SHAKEn minimization). Then set the
    # tolerance and get the simulation platform
    if mdin.cntrl_nml['ntt'] == 3:
        integrator = mm.LangevinIntegrator(mdin.cntrl_nml['temp0']*u.kelvin,
                                           mdin.cntrl_nml['gamma_ln']/u.picosecond,
                                           mdin.cntrl_nml['dt']*u.picosecond
                                          )
        if scriptfile is not None:
            scriptfile.write('integrator = LangevinIntegrator(%s*u.kelvin, '
                             '%s/u.picosecond,\n                    %s*u.picosecond\n)\n' %
                             (mdin.cntrl_nml['temp0'], mdin.cntrl_nml['gamma_ln'],
                              mdin.cntrl_nml['dt'])
                            )
    else:
        if mdin.cntrl_nml['ntt'] == 1:
            raise SimulationError('ntt must be 2 or 3 for OpenMM (Andersen '
                                  'thermostat or Langevin dynamics')
        elif mdin.cntrl_nml['ntt'] == 2:
            if scriptfile is not None:
                scriptfile.write('system.addForce(AndersenThermostat(%s*u.kelvin, %s\n)\n' %
                                 (mdin.cntrl_nml['temp0'], mdin.cntrl_nml['vrand']))
            system.addForce(
                mm.AndersenThermostat(mdin.cntrl_nml['temp0']*u.kelvin, mdin.cntrl_nml['vrand'])
            )
        if scriptfile is not None:
            scriptfile.write('integrator = VerletIntegrator(%s*u.picosecond)\n'
                             % mdin.cntrl_nml['dt'])
        integrator = mm.VerletIntegrator(mdin.cntrl_nml['dt']*u.picosecond)

    # Set the constraint tolerance on the integrator
    if constraints is not None:
        integrator.setConstraintTolerance(mdin.cntrl_nml['tol'])
        if scriptfile is not None:
            scriptfile.write('integrator.setConstraintTolerance(%s)\n' % mdin.cntrl_nml['tol'])

    # Define the platform if it was requested, or just take the default
    if plat is not None:
        if plat not in ('CUDA', 'Reference', 'CPU', 'OpenCL'):
            raise SimulationError('Platform must be one of CUDA, Reference, CPU, or OpenCL')
        platform = mm.Platform.getPlatformByName(plat)
        if scriptfile is not None:
            scriptfile.write('platform = Platform.getPlatformByName("%s")\n' % plat)
    else:
        if scriptfile is not None:
            scriptfile.write('platform = None\n')
        platform = None

    if prec not in ('mixed', 'double', 'single'):
        raise SimulationError('Invalid precision model [%s]. Must be single, '
                              'double, or mixed.' % prec)

    # Create the Simulation object
    timer.add_timer('simulation', 'Creating the Simulation object')
    timer.start_timer('simulation')
    LOGGER.info('Setting up the simulation...')
    if scriptfile is not None:
        scriptfile.write('simulation = Simulation(parm.topology, system, integrator, platform)\n\n')
    simulation = Simulation(parm.topology, system, integrator, platform)

    # Now set the default precision model. This needs to be done based on the
    # chosen platform. To handle the case where we chose the default platform,
    # we need to get the platform from the current context
    if platform is None:
        platform = simulation.context.getPlatform()
        if scriptfile is not None:
            scriptfile.write('platform = simulation.context.getPlatform()\n')
    if platform.getName() == 'CUDA':
        if scriptfile is not None:
            scriptfile.write('platform.setPropertyValue(simulation.context, '
                             "'CudaPrecision', '%s')\n\n" % prec)
        platform.setPropertyValue(simulation.context, 'CudaPrecision', prec)
    elif platform.getName() == 'OpenCL':
        if scriptfile is not None:
            scriptfile.write('platform.setPropertyValue(simulation.context, '
                             "'OpenCLPrecision', '%s')\n\n" % prec)
        platform.setPropertyValue(simulation.context, 'OpenCLPrecision', prec)
    elif platform.getName() == 'Reference' and prec != 'double':
        warnings.warn('The Reference platform only uses double precision', SimulationWarning)
    elif platform.getName() == 'CPU' and prec != 'single':
        warnings.warn('The CPU platform only uses single precision', SimulationWarning)

    # Set the particle positions and box vectors (if applicable) from either the
    # given restart file or from the topology file object. Also use particle
    # velocities if requested via irest=1
    LOGGER.info('Setting up initial coordinates and velocities...')
    if inpcrd is not None:
        position_container = Rst7.open(inpcrd)
        if position_container.natom != parm.ptr('natom'):
            raise SimulationError('inpcrd [%s] and prmtop mismatch in number of atoms' % inpcrd)
    else:
        position_container = parm
    if scriptfile is not None:
        scriptfile.write('# Set the positions and box vectors\n'
                         'simulation.context.setPositions(parm.positions)\n')
    simulation.context.setPositions(position_container.positions)
    if parm.ptr('ifbox') > 0:
        # Only set box vectors if box is present
        if scriptfile is not None:
            scriptfile.write('reduced_box_vecs = reduce_box_vectors(*position_container.box_vectors)*u.angstroms\n'
                             'simulation.context.setPeriodicBoxVectors(*reduced_box_vecs)\n\n')
        reduced_box_vecs = reduce_box_vectors(*position_container.box_vectors)*u.angstroms
        simulation.context.setPeriodicBoxVectors(*reduced_box_vecs)

    # Velocities
    if runmd and mdin.cntrl_nml['irest'] == 1 and position_container.hasvels:
        if scriptfile is not None:
            scriptfile.write('# Set velocities\n')
            scriptfile.write('simulation.context.setVelocities(parm.velocities/10)\n')
        simulation.context.setVelocities(position_container.velocities/10) # /10 for Å/ps to nm/ps
    elif runmd:
        if scriptfile is not None:
            scriptfile.write('# Set velocities\n')
            scriptfile.write('simulation.context.setVelocitiesToTemperature(\n'
                             '           %s*u.kelvin)\n' % mdin.cntrl_nml['tempi'])
        simulation.context.setVelocitiesToTemperature(mdin.cntrl_nml['tempi']*u.kelvin)
    # Add the energy reporter
    if scriptfile is not None:
        scriptfile.write('# Add the state data reporters\n')
        density = runmd and mdin.cntrl_nml['ntp'] > 0
        scriptfile.write('rep = StateDataReporter("%s", %s, volume=%s, density=%s)\n' %
                         (outputfile, mdin.cntrl_nml['ntpr'], parm.ptr('ifbox') > 0, density)
                        )
        scriptfile.write('simulation.reporters.append(rep)\n')
        scriptfile.write('rep = ProgressReporter("%s", %s, %s, volume=%s, density=%s)\n' %
                         (mdinfo, mdin.cntrl_nml['ntpr'], mdin.cntrl_nml['nstlim'],
                          parm.ptr('ifbox') > 0, density)
                        )
        scriptfile.write('simulation.reporters.append(rep)\n')
    rep = StateDataReporter(outputfile, mdin.cntrl_nml['ntpr'], volume=parm.ptr('ifbox') > 0,
                            density=density)
    simulation.reporters.append(rep)
    rep = ProgressReporter(mdinfo, mdin.cntrl_nml['ntpr'], mdin.cntrl_nml['nstlim'],
                           volume=mdin.cntrl_nml['ntb'] > 0, density=density)
    simulation.reporters.append(rep)

    timer.stop_timer('simulation')

    # Now see if we wanted minimization or not
    if mdin.cntrl_nml['imin'] > 0:
        # We're doing minimization
        if mdin.cntrl_nml['imin'] == 1:
            timer.add_timer('minimization', 'Structure minimization')
            timer.start_timer('minimization')
            LOGGER.info('Minimizing...')
            if scriptfile is not None:
                scriptfile.write('rep = EnergyMinimizerReporter("%s", volume=%s)\n' %
                                 (outputfile, parm.ptr('ifbox') > 0))
                scriptfile.write('rep.report(simulation)\n')
                scriptfile.write('# Minimize the energy\n')
                scriptfile.write('simulation.minimizeEnergy(\n'
                                 '    tolerance=%s*u.kilocalories_per_mole,\n'
                                 '    maxIterations=%s\n)\n' %
                                 (mdin.cntrl_nml['drms'], mdin.cntrl_nml['maxcyc']))
                scriptfile.write('rep.report(simulation)\n\n')
                scriptfile.write('# Write the restart file\n')
                scriptfile.write('restrt_reporter = RestartReporter("%s", 1,\n'
                                 '    False, %s, write_velocities=False\n)\n' %
                                 (restart, mdin.cntrl_nml['ntxo']==2))
                scriptfile.write('restrt_reporter.report(simulation,\n'
                                 '    simulation.context.getState(getPositions=True, '
                                 'enforcePeriodicBox=%s)\n)\n' % bool(mdin.cntrl_nml['ntb']))
            rep = EnergyMinimizerReporter(outputfile, volume=parm.ptr('ifbox') > 0)
            rep.report(simulation)
            simulation.minimizeEnergy(tolerance=mdin.cntrl_nml['drms']*u.kilocalories_per_mole,
                                      maxIterations=mdin.cntrl_nml['maxcyc'])
            rep.report(simulation)
            # Write a restart file with the new coordinates
            restrt_reporter = RestartReporter(restart, 1, False,
                                              mdin.cntrl_nml['ntxo'] == 2, write_velocities=False)
            restrt_reporter.report(simulation, simulation.context.getState(getPositions=True,
                                   enforcePeriodicBox=bool(mdin.cntrl_nml['ntb'])))
            timer.stop_timer('minimization')
        elif mdin.cntrl_nml['imin'] == 5:
            timer.add_timer('minimization', 'Structure minimization')
            timer.start_timer('minimization')
            LOGGER.info('Minimizing...')
            try:
                inptraj = NetCDFTraj.open_old(inptrajname)
                nframes = inptraj.frame
                if scriptfile is not None:
                    scriptfile.write('inptraj = NetCDFTraj.open_old("%s")\n' % inptrajname)
                    scriptfile.write('nframes = inptraj.frame\n')
            except RuntimeError:
                inptraj = AmberMdcrd(inptrajname, parm.ptr('natom'), mdin.cntrl_nml['ntb'] > 0,
                                     mode='r')
                nframes = len(inptraj.coordinates())
                if scriptfile is not None:
                    scriptfile.write('inptraj = AmberMdcrd("%s", %s, %s, mode="r")\n' %
                                     (inptrajname, "parm.ptr('natom')", mdin.cntrl_nml['ntb']>0))
                    scriptfile.write('nframes = inptraj.frame\n')
            # Create a reporter to handle the minimized coordinates
            if mdin.cntrl_nml['ioutfm'] == 0:
                crd_reporter = MdcrdReporter(trajectory, 1)
                if scriptfile is not None:
                    scriptfile.write('crd_reporter = MdcrdReporter("%s", 1)\n' % trajectory)
            elif dcd:
                crd_reporter = DCDReporter(trajectory, 1)
                if scriptfile is not None:
                    scriptfile.write('crd_reporter = DCDReporter("%s", 1)\n' % trajectory)
            else:
                crd_reporter = NetCDFReporter(trajectory, 1)
                if scriptfile is not None:
                    scriptfile.write('crd_reporter = NetCDFReporter("%s", 1)\n' % trajectory)
            if scriptfile is not None:
                scriptfile.write('f = open("%s", "w", 0)\n'
                                 'rep = EnergyMinimizerReporter(f, volume=%s)\n'
                                 'for frame in range(nframes)\n'
                                 '    crds = inptraj.coordinates(frame)\n'
                                 '    simulation.context.setPositions(\n'
                                 '         tuple([Vec3(crds[3*i], crds[3*i+1], '
                                 'crds[3*i+2])\n'
                                 '           for i in range(parm.ptr("natom"))]) * '
                                 'u.angstroms\n    )\n' % outputfile)
            f = open(outputfile, 'w', 0)
            rep = EnergyMinimizerReporter(f, volume=parm.ptr('ifbox') > 0)
            for frame in range(nframes):
                crds = inptraj.coordinates(frame)
                simulation.context.setPositions(
                        tuple([Vec3(crds[3*i], crds[3*i+1], crds[3*i+2])
                        for i in range(parm.ptr('natom'))]) * u.angstroms
                )
                rep.report(simulation, frame=frame+1)
                if mdin.cntrl_nml['maxcyc'] > 1:
                    if scriptfile is not None:
                        scriptfile.write('    simulation.minimizeEnergy(\n'
                                         '          tolerance=%s*u.kilocalories_per_mole,\n'
                                         '          maxIterations=%s,\n'
                                         '    )\n'
                                         '    crd_reporter.report(simulation,\n'
                                         '        simulation.context.getState('
                                         'getPositions=True),\n    )\n'
                                         '    rep.report(simulation)\n' %
                                         (mdin.cntrl_nml['drms'], mdin.cntrl_nml['maxcyc']))
                    simulation.minimizeEnergy(
                        tolerance=mdin.cntrl_nml['drms']* u.kilocalories_per_mole,
                        maxIterations=mdin.cntrl_nml['maxcyc']
                    )
                    crd_reporter.report(simulation, simulation.context.getState(getPositions=True))
                    rep.report(simulation)
                f.write('=' * 80 + '\n')
            f.close()
            if scriptfile is not None:
                scriptfile.write('    f.write("=" * 80 + "\\n")\n')
                scriptfile.write('f.close()\n')
            timer.stop_timer('minimization')
    # Otherwise we want to do dynamics
    else:
        timer.add_timer('md', 'Molecular Dynamics')
        timer.start_timer('md')
        # Coordinate trajectory reporters
        if mdin.cntrl_nml['ntwx'] > 0:
            if mdin.cntrl_nml['ioutfm'] == 0:
                simulation.reporters.append(
                        MdcrdReporter(trajectory, mdin.cntrl_nml['ntwx'])
                )
                if scriptfile is not None:
                    scriptfile.write('simulation.reporters.append(\n'
                                     '      MdcrdReporter("%s", %s)\n)\n' %
                                     (trajectory, mdin.cntrl_nml['ntwx']))
            elif dcd:
                simulation.reporters.append(DCDReporter(trajectory, mdin.cntrl_nml['ntwx']))
                if scriptfile is not None:
                    scriptfile.write('simulation.reporters.append(\n'
                                     '      DCDReporter("%s", %s)\n)\n' %
                                     (trajectory, mdin.cntrl_nml['ntwx']))
            else:
                simulation.reporters.append(
                   NetCDFReporter(trajectory, mdin.cntrl_nml['ntwx'], crds=True,
                                  vels=mdin.cntrl_nml['ntwv'] < 0, frcs=mdin.cntrl_nml['ntwf'] < 0)
                )
                if scriptfile is not None:
                    scriptfile.write('simulation.reporters.append(\n'
                                     '      NetCDFReporter("%s", %s,\n'
                                     '          crds=True, vels=%s, frcs=%s\n'
                                     '      )\n)\n' % (trajectory, mdin.cntrl_nml['ntwx'],
                                     mdin.cntrl_nml['ntwv'] < 0, mdin.cntrl_nml['ntwf'] < 0))
        # Velocity trajectory reporters
        if mdin.cntrl_nml['ntwv'] > 0:
            if mdin.cntrl_nml['ioutfm'] == 0:
                simulation.reporters.append(
                   MdcrdReporter(mdvel, mdin.cntrl_nml['ntwv'], crds=False, vels=True, frcs=False)
                )
                if scriptfile is not None:
                    scriptfile.write('simulation.reporters.append(\n'
                                     '      MdcrdReporter("%s", %s, crds=False,\n'
                                     '         vels=True, frcs=False)\n'
                                     ')\n' % (trajectory, mdin.cntrl_nml['ntwv']))
            elif dcd:
                # Should never make it here due to earlier checks; catch anyway
                raise SimulationError('Cannot write velocities to DCD files.')
            else:
                # Add forces too only if ntwf < 0 AND no coordinate traj is
                # being printed (otherwise forces will be printed to the
                # coordinate file
                _frcs = mdin.cntrl_nml['ntwf'] < 0 and mdin.cntrl_nml['ntwx'] == 0
                simulation.reporters.append(
                    NetCDFReporter(mdvel, mdin.cntrl_nml['ntwv'], crds=False, vels=True, frcs=_frcs)
                )
                if scriptfile is not None:
                    scriptfile.write('simulation.reporters.append(\n'
                                     '      NetCDFReporter("%s", %s,\n'
                                     '          crds=False, vels=True, frcs=%s\n'
                                     '      )\n)\n' % (mdvel, mdin.cntrl_nml['ntwv'],
                                      mdin.cntrl_nml['ntwf'] < 0 and mdin.cntrl_nml['ntwx'] == 0))
        # Force trajectory reporters
        if mdin.cntrl_nml['ntwf'] > 0:
            if mdin.cntrl_nml['ioutfm'] == 0:
                simulation.reporters.append(
                    MdcrdReporter(mdfrc, mdin.cntrl_nml['ntwf'], crds=False, vels=False, frcs=True)
                )
                if scriptfile is not None:
                    scriptfile.write('simulation.reporters.append(\n'
                                     '      MdcrdReporter("%s", %s, crds=False,\n'
                                     '         vels=False, frcs=True)\n'
                                     ')\n' % (trajectory, mdin.cntrl_nml['ntwf']))
            elif dcd:
                # Should never make it here due to earlier checks; catch anyway
                raise SimulationError('Cannot write forces to DCD files.')
            else:
                simulation.reporters.append(
                    NetCDFReporter(mdfrc, mdin.cntrl_nml['ntwf'], crds=False, vels=False, frcs=True)
                )
                if scriptfile is not None:
                    scriptfile.write('simulation.reporters.append(\n'
                                     '      NetCDFReporter("%s", %s,\n'
                                     '          crds=False, vels=False, frcs=True\n'
                                     '      )\n)\n' % (mdfrc, mdin.cntrl_nml['ntwv']))
        # Restart file reporter
        restrt_reporter = RestartReporter(restart, abs(mdin.cntrl_nml['ntwr']),
                                          mdin.cntrl_nml['ntwr'] < 0, mdin.cntrl_nml['ntxo'] == 2)
        if scriptfile is not None:
            scriptfile.write('restrt_reporter = RestartReporter("%s", %s,\n      %s, %s)\n' %
                             (restart, abs(mdin.cntrl_nml['ntwr']), mdin.cntrl_nml['ntwr'] < 0,
                              mdin.cntrl_nml['ntxo'] == 2))
        if mdin.cntrl_nml['ntwr'] != 0:
            simulation.reporters.append(restrt_reporter)
            if scriptfile is not None:
                scriptfile.write('simulation.reporters.append(restrt_reporter)\n')

        LOGGER.info('Running MD...')
        if scriptfile is not None:
            scriptfile.write('simulation.step(%s)\n'
                             'final_state = simulation.context.getState(getPositions=True,\n'
                             '     getVelocities=True, enforcePeriodicBox=%s)\n'
                             'restrt_reporter.report(simulation, final_state)\n' %
                             (mdin.cntrl_nml['nstlim'], mdin.cntrl_nml['ntb'] > 0))
        if runsim:
            simulation.step(mdin.cntrl_nml['nstlim'])
            # Now write the final restart file, hijacking the restart file
            # reporter
            final_state = simulation.context.getState(getPositions=True, getVelocities=True,
                                                      enforcePeriodicBox=mdin.cntrl_nml['ntb'] > 0)
            restrt_reporter.report(simulation, final_state)
            timer.stop_timer('md')
            nsperday = (mdin.cntrl_nml['nstlim'] * mdin.cntrl_nml['dt'] /
                        1000.0 / (timer.timers['md'] / 3600 / 24))

    timer.done()
    LOGGER.info('Done')
    for t in timer.timer_names:
        timer.print_(t, sys.stdout)

    if mdin.cntrl_nml['imin'] == 0 and runsim:
        LOGGER.info('MD timing: %.3f ns/day' % nsperday)

    if scriptfile is not None:
        scriptfile.close()

def energy(parm, args, output=sys.stdout):
    """
    Computes a single-point energy using OpenMM and prints the result to the
    desired output destination (defaults to stdout)
    """
    import tempfile
    global HAS_OPENMM
    if not HAS_OPENMM:
        raise SimulationError('Could not import OpenMM!')

    cutoff = args.get_key_float('cutoff', None)
    igb = args.get_key_int('igb', 5)
    saltcon = args.get_key_float('saltcon', 0.0)
    do_ewald = args.has_key('Ewald')
    plat = args.get_key_string('platform', None)
    prec = args.get_key_string('precision', 'mixed')
    decomp = args.has_key('decompose')
    applayer = args.has_key('applayer')
    vdw_longrange = not args.has_key('nodisper')
    # Get any unmarked arguments
    unmarked_cmds = args.unmarked()
    if len(unmarked_cmds) > 0:
        warnings.warn("Un-handled arguments: " + ' '.join(unmarked_cmds), UnhandledArgumentWarning)

    gbmeth = None
    if parm.ptr('ifbox') == 0:
        if cutoff is None or cutoff >= 500:
            nbmeth = ff.NoCutoff
            cutoff = 1000.0
        else:
            nbmeth = ff.CutoffNonPeriodic
        if not igb in (0, 1, 2, 5, 6, 7, 8):
            raise SimulationError('Bad igb value. Must be 0, 1, 2, 5, 6, 7, or 8')
        if igb == 1:
            gbmeth = HCT
        elif igb == 2:
            gbmeth = OBC1
        elif igb == 5:
            gbmeth = OBC2
        elif igb == 7:
            gbmeth = GBn
        elif igb == 8:
            gbmeth = GBn2
    else:
        cutoff = cutoff if cutoff is not None else 8.0
        if do_ewald:
            nbmeth = ff.Ewald
        else:
            nbmeth = ff.PME

    parm_ = parm
    if applayer:
        # Write out a temporary topology file, load an amberprmtopfile, then
        # delete that file
        fd, tmp = tempfile.mkstemp(suffix='.parm7')
        os.close(fd)
        try:
            parm.write_parm(tmp)
            parm_ = amberprmtopfile.AmberPrmtopFile(tmp)
            os.unlink(tmp)
        except IOError:
            raise SimulationError('Could not create temporary file for app ' # pragma: no cover
                                  'layer energy calculation.')

    # Time to create the OpenMM system
    system = parm_.createSystem(nonbondedMethod=nbmeth, nonbondedCutoff=cutoff*u.angstrom,
                                constraints=None, rigidWater=True, removeCMMotion=False,
                                implicitSolvent=gbmeth, implicitSolventSaltConc=saltcon*u.molar,
                                soluteDielectric=1.0, solventDielectric=78.5,
                                ewaldErrorTolerance=5e-5)

    # If we used the app layer, we need to assign force groups to enable energy
    # decomposition
    if applayer:
        for force in system.getForces():
            if isinstance(force, mm.HarmonicBondForce):
                force.setForceGroup(parm.BOND_FORCE_GROUP)
            elif isinstance(force, mm.HarmonicAngleForce):
                force.setForceGroup(parm.ANGLE_FORCE_GROUP)
            elif isinstance(force, mm.PeriodicTorsionForce):
                force.setForceGroup(parm.DIHEDRAL_FORCE_GROUP)
            else:
                # Treat the rest as nonbonded forces
                force.setForceGroup(parm.NONBONDED_FORCE_GROUP)
        # For periodic simulations, we need to set the box info
        if parm.ptr('ifbox') > 0:
            system.setDefaultPeriodicBoxVectors(*parm.box_vectors)

    # Now see if we want to turn on or off the dispersion correction
    for force in system.getForces():
        if isinstance(force, mm.NonbondedForce):
            force.setUseDispersionCorrection(vdw_longrange)
        if isinstance(force, mm.CustomNonbondedForce):
            force.setUseLongRangeCorrection(vdw_longrange)

    # Create a dummy integrator
    integrator = mm.VerletIntegrator(2.0)

    # Define the platform if it was requested, or just take the default
    if plat is not None:
        if plat not in ('CUDA', 'Reference', 'CPU', 'OpenCL'):
            raise SimulationError('Platform must be one of CUDA, Reference, CPU, or OpenCL')
        platform = mm.Platform.getPlatformByName(plat)
    else:
        platform = None

    if prec not in ('mixed', 'double', 'single'):
        raise SimulationError('Invalid precision model [%s]. Must be single, '
                              'double, or mixed.' % prec)

    # Create the Simulation object
    if platform is None:
        context = mm.Context(system, integrator)
    else:
        context = mm.Context(system, integrator, platform)

    # Set the particle positions
    context.setPositions(parm.positions)

    # Now set the default precision model. This needs to be done based on the
    # chosen platform. To handle the case where we chose the default platform,
    # we need to get the platform from the current context
    if platform is None:
        platform = context.getPlatform()
    if platform.getName() == 'CUDA':
        platform.setPropertyValue(context, 'CudaPrecision', prec)
    elif platform.getName() == 'OpenCL':
        platform.setPropertyValue(context, 'OpenCLPrecision', prec)
    elif platform.getName() == 'Reference' and prec != 'double':
        warnings.warn('The Reference platform only uses double precision', SimulationWarning)
    elif platform.getName() == 'CPU' and prec != 'single':
        warnings.warn('The CPU platform only uses single precision', SimulationWarning)

    # Now get the energy
    has_pbc = parm.ptr('ifbox') > 0
    nrg = u.kilocalories/u.mole
    if decomp:
        # Now we have to loop through every force group and compute them individually

        # Bonds first
        state = context.getState(getEnergy=True, enforcePeriodicBox=has_pbc,
                                 groups=2**parm.BOND_FORCE_GROUP)
        bond = state.getPotentialEnergy().value_in_unit(nrg)
        # Now angles
        state = context.getState(getEnergy=True, enforcePeriodicBox=has_pbc,
                                 groups=2**parm.ANGLE_FORCE_GROUP)
        angle = state.getPotentialEnergy().value_in_unit(nrg)
        # Now dihedrals
        state = context.getState(getEnergy=True, enforcePeriodicBox=has_pbc,
                                 groups=2**parm.DIHEDRAL_FORCE_GROUP)
        dihedral = state.getPotentialEnergy().value_in_unit(nrg)
        # Now the CHARMM-specific terms
        if parm.chamber:
            # Now Urey-Bradley energy
            state = context.getState(getEnergy=True, enforcePeriodicBox=has_pbc,
                                     groups=2**parm.UREY_BRADLEY_FORCE_GROUP)
            ub = state.getPotentialEnergy().value_in_unit(nrg)
            # Now improper
            state = context.getState(getEnergy=True, enforcePeriodicBox=has_pbc,
                                     groups=2**parm.IMPROPER_FORCE_GROUP)
            imp = state.getPotentialEnergy().value_in_unit(nrg)
            # Now cmap
            if parm.has_cmap:
                state = context.getState(getEnergy=True, enforcePeriodicBox=has_pbc,
                                         groups=2**parm.CMAP_FORCE_GROUP)
                cmap = state.getPotentialEnergy().value_in_unit(nrg)
            else:
                cmap = 0

        # Now non-bonded. No real way to decompose this.
        state = context.getState(getEnergy=True, enforcePeriodicBox=has_pbc,
                                 groups=2**parm.NONBONDED_FORCE_GROUP)
        nonbond = state.getPotentialEnergy().value_in_unit(nrg)

        if parm.chamber:
            output.write('Bond         = %20.7f     Angle        = %20.7f\n'
                         'Dihedral     = %20.7f     Urey-Bradley = %20.7f\n'
                         'Improper     = %20.7f     ' % (bond, angle, dihedral, ub, imp))
            if parm.has_cmap:
                output.write('CMAP         = %20.7f\n' % cmap)
            output.write('Nonbond      = %20.7f\nTOTAL        = %20.7f\n' % (nonbond,
                         bond+angle+dihedral+ub+imp+cmap+nonbond))
        else:
            output.write('Bond     = %20.7f     Angle    = %20.7f\n'
                         'Dihedral = %20.7f     Nonbond  = %20.7f\n'
                         'TOTAL    = %20.7f\n' % (bond, angle, dihedral, nonbond,
                         (bond+angle+dihedral+nonbond)))
    else:
        state = context.getState(getEnergy=True, enforcePeriodicBox=parm.ptr('ifbox') > 0)

        output.write('Potential Energy = %.7f kcal/mol\n' %
                     state.getPotentialEnergy().value_in_unit(nrg))

def minimize(parm, igb, saltcon, cutoff, restraintmask, weight,
             script, platform, precision, norun, tol, maxcyc):
    """ Minimizes a snapshot. Use the existing System if it exists """
    global HAS_OPENMM
    if not HAS_OPENMM:
        raise SimulationError('Could not import OpenMM!')

    # Open the script file
    if script is not None:
        scriptfile = open(script, 'w')
        scriptfile.write(_SCRIPT_HEADER % (parm, parm.crdname))
    else:
        scriptfile = None

    gbmeth, kappa = None, 0.0
    if parm.ptr('ifbox') == 0:
        if cutoff is None or cutoff >= 500:
            nbmeth = ff.NoCutoff
            cutoff = 1000.0
        else:
            nbmeth = ff.CutoffNonPeriodic
        if not igb in (0, 1, 2, 5, 6, 7, 8):
            raise SimulationError('Bad igb value. Must be 0, 1, 2, 5, 6, 7, or 8')
        if igb == 1:
            gbmeth = HCT
        elif igb == 2:
            gbmeth = OBC1
        elif igb == 5:
            gbmeth = OBC2
        elif igb == 7:
            gbmeth = GBn
        elif igb == 8:
            gbmeth = GBn2
        # Other choices are vacuum
        kappa = 0.73 * sqrt(saltcon * 0.10806)
    else:
        if cutoff is None: cutoff = 8.0
        nbmeth = ff.PME

    # Time to create the OpenMM system
    system = parm.createSystem(nonbondedMethod=nbmeth, nonbondedCutoff=cutoff*u.angstrom,
                               constraints=None, rigidWater=True, removeCMMotion=False,
                               implicitSolvent=gbmeth, implicitSolventKappa=kappa*(1.0/u.angstrom),
                               soluteDielectric=1.0, solventDielectric=78.5,
                               ewaldErrorTolerance=5e-5)

    if script is not None:
        scriptfile.write('# Create the system\n'
                         'system = parm.createSystem(nonbondedMethod=%s,\n'
                         '                   nonbondedCutoff=%s*u.angstrom,\n'
                         '                   constraints=None, rigidWater=True,\n'
                         '                   removeCMMotion=False, implicitSolvent=%s,\n'
                         '                   implicitSolventKappa=%s*(1.0/u.angstrom),\n'
                         '                   soluteDielectric=1.0,\n'
                         '                   solventDielectric=78.5,\n'
                         '                   ewaldErrorTolerance=5e-5\n'
                         ')\n\n'
                         '# Create a dummy integrator\n'
                         'integrator = VerletIntegrator(2.0*u.femtoseconds)\n\n' %
                         (nbmeth, cutoff, gbmeth, kappa))

    # See if we need to add restraints
    if restraintmask is not None:
        system.addForce(
            positional_restraints(restraintmask, u.kilocalorie_per_mole/u.angstrom/u.angstrom,
                                  parm.coordinates, scriptfile=scriptfile,)
        )

    if script is not None:
        scriptfile.write('system.addForce(frc)\n\n')

    # Create a dummy integrator and the simulation
    integrator = mm.VerletIntegrator(2.0*u.femtoseconds)
    if platform is not None:
        plat = mm.Platform.getPlatformByName(platform)
        simulation = Simulation(parm.topology, system, integrator, plat)
        if script is not None:
            scriptfile.write('# Create the platform\n'
                             'plat = Platform.getPlatformByName(%s)\n\n'
                             '# Create the simulation\n'
                             'simulation = Simulation(parm.topology, system, '
                             'integrator, plat)\n\n' % (platform))
    else:
        simulation = Simulation(parm.topology, system, integrator)
        if script is not None:
            scriptfile.write('# Create the simulation\n'
                             'simulation = Simulation(parm.topology, system, integrator)\n\n')

    # Now assign the coordinates
    simulation.context.setPositions(parm.positions)
    if script is not None:
        scriptfile.write('# Setting the positions\n'
                         'simulation.context.setPositions(parm.positions)\n\n'
                         '# Minimize the energy\n'
                         'simulation.minimizeEnergy(tolerance=%s, maxIterations=%s)\n'
                         '# Store the positions back in the parmtop\n'
                         'state = simulation.context.getState(getPositions=True\n'
                         '                                    enforcePeriodicBox=%s)\n'
                         'parm.positions = state.getPositions()\n' %
                         (bool(parm.ptr('ifbox') > 0), tol, maxcyc or 0))
        scriptfile.close()

    # Go ahead and minimize now and set the coordinates from the results of this
    # minimization if we wanted to run the calculation
    if not norun:
        simulation.minimizeEnergy(tolerance=tol, maxIterations=maxcyc or 0)
        # Now get the coordinates
        state = simulation.context.getState(getPositions=True,
                                            enforcePeriodicBox=parm.ptr('ifbox') > 0)
        parm.positions = state.getPositions()
