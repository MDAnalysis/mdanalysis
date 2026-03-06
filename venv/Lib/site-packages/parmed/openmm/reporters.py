"""
This module contains a handful of extra reporter classes for OpenMM
simulations
"""
from math import isnan, isinf
from time import time

from ..amber.asciicrd import AmberMdcrd
from ..geometry import box_vectors_to_lengths_and_angles
from ..amber.netcdffiles import NetCDFTraj
from ..amber.readparm import Rst7
from .. import unit as u
from ..utils.decorators import needs_openmm
VELUNIT = u.angstrom / u.picosecond
FRCUNIT = u.kilocalorie_per_mole / u.angstrom

class StateDataReporter:
    """
    This class acts as a state data reporter for OpenMM simulations, but it is a
    little more generalized. Notable differences are:

      -  It allows the units of the output to be specified, with defaults being
         those used in Amber (e.g., kcal/mol, angstroms^3, etc.)
      -  It will write to any object containing a 'write' method; not just
         files. This allows, for instance, writing to a GUI window that
         implements the desired 'write' attribute.

    Most of this code is copied from the OpenMM StateDataReporter class, with
    the above-mentioned changes made.

    Parameters
    ----------
    f : str or file-like
        Destination to write the state data (file name or file object)
    reportInterval : int
        Number of steps between state data reports
    step : bool, optional
        Print out the step number (Default True)
    time : bool, optional
        Print out the simulation time (Defaults True)
    potentialEnergy : bool, optional
        Print out the potential energy of the structure (Default True)
    kineticEnergy : bool, optional
        Print out the kinetic energy of the structure (Default True)
    totalEnergy : bool, optional
        Print out the total energy of the system (Default True)
    temperature : bool, optional
        Print out the temperature of the system (Default True)
    volume : bool, optional
        Print out the volume of the unit cell. If the system is not periodic,
        the value is meaningless (Default False)
    density : bool, optional
        Print out the density of the unit cell. If the system is not periodic,
        the value is meaningless (Default False)
    separator : str, optional
        The string to separate data fields (Default ',')
    systemMass : float, optional
        If not None, the density will be computed from this mass, since setting
        a mass to 0 is used to constrain the position of that particle. (Default
        None)
    energyUnit : unit, optional
        The units to print energies in (default unit.kilocalories_per_mole)
    timeUnit : unit, optional
        The units to print time in (default unit.picoseconds)
    volumeUnit : unit, optional
        The units print volume in (default unit.angstroms**3)
    densityUnit : unit, optional
        The units to print density in (default
        unit.grams/unit.item/unit.milliliter)
    """

    @needs_openmm
    def __init__(self, f, reportInterval, step=True, time=True,
                 potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                 temperature=True, volume=False, density=False, separator=',',
                 systemMass=None, energyUnit=u.kilocalories_per_mole,
                 timeUnit=u.picoseconds, volumeUnit=u.angstroms**3,
                 densityUnit=u.grams/u.item/u.milliliter):
        self._reportInterval = reportInterval
        self._openedFile = not hasattr(f, 'write')
        if self._openedFile:
            self._out = open(f, 'w')
        else:
            self._out = f
        self._step = step
        self._time = time
        self._potentialEnergy = potentialEnergy
        self._kineticEnergy = kineticEnergy
        self._totalEnergy = totalEnergy
        self._temperature = temperature
        self._volume = volume
        self._density = density
        self._separator = separator
        self._totalMass = systemMass
        self._hasInitialized = False
        self._needsPositions = False
        self._needsVelocities = False
        self._needsForces = False
        self._needEnergy = (potentialEnergy or kineticEnergy or totalEnergy or temperature)
        self._energyUnit = energyUnit
        self._densityUnit = densityUnit
        self._timeUnit = timeUnit
        self._volumeUnit = volumeUnit

    def describeNextReport(self, simulation):
        """
        Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : :class:`app.Simulation`
            The simulation to generate a report for

        Returns
        -------
        nsteps, pos, vel, frc, ene : int, bool, bool, bool, bool
            nsteps is the number of steps until the next report
            pos, vel, frc, and ene are flags indicating whether positions,
            velocities, forces, and/or energies are needed from the Context
        """
        steps_left = simulation.currentStep % self._reportInterval
        steps = self._reportInterval - steps_left
        return (steps, self._needsPositions, self._needsVelocities, self._needsForces, self._needEnergy)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : :class:`app.Simulation`
            The Simulation to generate a report for
        state : :class:`mm.State`
            The current state of the simulation
        """
        if not self._hasInitialized:
            self._initializeConstants(simulation)
            headers = self._constructHeaders()
            header = f'"{self._separator}"'.join(headers)
            self._out.write(f'#"{header}"\n')
            self._hasInitialized = True

        # Check for errors.
        self._checkForErrors(simulation, state)

        # Query for the values
        values = self._constructReportValues(simulation, state)

        # Write the values.
        self._out.write(self._separator.join(str(v) for v in values) + '\n')
        if hasattr(self._out, 'flush'):
            self._out.flush()

    def _constructReportValues(self, simulation, state):
        """
        Query the simulation for the current state of our observables of
        interest.

        Parameters
        ----------
        simulation : Simulation
            The Simulation object to generate a report for
        state : State
            The current state of the simulation object

        Returns: A list of values summarizing the current state of
        the simulation, to be printed or saved. Each element in the list
        corresponds to one of the columns in the resulting CSV file.
        """
        values = []
        if self._volume or self._density:
            volume = state.getPeriodicBoxVolume()
        if self._density:
            density = self._totalMass / volume
        ke = state.getKineticEnergy()
        pe = state.getPotentialEnergy()
        if self._temperature:
            temp = 2 * ke / (self._dof * u.MOLAR_GAS_CONSTANT_R)
        time = state.getTime()
        if self._step:
            values.append(simulation.currentStep)
        if self._time:
            values.append(time.value_in_unit(self._timeUnit))
        if self._potentialEnergy:
            values.append(pe.value_in_unit(self._energyUnit))
        if self._kineticEnergy:
            values.append(ke.value_in_unit(self._energyUnit))
        if self._totalEnergy:
            values.append((pe + ke).value_in_unit(self._energyUnit))
        if self._temperature:
            values.append(temp.value_in_unit(u.kelvin))
        if self._volume:
            values.append(volume.value_in_unit(self._volumeUnit))
        if self._density:
            values.append(density.value_in_unit(self._densityUnit))
        return values

    def _initializeConstants(self, simulation):
        """
        Initialize a set of constants required for the reports

        Parameters
        ----------
        simulation : Simulation
            The simulation to generate a report for
        """
        import openmm as mm
        system = simulation.system
        frclist = system.getForces()
        if self._temperature:
            # Compute the number of degrees of freedom.
            dof = 0
            for i in range(system.getNumParticles()):
                if system.getParticleMass(i) > 0*u.dalton:
                    dof += 3
            dof -= system.getNumConstraints()
            if any(isinstance(frc, mm.CMMotionRemover) for frc in frclist):
                dof -= 3
            self._dof = dof
        if self._density:
            if self._totalMass is None:
                # Compute the total system mass.
                self._totalMass = 0*u.dalton
                for i in range(system.getNumParticles()):
                    self._totalMass += system.getParticleMass(i)
            elif not u.is_quantity(self._totalMass):
                self._totalMass = self._totalMass*u.dalton

    def _constructHeaders(self):
        """Construct the headers for the CSV output

        Returns: a list of strings giving the title of each observable being
                 reported on.
        """
        headers = []
        if self._step:
            headers.append('Step')
        if self._time:
            headers.append('Time (ps)')
        if self._potentialEnergy:
            headers.append(f'Potential Energy ({self._energyUnit})')
        if self._kineticEnergy:
            headers.append(f'Kinetic Energy ({self._energyUnit})')
        if self._totalEnergy:
            headers.append(f'Total Energy ({self._energyUnit})')
        if self._temperature:
            headers.append('Temperature (K)')
        if self._volume:
            headers.append(f'Box Volume ({self._volumeUnit})')
        if self._density:
            headers.append(f'Density ({self._densityUnit})')
        return headers

    def _checkForErrors(self, simulation, state):
        """Check for errors in the current state of the simulation

        Parameters
        ----------
        simulation : :class:`app.Simulation`
            The Simulation to generate a report for
        state : :class:`State`
            The current state of the simulation
        """
        if self._needEnergy:
            energy = state.getKineticEnergy() + state.getPotentialEnergy()
            if isnan(energy._value):
                raise ValueError('Energy is NaN') # pragma: no cover
            if isinf(energy._value):
                raise ValueError('Energy is infinite') # pragma: no cover

    def __del__(self):
        if hasattr(self, '_openedFile') and self._openedFile:
            self._out.close()

    def finalize(self):
        """ Closes any open file """
        try:
            if self._out is not None and self._openedFile:
                self._out.close()
        except AttributeError: # pragma: no cover
            pass               # pragma: no cover

class NetCDFReporter(object):
    """ NetCDFReporter prints a trajectory in NetCDF format """

    @needs_openmm
    def __init__(self, file, reportInterval, crds=True, vels=False, frcs=False):
        """
        Create a NetCDFReporter instance.

        Parameters
        ----------
        file : str
            Name of the file to write the trajectory to
        reportInterval : int
            How frequently to write a frame to the trajectory
        crds : bool=True
            Should we write coordinates to this trajectory? (Default True)
        vels : bool=False
            Should we write velocities to this trajectory? (Default False)
        frcs : bool=False
            Should we write forces to this trajectory? (Default False)
        """
        if not crds and not vels and not frcs:
            raise ValueError('Either coordinates, velocities, or forces needed for NetCDFReporter')
        # Control flags
        self.crds, self.vels, self.frcs = crds, vels, frcs
        self._reportInterval = reportInterval
        self._out = None # not written yet
        self.fname = file

    def describeNextReport(self, simulation):
        """
        Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : :class:`app.Simulation`
            The Simulation to generate a report for

        Returns
        -------
        nsteps, pos, vel, frc, ene : int, bool, bool, bool, bool
            nsteps is the number of steps until the next report
            pos, vel, frc, and ene are flags indicating whether positions,
            velocities, forces, and/or energies are needed from the Context
        """
        stepsleft = simulation.currentStep % self._reportInterval
        steps = self._reportInterval - stepsleft
        return (steps, self.crds, self.vels, self.frcs, False)


    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : :class:`app.Simulation`
            The Simulation to generate a report for
        state : :class:`mm.State`
            The current state of the simulation
        """
        global VELUNIT, FRCUNIT
        if self.crds:
            crds = state.getPositions().value_in_unit(u.angstrom)
        if self.vels:
            vels = state.getVelocities().value_in_unit(VELUNIT)
        if self.frcs:
            frcs = state.getForces().value_in_unit(FRCUNIT)
        if self._out is None:
            # This must be the first frame, so set up the trajectory now
            if self.crds:
                atom = len(crds)
            elif self.vels:
                atom = len(vels)
            elif self.frcs:
                atom = len(frcs)
            self.uses_pbc = simulation.topology.getUnitCellDimensions() is not None
            self._out = NetCDFTraj.open_new(
                self.fname, atom, self.uses_pbc, self.crds, self.vels,
                self.frcs, title="ParmEd-created trajectory using OpenMM",
            )
        if self.uses_pbc:
            vecs = state.getPeriodicBoxVectors()
            lengths, angles = box_vectors_to_lengths_and_angles(*vecs)
            self._out.add_cell_lengths_angles(
                lengths.value_in_unit(u.angstrom), angles.value_in_unit(u.degree)
            )

        # Add the coordinates, velocities, and/or forces as needed
        if self.crds:
            self._out.add_coordinates(crds)
        if self.vels:
            # The velocities get scaled right before writing
            self._out.add_velocities(vels)
        if self.frcs:
            self._out.add_forces(frcs)
        # Now it's time to add the time.
        self._out.add_time(state.getTime().value_in_unit(u.picosecond))

    def __del__(self):
        try:
            if self._out is not None:
                self._out.close()
        except AttributeError:
            pass

    def finalize(self):
        """ Closes any open file """
        try:
            if self._out is not None:
                self._out.close()
        except AttributeError: # pragma: no cover
            pass               # pragma: no cover

class MdcrdReporter(object):
    """
    MdcrdReporter prints a trajectory in ASCII Amber format. This reporter will
    be significantly slower than binary file reporters (like DCDReporter or
    NetCDFReporter).

    Parameters
    ----------
    file : str
        Name of the file to write the trajectory to
    reportInterval : int
        Number of steps between writing trajectory frames
    crds : bool=True
        Write coordinates to this trajectory file?
    vels : bool=False
        Write velocities to this trajectory file?
    frcs : bool=False
        Write forces to this trajectory file?

    Notes
    -----
    You can only write one of coordinates, forces, or velocities to a mdcrd
    file.
    """

    @needs_openmm
    def __init__(self, file, reportInterval, crds=True, vels=False, frcs=False):
        # ASCII mdcrds can have either coordinates, forces, or velocities
        ntrue = 0
        if crds: ntrue += 1
        if vels: ntrue += 1
        if frcs: ntrue += 1
        if ntrue != 1:
            raise ValueError(
                'MdcrdReporter must print exactly one of either coordinates, velocities, or forces.'
            )
        # Control flags
        self.crds, self.vels, self.frcs = crds, vels, frcs
        self._reportInterval = reportInterval
        self._out = None # not written yet
        self.fname = file

    def describeNextReport(self, simulation):
        """
        Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : :class:`app.Simulation`
            The Simulation to generate a report for

        Returns
        -------
        nsteps, pos, vel, frc, ene : int, bool, bool, bool, bool
            nsteps is the number of steps until the next report
            pos, vel, frc, and ene are flags indicating whether positions,
            velocities, forces, and/or energies are needed from the Context
        """
        stepsleft = simulation.currentStep % self._reportInterval
        steps = self._reportInterval - stepsleft
        return (steps, self.crds, self.vels, self.frcs, False)

    def report(self, simulation, state):
        """
        Generate a report.

        Parameters:
            - simulation (Simulation) The Simulation to generate a report for
            - state (State) The current state of the simulation
        """
        from ..amber.asciicrd import VELSCALE
        global VELUNIT, FRCUNIT
        if self.crds:
            crds = state.getPositions().value_in_unit(u.angstrom)
        elif self.vels: # crds/vels/frcs are exclusive, elif works
            vels = state.getVelocities().value_in_unit(VELUNIT)
        elif self.frcs:
            frcs = state.getForces().value_in_unit(FRCUNIT)
        if self._out is None:
            # This must be the first frame, so set up the trajectory now
            if self.crds:
                self.atom = len(crds)
            elif self.vels:
                self.atom = len(vels)
            elif self.frcs:
                self.atom = len(frcs)
            self.uses_pbc = simulation.topology.getUnitCellDimensions() is not None
            self._out = AmberMdcrd(self.fname, self.atom, self.uses_pbc,
                    title="ParmEd-created trajectory using OpenMM", mode='w')

        # Add the coordinates, velocities, and/or forces as needed
        if self.crds:
            flatcrd = [0 for i in range(self.atom*3)]
            for i in range(self.atom):
                i3 = i*3
                flatcrd[i3], flatcrd[i3+1], flatcrd[i3+2] = crds[i]
            self._out.add_coordinates(flatcrd)
        if self.vels:
            # Divide by the scaling factor (works if vels is a list of Vec3's)
            # This is necessary since AmberMdcrd does not scale before writing
            # (since it expects coordinates)
            vels = [v / VELSCALE for v in vels]
            flatvel = [0 for i in range(self.atom*3)]
            for i in range(self.atom):
                i3 = i*3
                flatvel[i3], flatvel[i3+1], flatvel[i3+2] = vels[i]
            self._out.add_coordinates(flatvel)
        if self.frcs:
            flatfrc = [0 for i in range(self.atom*3)]
            for i in range(self.atom):
                i3 = i*3
                flatfrc[i3], flatfrc[i3+1], flatfrc[i3+2] = frcs[i]
            self._out.add_coordinates(flatfrc)
        # Now it's time to add the box lengths
        if self.uses_pbc:
            boxvecs = state.getPeriodicBoxVectors()
            lengths, angles = box_vectors_to_lengths_and_angles(*boxvecs)
            self._out.add_box(lengths.value_in_unit(u.angstroms))

    def __del__(self):
        try:
            if self._out is not None:
                self._out.close()
        except AttributeError:
            pass

    def finalize(self):
        """ Closes any open file """
        try:
            if self._out is not None:
                self._out.close()
        except AttributeError: # pragma: no cover
            pass               # pragma: no cover

class RestartReporter(object):
    """
    Use a reporter to handle writing restarts at specified intervals.

    Parameters
    ----------
    file : str
        Name of the file to write the restart to.
    reportInterval : int
        Number of steps between writing restart files
    write_multiple : bool=False
        Either write a separate restart each time (appending the step number in
        the format .# to the file name given above) if True, or overwrite the
        same file each time if False
    netcdf : bool=False
        Use the Amber NetCDF restart file format
    write_velocities : bool=True
        Write velocities to the restart file. You can turn this off for passing
        in, for instance, a minimized structure.
    """

    @needs_openmm
    def __init__(self, file, reportInterval, write_multiple=False, netcdf=False,
                 write_velocities=True):
        self.fname = file
        self._reportInterval = reportInterval
        self.write_multiple = write_multiple
        self.netcdf = netcdf
        self.write_velocities = write_velocities
        self.rst7 = None

    def describeNextReport(self, simulation):
        """
        Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : :class:`app.Simulation`
            The Simulation to generate a report for

        Returns
        -------
        nsteps, pos, vel, frc, ene : int, bool, bool, bool, bool
            nsteps is the number of steps until the next report
            pos, vel, frc, and ene are flags indicating whether positions,
            velocities, forces, and/or energies are needed from the Context
        """
        stepsleft = simulation.currentStep % self._reportInterval
        steps = self._reportInterval - stepsleft
        return (steps, True, True, False, False)

    def report(self, sim, state):
        """Generate a report.

        Parameters
        ----------
        sim : :class:`app.Simulation`
            The Simulation to generate a report for
        state : :class:`mm.State`
            The current state of the simulation
        """
        global VELUNIT
        crds = state.getPositions().value_in_unit(u.angstrom)
        if self.rst7 is None:
            self.uses_pbc = sim.topology.getUnitCellDimensions() is not None
            self.atom = len(crds)
            # First time written
            self.rst7 = Rst7(natom=self.atom, title='Restart file written by ParmEd with OpenMM')
        self.rst7.time = state.getTime().value_in_unit(u.picosecond)
        flatcrd = [0.0 for i in range(self.atom*3)]
        for i in range(self.atom):
            i3 = i*3
            flatcrd[i3], flatcrd[i3+1], flatcrd[i3+2] = crds[i]
        self.rst7.coordinates = flatcrd

        if self.write_velocities:
            vels = state.getVelocities().value_in_unit(VELUNIT)
            flatvel = [0.0 for i in range(self.atom*3)]
            for i in range(self.atom):
                i3 = i*3
                flatvel[i3], flatvel[i3+1], flatvel[i3+2] = vels[i]
            self.rst7.vels = flatvel

        if self.uses_pbc:
            boxvecs = state.getPeriodicBoxVectors()
            lengths, angles = box_vectors_to_lengths_and_angles(*boxvecs)
            lengths = lengths.value_in_unit(u.angstrom)
            angles = angles.value_in_unit(u.degree)
            self.rst7.box = [lengths[0], lengths[1], lengths[2], angles[0], angles[1], angles[2]]

        if self.write_multiple:
            fname = self.fname + '.%d' % sim.currentStep
        else:
            fname = self.fname

        self.rst7.write(fname, self.netcdf)

    def finalize(self):
        """ No-op here """
        pass

class ProgressReporter(StateDataReporter):
    """
    A class that prints out a progress report of how much MD (or minimization)
    has been done, how fast the simulation is running, and how much time is left
    (similar to the mdinfo file in Amber)

    Parameters
    ----------
    f : str
        The file name of the progress report file (overwritten each time)
    reportInterval : int
        The step interval between which to write frames
    totalSteps : int
        The total number of steps that will be run in the simulation (used to
        estimate time remaining)
    potentialEnergy : bool, optional
        Whether to print the potential energy (default True)
    kineticEnergy : bool, optional
        Whether to print the kinetic energy (default True)
    totalEnergy : bool, optional
        Whether to print the total energy (default True)
    temperature : bool, optional
        Whether to print the system temperature (default True)
    volume : bool, optional
        Whether to print the system volume (default False)
    density : bool, optional
        Whether to print the system density (default False)
    systemMass : float or :class:`unit.Quantity`
        The mass of the system used when reporting density (useful in instances
        where masses are set to 0 to constrain their positions)

    Notes
    -----
    In addition to the above, :class:`ProgressReporter` also accepts arguments for
    :class:`StateDataReporter`
    """

    @needs_openmm
    def __init__(self, f, reportInterval, totalSteps, potentialEnergy=True,
                 kineticEnergy=True, totalEnergy=True, temperature=False,
                 volume=False, density=False, systemMass=None, **kwargs):
        # Make sure we got a file name rather than a file-like object.
        # Immediately close the file after opening. This erases it, which isn't
        # a bad thing, but also prepares it for reports
        kwargs['time'] = kwargs['step'] = True
        super(ProgressReporter, self).__init__(
            f, reportInterval, potentialEnergy=potentialEnergy,
            kineticEnergy=kineticEnergy, totalEnergy=totalEnergy,
            temperature=temperature, volume=volume, density=density,
            systemMass=systemMass, **kwargs
        )
        if not self._openedFile:
            raise ValueError('ProgressReporter requires a file name (not file object)')
        self._out.close()
        del self._out
        self.fname = f
        self._totalSteps = totalSteps

        # For timing
        self._startTime = None
        self._firstStep = None
        self._lastReportTime = None
        self._timeStep = None

    def describeNextReport(self, simulation):
        """
        Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : :class:`app.Simulation`
            The Simulation to generate a report for

        Returns
        -------
        nsteps, pos, vel, frc, ene : int, bool, bool, bool, bool
            nsteps is the number of steps until the next report
            pos, vel, frc, and ene are flags indicating whether positions,
            velocities, forces, and/or energies are needed from the Context
        """
        if self._startTime is None:
            # First time this is called, initialize the timers
            self._startTime = time()
            self._lastReportTime = time()
            self._timeStep = simulation.integrator.getStepSize()
            self._timeStep = self._timeStep.value_in_unit(u.nanosecond)
            self._firstStep = simulation.currentStep
        stepsleft = simulation.currentStep % self._reportInterval
        steps = self._reportInterval - stepsleft
        return (steps, False, False, False, self._needEnergy)

    def report(self, simulation, state):
        """
        Generate a report and predict the time to completion (and
        current/overall MD performance)
        """
        if not self._hasInitialized:
            self._initializeConstants(simulation)
            self._hasInitialized = True

        # Check for errors.
        self._checkForErrors(simulation, state)

        # Query for the values
        values = self._constructReportValues(simulation, state)

        now = time()

        total_time = now - self._startTime
        partial_time = now - self._lastReportTime
        self._lastReportTime = now

        total_nsperday = (values['step'] - self._firstStep) * self._timeStep / total_time
        partial_nsperday = (self._reportInterval*self._timeStep) / partial_time

        # Get total and partial ns/day (currently ns/second)
        total_nsperday *= 3600 * 24
        partial_nsperday *= 3600 * 24

        # Estimated time to completion (based on performance of last N steps
        remaining_steps = self._totalSteps - values['step'] + self._firstStep
        etc = partial_time / self._reportInterval * remaining_steps

        etc, unitstr = _format_time(etc)

        # Write the values.
        with open(self.fname, 'w') as f:
            f.write('-+' * 39 + '\n')
            f.write('\n')
            f.write(f'  On step {values["step"] - self._firstStep} of {self._totalSteps}\n')
            f.write('\n')
            if self._totalEnergy:
                f.write(f'  Total Energy     = {values["totalEnergy"]:12.4f}\n')
            if self._potentialEnergy:
                f.write(f'  Potential Energy = {values["potentialEnergy"]:12.4f}\n')
            if self._kineticEnergy:
                f.write(f'  Kinetic Energy   = {values["kineticEnergy"]:12.4f}\n')
            if self._volume:
                f.write(f'  Volume           = {values["volume"]:12.4f}\n')
            if self._density:
                f.write(f'  Density          = {values["density"]:12.4f}\n')
            if self._temperature:
                f.write(f'  Temperature      = {values["temperature"]:12.4f}\n')
            f.write('\n')
            f.write(
                f' Time for last {self._reportInterval:8d} steps: {partial_time:10.4f} s. ({partial_nsperday:.3f} ns/day)\n'
            )
            f.write(f' Time for all {values["step"] - self._firstStep:9d} steps: {total_time:10.4f} s. ({total_nsperday:.3f} ns/day)\n')
            f.write('\n')
            f.write(f' Estimated time to completion: {etc:.3f} {unitstr}\n')
            f.write('\n')
            f.write('-+' * 39 + '\n')

        if remaining_steps == 0:
            self._startTime = None

    def _constructReportValues(self, simulation, state):
        """
        Query the simulation for the current state of our observables of
        interest.

        Parameters:
            - simulation (Simulation) The Simulation to generate a report for
            - state (State) The current state of the simulation

        Returns: A list of values summarizing the current state of the
            simulation, to be printed or saved. Each element in the list
            corresponds to one of the columns in the resulting CSV file.
        """
        values = dict()
        values['step'] = simulation.currentStep
        values['time'] = state.getTime()
        volume = state.getPeriodicBoxVolume()
        pe = state.getPotentialEnergy().value_in_unit(self._energyUnit)
        ke = state.getKineticEnergy()
        if self._temperature:
            temp = 2 * ke / (self._dof * u.MOLAR_GAS_CONSTANT_R)
        ke = ke.value_in_unit(self._energyUnit)
        if self._potentialEnergy:
            values['potentialEnergy'] = pe
        if self._kineticEnergy:
            values['kineticEnergy'] = ke
        if self._totalEnergy:
            values['totalEnergy'] = pe + ke
        if self._temperature:
            values['temperature'] = temp.value_in_unit(u.kelvin)
        if self._volume:
            values['volume'] = volume.value_in_unit(self._volumeUnit)
        if self._density:
            dens = self._totalMass / volume
            # md_unit_system does the right thing that every unit system would
            # want
            values['density'] = dens.value_in_unit(self._densityUnit)

        return values

    def __del__(self):
        """ We already closed the file. """

class EnergyMinimizerReporter(StateDataReporter):
    """
    This class acts as a simple energy reporter class for minimizations. This
    is not meant to be used as a reporter for OpenMM's molecular dynamics
    routines, but instead passed a simulation object for printing out
    single-point energies.

    Parameters
    ----------
    f : str or file-like
        File name or object to write energies to
    volume : bool, optional
        If True, write the system volume (default is False)
    """

    def __init__(self, f, volume=False, **kwargs):
        super().__init__(f, 1, **kwargs)
        self._volume = volume

    def describeNextReport(self, *args, **kwargs):
        """ Disable this reporter inside MD """
        raise NotImplementedError('EnergyMinimizerReporter is not intended for '
                                  'use in reporting on molecular dynamics')

    def report(self, simulation, frame=None):
        """ Print out the current energy """
        has_pbc = simulation.topology.getUnitCellDimensions() is not None
        state = simulation.context.getState(getEnergy=True, enforcePeriodicBox=has_pbc)
        if frame is not None:
            self._out.write(f'Frame: {frame:10d}\n')
        e = state.getPotentialEnergy().value_in_unit(self._energyUnit)
        self._out.write(f'   Potential Energy = {e:12.4f} {self._energyUnit}\n')
        if has_pbc and (self._volume or self._density):
            vol = state.getPeriodicBoxVolume().value_in_unit(self._volumeUnit)
        if has_pbc and self._volume:
            self._out.write(f'   Volume = {vol:12.4f} {self._volumeUnit}\n')
        self._out.write('\n')

    def finalize(self):
        """ Closes any open file """
        try:
            if self._out is not None:
                self._out.close()
        except AttributeError: # pragma: no cover
            pass               # pragma: no cover

# Private helper functions
def _format_time(etc):
    """ Formats how time is printed in ProgressReporter """
    if etc > 3600:
        etc /= 3600
        unitstr = 'hr.'
    elif etc > 60:
        etc /= 60
        unitstr = 'min.'
    else:
        unitstr = 'sec.'
    return etc, unitstr
