import numpy as np

import simtk.unit as u
import simtk.openmm as mm
import simtk.openmm.app as app

from openmmtools.forcefactories import restrain_atoms_by_dsl
from openmmtools.states import ThermodynamicState, SamplerState

from openmmtools.utils import with_timer


class Esmacs:

    _TIMESTEP = 2 * u.femtosecond
    _FRICTION_COEFFICIENT = 5 / u.picosecond
    _K = 4.0 * u.kilocalorie / (u.angstrom ** 2 * u.mole)

    def __init__(self, thermodynamic_state: ThermodynamicState, sampler_state: SamplerState, topology):

        self.topology = topology
        self.sampler_state = sampler_state
        self.thermodynamic_state = thermodynamic_state

    @classmethod
    def from_amber(cls, prmtop, inpcrd, temperature=50*u.kelvin):
        prmtop = app.AmberPrmtopFile(prmtop)
        inpcrd = app.AmberInpcrdFile(inpcrd)
        system = prmtop.createSystem(nonbondedMethod=app.PME,
                                     constraints=app.HBonds,
                                     nonbondedCutoff=10*u.angstroms,
                                     switchDistance=8*u.angstroms)

        thermodynamic_state = ThermodynamicState(system, temperature=temperature)
        sampler_state = SamplerState(positions=inpcrd.getPositions(asNumpy=True), box_vectors=inpcrd.boxVectors)

        return Esmacs(thermodynamic_state, sampler_state, prmtop.topology)

    def apply_restraint(self):
        restrain_atoms_by_dsl(self.thermodynamic_state, self.sampler_state, self.topology,
                              atoms_dsl="protein and not type H")

    def remove_restraint(self):
        # Need to remove the restraining force because the positions of atoms have changed.
        system = self.thermodynamic_state.system
        forces = system.getForces()

        for index, force in enumerate(forces):
            if force.__class__.__name__ == 'CustomExternalForce':
                system.removeForce(index)
                break

        self.thermodynamic_state.system = system

    def get_context(self):
        integrator = mm.LangevinIntegrator(self.thermodynamic_state.temperature,
                                           self._FRICTION_COEFFICIENT, self._TIMESTEP)

        context = self.thermodynamic_state.create_context(integrator)

        self.sampler_state.apply_to_context(context)

        return context, integrator

    @with_timer
    def minimize_energy(self, max_iterations=100):
        self.apply_restraint()
        context, integrator = self.get_context()

        for K in self._K * 10 * np.append(np.logspace(0, 10, num=11, base=0.5), 0):
            context.setParameter('K', K)
            print('Minimizing for {} with K={}.'.format(max_iterations, K))
            mm.LocalEnergyMinimizer.minimize(context, maxIterations=max_iterations)

        self.sampler_state.update_from_context(context)

        self.remove_restraint()

        del context

    @with_timer
    def heat_system(self, destination_temperature=300*u.kelvin, num_steps=100, equilibrate=100):
        self.apply_restraint()
        context, integrator = self.get_context()

        context.setParameter('K', self._K)

        while self.thermodynamic_state.temperature < destination_temperature:
            self.thermodynamic_state.temperature += 1 * u.kelvin
            self.thermodynamic_state.apply_to_context(context)

            print('Heating to {}.'.format(self.thermodynamic_state.temperature))
            integrator.step(num_steps)

        print('Further NVT equilibration for {} at {}'.format(equilibrate, self.thermodynamic_state.temperature))
        integrator.step(equilibrate)

        self.sampler_state.update_from_context(context)
        self.remove_restraint()

        del context

    @with_timer
    def equilibrate_system(self, pressure=1*u.atmosphere, num_steps=100):

        self.thermodynamic_state.pressure = pressure

        self.apply_restraint()
        context, integrator = self.get_context()

        for K in self._K * 10 * np.logspace(0, 10, num=11, base=0.5):
            context.setParameter('K', K)
            print('NPT equilibration for {} with K={}.'.format(num_steps, K))
            integrator.step(num_steps)

        self.sampler_state.update_from_context(context)
        self.remove_restraint()

        del context

    @with_timer
    def simulate_system(self, equilibrate=100, production=100, dcd_frequency=10):

        context, integrator = self.get_context()

        # integrator.step(num_steps)

        dummy_integrator = mm.VerletIntegrator(1*u.femtosecond)
        simulation = app.Simulation(self.topology, context.getSystem(), dummy_integrator)
        simulation.context = context
        simulation.integrator = integrator

        simulation.reporters.append(app.DCDReporter('equilibration.dcd', dcd_frequency))
        simulation.step(equilibrate)

        simulation.reporters.clear()
        simulation.reporters.append(app.DCDReporter('production.dcd', dcd_frequency))
        simulation.step(production)

        self.sampler_state.update_from_context(context)

        del context

    @with_timer
    def run_protocol(self, short_run=False):
        if short_run:
            steps = [10, 100, 100, 100, 100, 100]
        else:
            steps = [1000, 100, 5000, 20000, 800000, 2000000]

        self.minimize_energy(max_iterations=steps[0])
        self.heat_system(destination_temperature=300*u.kelvin, num_steps=steps[1], equilibrate=steps[2])
        self.equilibrate_system(pressure=1*u.atmosphere, num_steps=steps[3])
        self.simulate_system(equilibrate=steps[4], production=steps[5], dcd_frequency=5_000)

    def run_raw_protocol(self):

        integrator = mm.LangevinIntegrator(50*u.kelvin, self._FRICTION_COEFFICIENT, self._TIMESTEP)

        simulation = app.Simulation(self.topology, self.system, integrator)
