import simtk.unit as u
import simtk.openmm as mm
import simtk.openmm.app as app

from openmmtools.forcefactories import restrain_atoms_by_dsl
from openmmtools.states import ThermodynamicState, SamplerState


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
                                     nonbondedCutoff=12 * u.angstroms,
                                     switchDistance=10 * u.angstroms)

        thermodynamic_state = ThermodynamicState(system, temperature=temperature)
        sampler_state = SamplerState(positions=inpcrd.getPositions(asNumpy=True), box_vectors=inpcrd.boxVectors)

        return Esmacs(thermodynamic_state, sampler_state, prmtop.topology)

    def apply_restraint(self):
        restrain_atoms_by_dsl(self.thermodynamic_state, self.sampler_state, self.topology,
                              atoms_dsl="protein and not type H")

    def remove_restraint(self):
        # Need to remove the restraining force because the positions of atoms have changed.
        system = self.thermodynamic_state.system
        system.removeForce(system.getNumForces()-1)
        self.thermodynamic_state.system = system

    def get_context(self):
        integrator = mm.LangevinIntegrator(self.thermodynamic_state.temperature,
                                           self._FRICTION_COEFFICIENT, self._TIMESTEP)

        context = self.thermodynamic_state.create_context(integrator)

        self.sampler_state.apply_to_context(context)

        return context, integrator

    def minimize_energy(self, max_iterations=1000):
        self.apply_restraint()

        context, integrator = self.get_context()

        restrain_scaling = 10.0

        for _ in range(10):
            context.setParameter('K', self._K * restrain_scaling)
            mm.LocalEnergyMinimizer.minimize(context, maxIterations=max_iterations)
            restrain_scaling *= 0.5

        context.setParameter('K', self._K * 0)
        mm.LocalEnergyMinimizer.minimize(context, maxIterations=max_iterations)

        self.sampler_state.update_from_context(context)
        self.remove_restraint()

    def heat_system(self, destination_temperature=300*u.kelvin, num_steps=100, equilibrate=5000):
        self.apply_restraint()

        context, integrator = self.get_context()

        context.setParameter('K', self._K)

        while self.thermodynamic_state.temperature <= destination_temperature:
            self.thermodynamic_state.temperature += 1 * u.kelvin
            self.thermodynamic_state.apply_to_context(context)
            integrator.step(num_steps)

        integrator.step(equilibrate)

        self.sampler_state.update_from_context(context)

        self.remove_restraint()

    def equilibrate_system(self, pressure=1*u.atmosphere, num_steps=100, equilibrate=5000):

        self.thermodynamic_state.pressure = pressure

        self.apply_restraint()

        context, integrator = self.get_context()

        restrain_scaling = 10.0

        for _ in range(10):
            context.setParameter('K', self._K * restrain_scaling)
            integrator.step(num_steps)
            restrain_scaling *= 0.5

        context.setParameter('K', self._K * 0)
        integrator.step(equilibrate)

        self.sampler_state.update_from_context(context)

        self.remove_restraint()

    def simulate_system(self, num_steps=100):

        context, integrator = self.get_context()

        integrator.step(num_steps)
        self.sampler_state.update_from_context(context)
