import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u

import numpy as np

from esmacs import restrain


class Equilibrate:

    _TIMESTEP = 2 * u.femtosecond
    _FRICTION_COEFFICIENT = 5 / u.picosecond

    def __init__(self, system, positions, topology, num_steps=100):
        """Simulation to equilibrate the system

        :param system: simtk.openmm.System
        :param positions: numpy.ndarray
        :param topology: simtk.openmm.Topology
        :param num_steps: int
        """
        self.system = system
        self.positions = positions
        self.topology = topology
        self.num_steps = num_steps

    def equilibrate(self, pressure, temperature, velocities):

        K = 4 * u.kilocalorie / (u.angstrom ** 2 * u.mole)
        self.positions = restrain.restrain_atoms_by_dsl(self.system, pressure=True, positions=self.positions,
                                                        constant=K, topology=self.topology,
                                                        atoms_dsl='protein and not type H')

        barostat = mm.MonteCarloBarostat(pressure, temperature)
        self.system.addForce(barostat)

        integrator = mm.LangevinIntegrator(temperature, self._FRICTION_COEFFICIENT, self._TIMESTEP)
        simulation = app.Simulation(self.topology, self.system, integrator)

        simulation.context.setPositions(self.positions)
        simulation.context.setVelocities(velocities)

        print('Equilibrating with harmonic restraint constants:')
        for scaled_k in K * 10 * np.append(np.logspace(0, 10, num=11, base=0.5), 0):
            print(scaled_k)
            simulation.context.setParameter('K', scaled_k)
            simulation.step(self.num_steps)

        simulation.step(self.num_steps)
