import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u

import numpy as np

from esmacs import restrain


class Minimizer:
    def __init__(self, system, positions, topology, num_steps=1000):

        self.system = system
        self.positions = positions
        self.topology = topology
        self.num_steps = num_steps

    def minimize(self):

        K = 4 * u.kilocalorie / (u.angstrom ** 2 * u.mole)
        restrain.restrain_atoms_by_dsl(self.system, pressure=False, positions=self.positions,
                                       constant=K, topology=self.topology,
                                       atoms_dsl='protein and not type H')

        integrator = mm.VerletIntegrator(1*u.femtosecond)
        simulation = app.Simulation(self.topology, self.system, integrator)

        simulation.context.setPositions(self.positions)
        simulation.context.setVelocitiesToTemperature(50*u.kelvin)

        print('Minimizing with harmonic restraint constants:')
        for scaled_k in K * 10 * np.append(np.logspace(0, 10, num=11, base=0.5), 0):
            print(scaled_k)
            simulation.context.setParameter('K', scaled_k)
            simulation.minimizeEnergy(maxIterations=self.num_steps)

        self.system.removeForce(self.system.getNumForces()-1)

        return simulation.context
