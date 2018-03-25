import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u

import numpy as np

from esmacs import restrain


class Heater:

    _TIMESTEP = 2 * u.femtosecond
    _FRICTION_COEFFICIENT = 5 / u.picosecond

    def __init__(self, system, positions, topology, num_steps=100):

        self.system = system
        self.positions = positions
        self.topology = topology
        self.num_steps = num_steps

    def heat(self, t_i, t_f):

        K = 4 * u.kilocalorie / (u.angstrom ** 2 * u.mole)
        restrain.restrain_atoms_by_dsl(self.system, pressure=False, positions=self.positions,
                                       constant=K, topology=self.topology,
                                       atoms_dsl='protein and not type H')

        integrator = mm.LangevinIntegrator(t_i*u.kelvin, self._FRICTION_COEFFICIENT, self._TIMESTEP)
        simulation = app.Simulation(self.topology, self.system, integrator)

        simulation.context.setPositions(self.positions)
        simulation.context.setVelocitiesToTemperature(t_i)

        print('Heating system:')
        for temperature in np.arange(t_i, t_f, dtype=np.float)*u.kelvin:
            print(temperature)
            integrator.setTemperature(temperature)
            simulation.step(self.num_steps)

        integrator.setTemperature(t_f)
        simulation.step(5000)

        return simulation.context
