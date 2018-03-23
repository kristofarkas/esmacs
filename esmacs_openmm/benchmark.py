#!/bin/env python

"""
Benchmark.
"""

from simtk.openmm import app
from simtk import unit, openmm

from parmed.amber import AmberParm

import time


def main():

    prmtop_filename = '/lustre/atlas/scratch/farkaspall/chm126/inspire-data/nilotinib/e255k/build/e255k-complex.top'
    crd_filename = '/lustre/atlas/scratch/farkaspall/chm126/inspire-data/nilotinib/e255k/build/e255k-complex.inpcrd'

    prmtop = AmberParm(prmtop_filename, crd_filename)

    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=10*unit.angstrom, constraints=app.HBonds)

    temperature = 300 * unit.kelvin
    collision_rate = 5 / unit.picosecond
    timestep = 2 * unit.femtosecond

    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)

    system.setDefaultPeriodicBoxVectors(*prmtop.box_vectors)

    context = openmm.Context(system, integrator)

    context.setPositions(prmtop.positions)
    context.setVelocitiesToTemperature(temperature)

    print('System contains {} atoms.'.format(system.getNumParticles()))
    print('Using platform "{}".'.format(context.getPlatform().getName()))
    print('Initial potential energy is {}'.format(context.getState(getEnergy=True).getPotentialEnergy()))

    # Warm up the integrator to compile kernels, etc
    print('Warming up integrator to trigger kernel compilation...')
    integrator.step(10)

    # Time integration
    print('Benchmarking...')
    nsteps = 5000
    initial_time = time.time()
    integrator.step(nsteps)
    final_time = time.time()
    elapsed_time = (final_time - initial_time) * unit.seconds
    simulated_time = nsteps * timestep
    performance = (simulated_time / elapsed_time)
    print('Completed {} steps in {}.'.format(nsteps, elapsed_time))
    print('Performance is {} ns/day'.format(performance / (unit.nanoseconds/unit.day)))
    print('Final potential energy is {}'.format(context.getState(getEnergy=True).getPotentialEnergy()))


if __name__ == '__main__':
    main()
