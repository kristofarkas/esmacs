import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u

import numpy as np
import mdtraj
import time
import os

from esmacs import restrain

# Parameters

Ti = 50 * u.kelvin
T = 300 * u.kelvin
P = 1 * u.atmosphere

ts = 2 * u.femtosecond

kT = T * u.BOLTZMANN_CONSTANT_kB * u.AVOGADRO_CONSTANT_NA
sigma = 3 * u.angstrom

# System

system_path = '/home/kristof/Research/INSPIRE/nilotinib/e255k/build'
prmtop = app.AmberPrmtopFile(os.path.join(system_path, 'complex.top'))
inpcrd = app.AmberInpcrdFile(os.path.join(system_path, 'complex.inpcrd'))
pdb = app.PDBFile(os.path.join(system_path, 'complex.pdb'))

system = prmtop.createSystem(nonbondedMethod=app.PME,
                             constraints=app.HBonds,
                             nonbondedCutoff=12*u.angstroms,
                             switchDistance=10*u.angstroms)
topology = mdtraj.Topology.from_openmm(prmtop.topology)

K = 4 * u.kilocalorie/(u.angstrom**2*u.mole)  # kT / sigma**2

restrain.restrain_atoms_by_dsl(system,
                               pressure=False,
                               positions=pdb.positions,
                               atoms_dsl='not water and not type H',
                               topology=topology,
                               constant=K)


# Initialise simulation and positions

integrator = mm.LangevinIntegrator(Ti, 5/u.picoseconds, ts)

simulation = app.Simulation(prmtop.topology, system, integrator)

simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(Ti)

# Minimise

print('Minimizing, from initial energy:', simulation.context.getState(getEnergy=True).getPotentialEnergy())

initial_time = time.time()

for scaled_k in K*10*np.logspace(0, 10, num=11, base=0.5):
    simulation.context.setParameter('K', scaled_k)
    simulation.minimizeEnergy(maxIterations=100)
    print('K:', scaled_k, simulation.context.getState(getEnergy=True).getPotentialEnergy())

simulation.context.setParameter('K', 0)
simulation.minimizeEnergy(maxIterations=1000)

print('Final energy after minimization:', simulation.context.getState(getEnergy=True).getPotentialEnergy())

final_time = time.time()
elapsed_time = (final_time - initial_time) * u.seconds
print('Completed minimization in %8.3f s' % (elapsed_time / u.seconds))

print('Saving checkpoint.')
simulation.saveCheckpoint('minimized.chk')

