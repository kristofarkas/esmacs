import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u

import numpy as np
import mdtraj
import time

# Parameters

Ti = 50 * u.kelvin
T = 300 * u.kelvin
P = 1 * u.atmosphere

ts = 2 * u.femtosecond

kT = T * u.BOLTZMANN_CONSTANT_kB * u.AVOGADRO_CONSTANT_NA
sigma = 3 * u.angstrom

# System

prmtop = app.AmberPrmtopFile('complex.top')
pdb = app.PDBFile('complex.pdb')
system = prmtop.createSystem(nonbondedMethod=app.PME,
                             constraints=app.HBonds,
                             nonbondedCutoff=12*u.angstroms,
                             switchDistance=10*u.angstroms)
topology = mdtraj.Topology.from_openmm(prmtop.topology)

atoms_to_restrain = topology.select('not water and not type H')
K = kT / sigma**2

restraint_force = mm.CustomExternalForce('K*periodicdistance(x, y, x, x0, y0, z0)^2')
restraint_force.addGlobalParameter('K', K)
restraint_force.addPerParticleParameter('x0')
restraint_force.addPerParticleParameter('y0')
restraint_force.addPerParticleParameter('z0')

for atom_index in atoms_to_restrain:
    position = pdb.positions[atom_index]
    restraint_force.addParticle(int(atom_index), position.value_in_unit_system(u.md_unit_system))

system.addForce(restraint_force)

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

print('Saving PDB.')
positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(prmtop.topology, positions, open('minimized.pdb', 'w'), keepIds=True)

