import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u

import numpy as np
import mdtraj
import time

# System

prmtop = app.AmberPrmtopFile('complex.top')
system = prmtop.createSystem(nonbondedMethod=app.PME,
                             constraints=app.HBonds,
                             nonbondedCutoff=12*u.angstroms,
                             switchDistance=10*u.angstroms)

topology = mdtraj.Topology.from_openmm(prmtop.topology)

total_steps = 10_000 # Reducing for testing purposes from 3M


integrator = mm.LangevinIntegrator(300*u.kelvin, 5/u.picosecond, 0.002*u.picoseconds)
barostat = mm.MonteCarloBarostat(1.0*u.bar, 300*u.kelvin)


atoms_to_restrain = topology.select('not water and not type H')
default_k = 4.0*u.kilocalories_per_mole/u.angstroms**2

harmonic_restraint = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
harmonic_restraint.addGlobalParameter('k', default_k)
harmonic_restraint.addPerParticleParameter("x0")
harmonic_restraint.addPerParticleParameter("y0")
harmonic_restraint.addPerParticleParameter("z0")

system.addForce(harmonic_restraint)
system.addForce(barostat)

# Initialise simulation and positions

simulation = app.Simulation(prmtop.topology, system, integrator)

simulation.reporters.append(DCDReporter('equilibrate.dcd', 10))

simulation.loadState('heated.xml')

state = simulation.context.getState(getPositions=True)
positions = state.getPositions()

for restraint in range(harmonic_restraint.getNumParticles()):
    atomindex, _ = harmonic_restraint.getParticleParameters(restraint)
    position = positions[atomindex]
    harmonic_restraint.setParticleParameters(restraint ,atomindex, position.value_in_unit_system(u.md_unit_system))

harmonic_restraint.updateParametersInContext(simulation.context)

# Equilibrate

print('Equilibrate...')
initial_time = time.time()

for scaled_k in default_k*10*np.logspace(0, 10, num=11, base=0.5):
    simulation.context.setParameter('k', scaled_k)
    simulation.step(int(total_steps/60))
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy(), temperature/u.kelvin)


simulation.context.setParameter('k', 0)
simulation.step(int(total_steps*0.15666666))

final_time = time.time()
elapsed_time = (final_time - initial_time) * u.seconds
print('Completed equilibration in %8.3f s' % (elapsed_time / u.seconds))
print('Potential energy is %.3f kcal/mol' % (simulation.context.getState(getEnergy=True).getPotentialEnergy() / u.kilocalories_per_mole))

print('Saving state.')
simulation.saveState('equilibrated.xml')
