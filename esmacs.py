import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u

import parmed as pmd
import numpy as np
import mdtraj
import time

# System

prmtop = app.AmberPrmtopFile('complex.top')
parmed_parm = pmd.load_file('complex.top')
pdb = app.PDBFile('complex.pdb')
system = prmtop.createSystem(nonbondedMethod=app.PME,
                             constraints=app.HBonds,
                             nonbondedCutoff=12*u.angstroms,
                             switchDistance=10*u.angstroms)
topology = mdtraj.Topology.from_openmm(prmtop.topology)

total_steps = 1000 # Reducing for testing purposes from 3M

# Integrator, forces

integrator = mm.LangevinIntegrator(50*u.kelvin, 5/u.picosecond, 0.002*u.picoseconds)

barostat = mm.MonteCarloBarostat(1.0*u.bar, 300*u.kelvin)

atoms_to_restrain = topology.select('not water and not type H')
default_k = 0.4*u.kilocalories_per_mole/u.angstroms**2

harmonic_restraint = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
harmonic_restraint.addGlobalParameter('k', default_k)
harmonic_restraint.addPerParticleParameter("x0")
harmonic_restraint.addPerParticleParameter("y0")
harmonic_restraint.addPerParticleParameter("z0")

for atomindex in atoms_to_restrain:
    position = pdb.positions[atomindex]
    harmonic_restraint.addParticle(int(atomindex), position.value_in_unit_system(u.md_unit_system))

system.addForce(harmonic_restraint)

# Initialise simulation and positions

simulation = app.Simulation(prmtop.topology, system, integrator)

simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(50*u.kelvin)

# Minimise

print('Minimizing, from inital energy:', simulation.context.getState(getEnergy=True).getPotentialEnergy())

initial_time = time.time()

for scaled_k in default_k*10*np.logspace(0, 10, num=11, base=0.5):
    simulation.context.setParameter('k', scaled_k)
    simulation.minimizeEnergy(maxIterations=10)
    print('k:', scaled_k, simulation.context.getState(getEnergy=True).getPotentialEnergy())
    
simulation.context.setParameter('k', 0)
simulation.minimizeEnergy()

print('Final energy after minimization:', simulation.context.getState(getEnergy=True).getPotentialEnergy())

final_time = time.time()
elapsed_time = (final_time - initial_time) * u.seconds
print('Completed minimization in %8.3f s' % (elapsed_time / u.seconds))

# Heating

print('Heating...')

state = simulation.context.getState(getPositions=True)
positions = state.getPositions()

for restraint in range(harmonic_restraint.getNumParticles()):
    atomindex, _ = harmonic_restraint.getParticleParameters(restraint)
    position = positions[atomindex]
    harmonic_restraint.setParticleParameters(restraint ,atomindex, position.value_in_unit_system(u.md_unit_system))

harmonic_restraint.updateParametersInContext(simulation.context)
simulation.context.setParameter('k', default_k)

initial_time = time.time()

for temperature in np.linspace(50, 300, 251)*u.kelvin:
    # integrator.setTemperature(temperature)
    simulation.step(2)
    
    print(pmd.openmm.energy_decomposition(parmed_parm, simulation.context))
    #print(simulation.context.getState(getEnergy=True).getKineticEnergy(),
    #      simulation.context.getState(getEnergy=True).getPotentialEnergy(),
    #      temperature/u.kelvin)

simulation.step(500)

final_time = time.time()
elapsed_time = (final_time - initial_time) * u.seconds
print('Completed heating in %8.3f s' % (elapsed_time / u.seconds))
print('Potential energy is %.3f kcal/mol' % (simulation.context.getState(getEnergy=True).getPotentialEnergy() / u.kilocalories_per_mole))

# Equilibrate

print('Equilibrate...')
initial_time = time.time()

simulation.system.addForce(barostat)

state = simulation.context.getState(getPositions=True)
positions = state.getPositions()

for restraint in range(harmonic_restraint.getNumParticles()):
    atomindex, _ = harmonic_restraint.getParticleParameters(restraint)
    position = positions[atomindex]
    harmonic_restraint.setParticleParameters(restraint ,atomindex, position.value_in_unit_system(u.md_unit_system))

harmonic_restraint.updateParametersInContext(simulation.context)

for scaled_k in default_k*10*np.logspace(0, 10, num=11, base=0.5):
    simulation.context.setParameter('k', scaled_k)
    simulation.step(int(total_steps/60))
    
simulation.context.setParameter('k', 0)
simulation.step(int(total_steps*0.15666666))

final_time = time.time()
elapsed_time = (final_time - initial_time) * u.seconds
print('Completed equilibration in %8.3f s' % (elapsed_time / u.seconds))
print('Potential energy is %.3f kcal/mol' % (simulation.context.getState(getEnergy=True).getPotentialEnergy() / u.kilocalories_per_mole))

# Production
print('Production...')
initial_time = time.time()

simulation.step(int(total_steps*2/3))

final_time = time.time()
elapsed_time = (final_time - initial_time) * u.seconds
simulated_time = total_steps * integrator.getStepSize()
performance = (simulated_time / elapsed_time)
print('Completed %8d steps in %8.3f s : performance is %8.3f ns/day' % (total_steps, elapsed_time / u.seconds, performance / (u.nanoseconds/u.day)))
print('Final potential energy is %.3f kcal/mol' % (simulation.context.getState(getEnergy=True).getPotentialEnergy() / u.kilocalories_per_mole))
