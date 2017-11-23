import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u

import numpy as np
import mdtraj
import time

# System

prmtop = app.AmberPrmtopFile('complex.top')
pdb = app.PDBFile('minimized.pdb')
system = prmtop.createSystem(nonbondedMethod=app.PME,
                             constraints=app.HBonds,
                             nonbondedCutoff=12*u.angstroms,
                             switchDistance=10*u.angstroms)

topology = mdtraj.Topology.from_openmm(prmtop.topology)

integrator = mm.LangevinIntegrator(50*u.kelvin, 5/u.picosecond, 0.002*u.picoseconds)

atoms_to_restrain = topology.select('not water and not type H')
default_k = 4.0*u.kilocalories_per_mole/u.angstroms**2

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

simulation.reporters.append(app.DCDReporter('heating.dcd', 50))


# Heating

print('Heating...')

initial_time = time.time()

for temperature in np.linspace(50, 300, 251)*u.kelvin:
    integrator.setTemperature(temperature)
    simulation.step(100)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy(), temperature/u.kelvin)

simulation.step(500)

final_time = time.time()
elapsed_time = (final_time - initial_time) * u.seconds
print('Completed heating in %8.3f s' % (elapsed_time / u.seconds))
print('Potential energy is %.3f kcal/mol' % (simulation.context.getState(getEnergy=True).getPotentialEnergy() / u.kilocalories_per_mole))

print('Saving state.')
simulation.saveState('heated.xml')
