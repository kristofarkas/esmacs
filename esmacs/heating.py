import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u

import numpy as np
import mdtraj
import time

from esmacs import restrain

# System

system_path = '/home/kristof/Research/INSPIRE/nilotinib/e255k/build'
prmtop = app.AmberPrmtopFile(os.path.join(system_path, 'complex.top'))
system = prmtop.createSystem(nonbondedMethod=app.PME,
                             constraints=app.HBonds,
                             nonbondedCutoff=12*u.angstroms,
                             switchDistance=10*u.angstroms)

topology = mdtraj.Topology.from_openmm(prmtop.topology)

integrator = mm.LangevinIntegrator(50*u.kelvin, 5/u.picosecond, 0.002*u.picoseconds)

K = 4 * u.kilocalorie/(u.angstrom**2*u.mole)  # kT / sigma**2

restrain.restrain_atoms_by_dsl(system,
                               pressure=False,
                               positions=pdb.positions,
                               atoms_dsl='not water and not type H',
                               topology=topology,
                               constant=K)


# Initialise simulation and positions

simulation = app.Simulation(prmtop.topology, system, integrator)

simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(50*u.kelvin)

simulation.reporters.append(app.DCDReporter('heating.dcd', 10))


# Heating

print('Heating, from inital energy:', 
     simulation.context.getState(getEnergy=True).getPotentialEnergy())

initial_time = time.time()

for temperature in np.linspace(50, 300, 251)*u.kelvin:
    integrator.setTemperature(temperature)
    simulation.step(1)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy(), temperature/u.kelvin)

simulation.step(500)

final_time = time.time()
elapsed_time = (final_time - initial_time) * u.seconds
print('Completed heating in %8.3f s' % (elapsed_time / u.seconds))
print('Potential energy is %.3f kcal/mol' % (simulation.context.getState(getEnergy=True).getPotentialEnergy() / u.kilocalories_per_mole))

print('Saving state.')
simulation.saveState('heated.xml')
