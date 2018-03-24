import numpy as np
from simtk import openmm, unit


# =============================================================================
# RESTRAIN ATOMS
# Original at https://github.com/choderalab/openmmtools/blob/master/openmmtools/forcefactories.py
# =============================================================================

def restrain_atoms_by_dsl(system, pressure, positions, topology, atoms_dsl, constant):
    # Determine indices of the atoms to restrain.
    restrained_atoms = topology.select(atoms_dsl).tolist()
    return restrain_atoms(system, pressure, positions, restrained_atoms, constant)


def restrain_atoms(system, pressure, positions, restrained_atoms, constant):
    """Apply a soft harmonic restraint to the given atoms.
    This modifies the ``ThermodynamicState`` object.
    Parameters
    ----------
    system : openmm.app.System
        The thermodynamic state with the system. This will be modified.
    positions : np.array or list
        List of the positions.
    pressure : bool
        The system is NPT or not.
    restrained_atoms: list
        List of indexes of atoms to restrain
    constant: harmonic restraint constant

    Returns
    -------
    positions: list
        In NPT ensemble the atoms have to be moved for the barostat to work.
    """
    K = constant  # thermodynamic_state.kT / sigma**2  # Spring constant.
    # system = thermodynamic_state.system  # This is a copy.

    # Check that there are atoms to restrain.
    if len(restrained_atoms) == 0:
        raise ValueError('No atoms to restrain.')

    # We need to translate the restrained molecule to the origin
    # to avoid MonteCarloBarostat rejections (see openmm#1854).
    if pressure:
        # First, determine all the molecule atoms. Reference platform is the cheapest to allocate?
        reference_platform = openmm.Platform.getPlatformByName('Reference')
        integrator = openmm.VerletIntegrator(1.0*unit.femtosecond)
        context = openmm.Context(system, integrator, reference_platform)
        molecules_atoms = context.getMolecules()
        del context, integrator

        # Make sure the atoms to restrain belong only to a single molecule.
        molecules_atoms = [set(molecule_atoms) for molecule_atoms in molecules_atoms]
        restrained_atoms_set = set(restrained_atoms)
        restrained_molecule_atoms = None
        for molecule_atoms in molecules_atoms:
            if restrained_atoms_set.issubset(molecule_atoms):
                # Convert set to list to use it as numpy array indices.
                restrained_molecule_atoms = list(molecule_atoms)
                break
        if restrained_molecule_atoms is None:
            raise ValueError('Cannot match the restrained atoms to any molecule. Restraining '
                             'two molecules is not supported when using a MonteCarloBarostat.')

        # Translate system so that the center of geometry is in
        # the origin to reduce the barostat rejections.
        distance_unit = positions.unit
        centroid = np.mean(positions[restrained_molecule_atoms, :] / distance_unit, axis=0)
        positions -= centroid * distance_unit

    # Create a CustomExternalForce to restrain all atoms.
    energy_expression = '(K/2)*periodicdistance(x, y, z, x0, y0, z0)^2'  # periodic distance
    restraint_force = openmm.CustomExternalForce(energy_expression)
    # Adding the spring constant as a global parameter allows us to turn it off if desired
    restraint_force.addGlobalParameter('K', K)
    restraint_force.addPerParticleParameter('x0')
    restraint_force.addPerParticleParameter('y0')
    restraint_force.addPerParticleParameter('z0')
    for index in restrained_atoms:
        parameters = positions[index].value_in_unit_system(unit.md_unit_system)
        restraint_force.addParticle(index, parameters)

    # Update thermodynamic state.
    system.addForce(restraint_force)

    return positions
