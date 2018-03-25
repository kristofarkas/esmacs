import os

import simtk.openmm.app as app
import simtk.unit as u

from esmacs import Minimizer, Heater, Equilibrate


def main():
    system_path = '/home/kristof/Research/INSPIRE/nilotinib/e255k/build'
    prmtop = app.AmberPrmtopFile(os.path.join(system_path, 'complex.top'))
    inpcrd = app.AmberInpcrdFile(os.path.join(system_path, 'complex.inpcrd'))
    system = prmtop.createSystem(nonbondedMethod=app.PME,
                                 constraints=app.HBonds,
                                 nonbondedCutoff=12*u.angstroms,
                                 switchDistance=10*u.angstroms)
    s1 = Minimizer(system=system, positions=inpcrd.positions, topology=prmtop.topology)

    context = s1.minimize()
    state = context.getState(getPositions=True)
    s2 = Heater(system=system, positions=state.getPositions(), topology=prmtop.topology)
    context = s2.heat(50, 300)

    state = context.getState(getPositions=True,	getVelocities=True)
    s3 = Equilibrate(system=system, positions=state.getPositions(asNumpy=True), topology=prmtop.topology)
    s3.equilibrate(1*u.atmosphere, 300*u.kelvin, state.getVelocities())


if __name__ == '__main__':
    main()
