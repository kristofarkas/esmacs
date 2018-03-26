import os
import argparse

from esmacs import Esmacs


def main():
    args = parse_args()

    os.makedirs(os.path.join(args.system, args.replica), exist_ok=True)
    os.chdir(os.path.join(args.system, args.replica))

    # root = '/lustre/atlas/scratch/farkaspall/chm126/inspire-data/nilotinib/{}/build'
    root = '/home/kristof/Research/INSPIRE/nilotinib/{}/build'
    prmtop = os.path.join(root, 'complex.top').format(args.system)
    inpcrd = os.path.join(root, 'complex.inpcrd').format(args.system)

    sim = Esmacs.from_amber(prmtop, inpcrd)
    sim.run_protocol()
    

def parse_args():
    parser = argparse.ArgumentParser(description='Run ESMACS simulation.')
    parser.add_argument('--systems', required=True)
    parser.add_argument('--replica', required=True)

    return parser.parse_args()


if __name__ == '__main__':
    main()
