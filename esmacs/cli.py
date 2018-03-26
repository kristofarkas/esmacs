#!/usr/bin/env python

import os
import argparse

from esmacs import Esmacs


def run_esmacs():
    args = parse_args()

    os.makedirs(os.path.join(args.system, args.replica), exist_ok=True)
    os.chdir(os.path.join(args.system, args.replica))

    root = os.path.join(args.root, args.drug, args.mutation)
    # root = '/lustre/atlas/scratch/farkaspall/chm126/inspire-data/nilotinib/{}/build'
    # root = '/home/kristof/Research/INSPIRE/nilotinib/{}/build'
    prmtop = os.path.join(root, 'build/complex.top')
    inpcrd = os.path.join(root, 'build/complex.inpcrd')

    sim = Esmacs.from_amber(prmtop, inpcrd)
    sim.run_protocol(short_run=args.short)


def parse_args():
    parser = argparse.ArgumentParser(description='Run ESMACS simulation.')
    parser.add_argument('--root', required=True)
    parser.add_argument('--drug', required=True)
    parser.add_argument('--mutation', required=True)
    parser.add_argument('--replica', required=True)
    parser.add_argument('--short', action='store_true')

    return parser.parse_args()
