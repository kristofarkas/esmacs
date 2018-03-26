#!/usr/bin/env python

import os
import argparse

from esmacs import Esmacs


def run_esmacs():
    args = parse_args()

    os.makedirs(os.path.join(args.system, args.replica), exist_ok=True)
    os.chdir(os.path.join(args.system, args.replica))

    # root = '/lustre/atlas/scratch/farkaspall/chm126/inspire-data/nilotinib/{}/build'
    root = '/home/kristof/Research/INSPIRE/nilotinib/{}/build'
    prmtop = os.path.join(root, 'complex.top').format(args.system)
    inpcrd = os.path.join(root, 'complex.inpcrd').format(args.system)

    sim = Esmacs.from_amber(prmtop, inpcrd)
    sim.run_protocol(short_run=args.short)


def parse_args():
    parser = argparse.ArgumentParser(description='Run ESMACS simulation.')
    parser.add_argument('--system', required=True)
    parser.add_argument('--replica', required=True)
    parser.add_argument('--short', default=False, type=bool)

    return parser.parse_args()
