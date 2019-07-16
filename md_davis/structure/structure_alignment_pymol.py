#! /usr/bin/env python

""" This script calculates and lists the RMSD between that first PDB
    structure provided with the rest of the PDB structures provided """

import sys
import argparse
import pymol


def initialize_pymol(window=False):
    """ The commands necessary to get pymol running """
    if window:
        pymol.pymol_argv = ['pymol']
    else:
        pymol.pymol_argv = ['pymol','-cQ']
    pymol.finish_launching()
    return pymol.cmd


def align_structures(target, pdblist):
    """ Align the two pdb files using PyMOL and report the RMSD """
    cmd = initialize_pymol(window=False)
    cmd.set('retain_order', 1)
    cmd.load(target)
    target = cmd.get_object_list()[-1]
    for mobile in pdblist:
        cmd.load(mobile)
        last_obj = cmd.get_object_list()[-1]
        cmd.align(last_obj, target)
        cmd.save(filename=last_obj + '.pdb', selection=last_obj)


def main():
    parser = argparse.ArgumentParser(description='Find the RMSD between two PDB structures.')
    parser.add_argument('target', help='Target PDB structure')
    parser.add_argument('pdblist', nargs='+', help='List of PDB structures')
    args = parser.parse_args()

    align_structures(args.target, args.pdblist)


if __name__ == "__main__":
    main()


