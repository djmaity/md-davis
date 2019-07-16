#! /usr/bin/env python

""" This script calculates and lists the RMSD between that first PDB
    structure provided with the rest of the PDB structures provided """

import sys
import argparse
# import os


def initialize_pymol():
    """ The commands necessary to get pymol running """

    # # autocompletion
    # import readline
    # import rlcompleter
    # readline.parse_and_bind('tab: complete')

    # # pymol environment
    # moddir='/opt/pymol-svn/modules'
    # sys.path.insert(0, moddir)
    # os.environ['PYMOL_PATH'] = os.path.join(moddir, 'pymol/pymol_path')

    # pymol launching: quiet (-q), without GUI (-c) and with arguments from command line
    import pymol
    pymol.pymol_argv = ['pymol','-cQ'] + sys.argv[1:]
    # pymol.pymol_argv = ['pymol'] + sys.argv[1:]
    pymol.finish_launching()
    return pymol.cmd


def calculate_rmsd(target, pdblist):
    """ Align the two pdb files using PyMOL and report the RMSD """
    cmd = initialize_pymol()
    cmd.load(target, object='target')
    rmsd = []
    for mobile in pdblist:
        cmd.delete('mobile')
        cmd.load(mobile, object='mobile')
        rmsd.append(cmd.align('mobile', 'target')[0])
    return rmsd


def main():
    parser = argparse.ArgumentParser(description='Find the RMSD between two PDB structures.')
    parser.add_argument('target', help='Target PDB structure')
    parser.add_argument('pdblist', nargs='+', help='List of PDB structures')
    args = parser.parse_args()

    rmsd_array = calculate_rmsd(args.target, args.pdblist)
    for pdb, rmsd in zip(args.pdblist, rmsd_array):
        print(f'{rmsd} Å is the RMSD of {pdb} with {args.target}')


if __name__ == "__main__":
    main()


