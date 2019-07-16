#! /usr/bin/env python
""" Argument parser to get a PDB filename from the command line. """

import argparse


def get_pdb_filename(description):
    """ Parse the PDB filename from the command line. """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('file', help='input file', metavar='<PDB File>')
    return parser.parse_args()



