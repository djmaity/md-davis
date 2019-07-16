#! /usr/bin/env python

# import argparse
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from my_parser import get_pdb_filename


def get_sequence(filename):
    """ Get the sequence from the PDB file. """
    structure = PDBParser(QUIET=True).get_structure('structure', filename)
    chain_labels = [chain.id for chain in structure.get_chains()]
    ppb=PPBuilder()
    output = {}
    for chain, pp in zip(chain_labels, ppb.build_peptides(structure)):
        output[chain] = pp.get_sequence().__str__()
    return output


def main():
    args = get_pdb_filename('Get the sequence of residues from the PDB file.')
    sequence = get_sequence(args.file)
    print(sequence)


if __name__ == '__main__':
    main()