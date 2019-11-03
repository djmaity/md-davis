#! /usr/bin/env python
"""
Get the sequence from PDB files

Usage:
  md_davis sequence [-fd] <PDB_file>
  md_davis sequence -h | --help

Options:
  -h, --help        Show this screen.
  -f, --fasta       Print in fasta format
  -d, --dict        Return a dictionary of chains and seqeuences
"""


from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import docopt


def get_sequence(filename, return_dict=False):
    """ Get the sequence from the PDB file """
    try:
        structure = PDBParser(QUIET=True).get_structure('structure', filename)
    except FileNotFoundError:
        print(f'PDB File {filename} not found')
        return
    chain_labels = [chain.id for chain in structure.get_chains()]
    ppb=PPBuilder()
    output = []
    output_dict = {}
    for chain, pp in zip(chain_labels, ppb.build_peptides(structure)):
        output.append(pp.get_sequence().__str__())
        output_dict[chain] = pp.get_sequence().__str__()
    if return_dict:
        return output_dict
    else:
        return '/'.join(output)


def main(argv=None):
    if argv:
        args = docopt.docopt(__doc__, argv=argv)
    else:
        args = docopt.docopt(__doc__)

    if args['--dict']:
        seq = get_sequence(args['<PDB_file>'], return_dict=True)
        if seq != None:
            print(seq)
    else:
        seq = get_sequence(args['<PDB_file>'])
        if seq != None:
            print(seq)


if __name__ == '__main__':
    main()
