#! /usr/bin/env python
"""
Get the sequence from a PDB file

Usage:
  md_davis sequence [-fd] <PDB_file>
  md_davis sequence -h | --help
  md_davis sequence --version

Options:
  -h --help     Show this screen.
  --version     Show version.
  -f --fasta    Print in fasta format
  -d --dict     Return a dictionary of chains and seqeuences
"""


from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import docopt
# import mdtraj


def get_sequence_dict(filename):
    """ Get the sequence from the PDB file and return a dictionary. """
    structure = PDBParser(QUIET=True).get_structure('structure', filename)
    chain_labels = [chain.id for chain in structure.get_chains()]
    ppb=PPBuilder()
    output = {}
    for chain, pp in zip(chain_labels, ppb.build_peptides(structure)):
        output[chain] = pp.get_sequence().__str__()
    return output


def get_fasta(filename):
    """ Get the sequence from the PDB file and return a fasta sequence. """
    structure = PDBParser(QUIET=True).get_structure('structure', filename)
    chain_labels = [chain.id for chain in structure.get_chains()]
    ppb=PPBuilder()
    output = []
    for chain, pp in zip(chain_labels, ppb.build_peptides(structure)):
        output.append(pp.get_sequence().__str__())
    return '/'.join(output)


# def get_sequence(filename):
#     """ Evaluate the dihedral angles using MDTraj and store it as HDF5 dataset """
#     top = mdtraj.load(filename, top=filename).topology
#     print(top)
#     return None


def main(arguments):
    if arguments['--dict']:
        print(get_sequence_dict(arguments['<PDB_file>']))
    else:
        print(get_fasta(arguments['<PDB_file>']))


if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    main(args)