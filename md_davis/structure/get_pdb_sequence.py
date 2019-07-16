#! /usr/bin/env python

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from .my_parser import get_pdb_filename
import mdtraj


def get_sequence_dict(filename):
    """ Get the sequence from the PDB file. """
    structure = PDBParser(QUIET=True).get_structure('structure', filename)
    chain_labels = [chain.id for chain in structure.get_chains()]
    ppb=PPBuilder()
    output = {}
    for chain, pp in zip(chain_labels, ppb.build_peptides(structure)):
        output[chain] = pp.get_sequence().__str__()
    return output


def get_fasta(filename):
    """ Get the sequence from the PDB file. """
    structure = PDBParser(QUIET=True).get_structure('structure', filename)
    chain_labels = [chain.id for chain in structure.get_chains()]
    ppb=PPBuilder()
    output = []
    for chain, pp in zip(chain_labels, ppb.build_peptides(structure)):
        output.append(pp.get_sequence().__str__())
    return '/'.join(output)


def get_sequence(group, trajectory, structure, atoms=None):
    """ Evaluate the dihedral angles using MDTraj and store it as HDF5 dataset """
    top = mdtraj.load(structure, top=structure, atom_indices=atoms).topology
    print(top)
    return None


def main():
    args = get_pdb_filename('Get the sequence of residues from the PDB file.')
    # print(get_sequence_dict(args.file))
    # print(get_fasta(args.file))

    print(get_sequence(args.file))



if __name__ == '__main__':
    main()