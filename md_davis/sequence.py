#!/usr/bin/env python3
"""
Usage: sequence.py [OPTIONS] PDB_FILE

  Get the sequence of amino acid residues from a PDB file

  Determines the sequence from the order of amino acid residues in the
  Protein Data Bank (PDB) file. The SEQRES records in the PDB file are not
  used, because it is sometimes different from the actual list of amino
  acids in the PDB file.

  It is assumed that a PDB file used for dynamics will not contain any
  missing residues. The sequence gets truncated at the first missing residue

Options:
  -l, --label TEXT
  -r, --return-type [dict|fasta]  `dict` returns a dictionary containing chain
                                  labels as keys and each chain's sequence as
                                  values. `fasta` returns each chain's sequence
                                  separately in FASTA format.

  -h, --help                      Show this message and exit.
"""

# import sys
# import warnings

# if not sys.warnoptions:
#     warnings.filterwarnings('ignore',
#                             category=RuntimeWarning,
#                             module='runpy')

import os
import click
from collections import OrderedDict
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder


def get_sequence(pdb_file, label=None, return_type=None):
    """Get the sequence of amino acid residues from a PDB file

    Determines the sequence from the order of amino acid residues in the
    Protein Data Bank (PDB) file. The SEQRES records in the PDB file are
    not used, because it is sometimes different from the actual list of
    amino acids in the PDB file.

    Args
    ----
        pdb_file : string
            Input PDB file.
        label : string, optional
            Label for the sequence when returning FASTA format.

    Parameters
    ----------
        return_type : {'dict', 'fasta'}, optional
            `dict` returns a dictionary containing chain labels as keys
            and each chain's sequence as values. `fasta` returns each
            chain's sequence separately in FASTA format.

    Returns
    -------
        sequence: string or dict
            Sequence of amino acids in the PDB file

    Warns
    -----
        It is assumed that a PDB file used for dynamics will not contain
        any missing residues. The sequence gets truncated at the first
        missing residue
    """
    if not label:
        label = os.path.splitext(os.path.basename(pdb_file))[0]

    structure = PDBParser(QUIET=True).get_structure('structure', pdb_file)
    chain_labels = [chain.id for chain in structure.get_chains()]
    ppb = PPBuilder()
    sequence = OrderedDict()
    for chain, pp in zip(chain_labels, ppb.build_peptides(structure)):
        sequence[chain] = pp.get_sequence().__str__()

    if return_type == 'dict':
        return dict(sequence)
    elif return_type == 'fasta':
        sequence_string = ''
        for chain, seq in sequence.items():
            sequence_string += f'>{label}'
            if chain != ' ':
                sequence_string += f'|{chain}'
            sequence_string += f'\n{seq}\n'
        return sequence_string
    else:
        return '/'.join(sequence.values())


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(name='sequence', context_settings=CONTEXT_SETTINGS)
@click.option('-l', '--label', default=None, type=str)
@click.option('-r', '--return-type', 'return_type', default=None,
              type=click.Choice(['dict', 'fasta'], case_sensitive=False),
              help="`dict` returns a dictionary containing chain labels "
                   "as keys and each chain's sequence as values. `fasta` "
                   "returns each chain's sequence separately in FASTA format.")
@click.argument('pdb_file')
def main(pdb_file, label, return_type):
    """Get the sequence of amino acid residues from a PDB file

    Determines the sequence from the order of amino acid residues in the
    Protein Data Bank (PDB) file. The SEQRES records in the PDB file are
    not used, because it is sometimes different from the actual list of
    amino acids in the PDB file.

    It is assumed that a PDB file used for dynamics will not contain
    any missing residues. The sequence gets truncated at the first
    missing residue
    """
    sequence = get_sequence(pdb_file=pdb_file, label=label,
                            return_type=return_type)
    if sequence is not None:
        click.echo(sequence)


if __name__ == '__main__':
    main()
