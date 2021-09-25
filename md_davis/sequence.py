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

import os
import warnings

import click
import mdtraj
import string
from md_davis.amino_acids import THREE_TO_ONE


def get_sequence(structure, label=None, return_type=None):
    """Get the sequence of amino acid residues from a PDB file

    Determines the sequence from the order of amino acid residues in the
    Protein Data Bank (PDB) file. The SEQRES records in the PDB file are
    not used, because it is sometimes different from the actual list of
    amino acids in the PDB file.

    Args
    ----
        structure : string
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
        label = os.path.splitext(os.path.basename(structure))[0]

    top = mdtraj.load(structure).topology
    output = {}
    sequence = {}
    for chain in top.chains:
        resi = []
        resn = []
        sequence_string = ''
        for residue in chain.residues:
            resi.append(residue.resSeq)
            resn.append(residue.name)
            if residue.name in THREE_TO_ONE:
                sequence_string += THREE_TO_ONE[residue.name]
            else:
                sequence_string += 'X'
                warnings.warn(f'Non standard residue {residue.name} found, it will be represented as X in the output')

        output[chain.index] = (resi, resn)
        sequence[chain.index] = sequence_string

    if return_type == 'dict':
        return dict(sequence)
    elif return_type == 'fasta':
        sequence_string = ''
        chain = string.ascii_uppercase + string.ascii_lowercase
        for ch, seq in sequence.items():
            sequence_string += f'>{label}|{chain[ch]}\n{seq}\n'
        return sequence_string
    elif return_type == 'modeller':
        return '/'.join(sequence.values())
    elif return_type == 'toml':
        toml_string = ''
        for ch, residues in output.items():
            toml_string += f"'Chain {ch}' = {str([[str(resi), str(resn)] for resi, resn in zip(*residues)])}\n"
        return toml_string
    else:
        return output


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(name='sequence', context_settings=CONTEXT_SETTINGS)
@click.option('-l', '--label', default=None, type=str)
@click.option('-r', '--return-type', 'return_type', default=None,
              type=click.Choice(['dict', 'fasta', 'modeller', 'toml'], case_sensitive=False),
              help="`dict` returns a dictionary containing chain labels "
                   "as keys and each chain's sequence as values. `fasta` "
                   "returns each chain's sequence separately in FASTA format.")
@click.argument('structure')
def main(structure, label, return_type):
    """Get the sequence of amino acid residues from a PDB file

    Determines the sequence from the order of amino acid residues in the
    Protein Data Bank (PDB) file. The SEQRES records in the PDB file are
    not used, because it is sometimes different from the actual list of
    amino acids in the PDB file.

    It is assumed that a PDB file used for dynamics will not contain
    any missing residues. The sequence gets truncated at the first
    missing residue
    """
    sequence = get_sequence(structure=structure, label=label,
                            return_type=return_type)
    if sequence is not None:
        click.echo(sequence)


if __name__ == '__main__':
    main()
