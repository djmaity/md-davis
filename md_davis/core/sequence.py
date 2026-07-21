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

from __future__ import annotations

import os
import string
import warnings
from os import PathLike
from typing import Literal

import mdtraj

from md_davis.core.amino_acids import THREE_TO_ONE


ResidueSequence = dict[int, tuple[list[int], list[str]]]


def get_sequence(
    structure: str | PathLike[str],
    label: str | None = None,
    return_type: Literal["dict", "fasta", "modeller", "toml"] | None = None,
) -> ResidueSequence | dict[int, str] | str:
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
        sequence_string = ""
        for residue in chain.residues:
            resi.append(residue.resSeq)
            resn.append(residue.name)
            if residue.name in THREE_TO_ONE:
                sequence_string += THREE_TO_ONE[residue.name]
            else:
                sequence_string += "X"
                warnings.warn(
                    f"Non standard residue {residue.name} found, it will be represented as X in the output",
                    stacklevel=2,
                )

        output[chain.index] = (resi, resn)
        sequence[chain.index] = sequence_string

    if return_type == "dict":
        return dict(sequence)
    elif return_type == "fasta":
        sequence_string = ""
        chain = string.ascii_uppercase + string.ascii_lowercase
        for ch, seq in sequence.items():
            sequence_string += f">{label}|{chain[ch]}\n{seq}\n"
        return sequence_string
    elif return_type == "modeller":
        return "/".join(sequence.values())
    elif return_type == "toml":
        toml_string = ""
        for ch, residues in output.items():
            toml_string += (
                f"'Chain {ch}' = {str([[str(resi), str(resn)] for resi, resn in zip(*residues, strict=False)])}\n"
            )
        return toml_string
    else:
        return output


def main(structure: str, label: str | None = None, return_type: str | None = None) -> None:
    """Get the sequence of amino acid residues from a PDB file

    Determines the sequence from the order of amino acid residues in the
    Protein Data Bank (PDB) file. The SEQRES records in the PDB file are
    not used, because it is sometimes different from the actual list of
    amino acids in the PDB file.

    It is assumed that a PDB file used for dynamics will not contain
    any missing residues. The sequence gets truncated at the first
    missing residue
    """
    sequence = get_sequence(structure=structure, label=label, return_type=return_type)
    if sequence is not None:
        print(sequence)


if __name__ == "__main__":
    import sys

    from md_davis.cli import main as cli_main

    raise SystemExit(cli_main(["sequence", *sys.argv[1:]]))
