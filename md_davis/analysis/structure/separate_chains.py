"""
Created on Mon Mar 19 07:07:17 IST 2018

This script separates the chains in PDB file and writes in to separate files
with chain id as suffix.

Usage: python separate_chains.py <PDB FILE>

@author: Dibyajyoti Maity
"""

import Bio.PDB

from md_davis.analysis.structure.parser import get_pdb_filename


def separate_chains(filename: str) -> None:
    """Separate chains in a PDB files and write it to individual files."""
    pdb_parser = Bio.PDB.PDBParser()
    io = Bio.PDB.PDBIO()
    structure = pdb_parser.get_structure("structure", filename)
    for model in structure:
        for chain in model:
            io.set_structure(chain)
            io.save(filename + f"_chain{chain.id}.pdb")


def main() -> None:
    args = get_pdb_filename("Separate chains in a PDB files and write it to individual files.")
    separate_chains(args.file)


if __name__ == "__main__":
    main()
