#! /usr/bin/env python

"""This script calculates and lists the RMSD between that first PDB
structure provided with the rest of the PDB structures provided"""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from md_davis.integrations.pymol import initialize_pymol


def align_structures(target: str, pdblist: Sequence[str]) -> None:
    """Align the two pdb files using PyMOL and report the RMSD"""
    cmd = initialize_pymol(window=False)
    cmd.set("retain_order", 1)
    cmd.load(target)
    target = cmd.get_object_list()[-1]
    for mobile in pdblist:
        cmd.load(mobile)
        last_obj = cmd.get_object_list()[-1]
        cmd.align(last_obj, target)
        cmd.save(filename=last_obj + ".pdb", selection=last_obj)


def main() -> None:
    parser = argparse.ArgumentParser(description="Find the RMSD between two PDB structures.")
    parser.add_argument("target", help="Target PDB structure")
    parser.add_argument("pdblist", nargs="+", help="List of PDB structures")
    args = parser.parse_args()

    align_structures(args.target, args.pdblist)


if __name__ == "__main__":
    main()
