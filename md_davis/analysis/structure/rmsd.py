#! /usr/bin/env python

"""This script calculates and lists the RMSD between that first PDB
structure provided with the rest of the PDB structures provided"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from md_davis.integrations.pymol import initialize_pymol


def calculate_rmsd(target: str, pdblist: Sequence[str]) -> list[float]:
    """Align the two pdb files using PyMOL and report the RMSD"""
    cmd = initialize_pymol(extra_args=sys.argv[1:])
    cmd.load(target, object="target")
    rmsd = []
    for mobile in pdblist:
        cmd.delete("mobile")
        cmd.load(mobile, object="mobile")
        rmsd.append(cmd.align("mobile", "target")[0])
    return rmsd


def main() -> None:
    parser = argparse.ArgumentParser(description="Find the RMSD between two PDB structures.")
    parser.add_argument("target", help="Target PDB structure")
    parser.add_argument("pdblist", nargs="+", help="List of PDB structures")
    args = parser.parse_args()

    rmsd_array = calculate_rmsd(args.target, args.pdblist)
    for pdb, rmsd in zip(args.pdblist, rmsd_array, strict=False):
        print(f"{rmsd} Å is the RMSD of {pdb} with {args.target}")


if __name__ == "__main__":
    main()
