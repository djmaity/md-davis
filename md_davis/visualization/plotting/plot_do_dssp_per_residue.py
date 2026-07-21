#!/usr/bin/env python

"""
This module parses the strings in .dat output of gmx do_dssp and
plots the secondary structure exhibited by each residue.

File name: plot_do_dssp_per_residue.py
Author: Dibyajyoti Maity
Date created: 20 June 2017
Date last modified: 18 July 2018
Python Version: 3.6
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from os import PathLike

import numpy as np
from numpy.typing import NDArray

from md_davis.visualization.matplotlib import plot_secondary_structure_counts


def get_arguments() -> argparse.Namespace:
    """Get input and output filename from the commandline"""
    parser = argparse.ArgumentParser(
        description="Calculate and plot the secondary structure exhibited by each residue in the output of gmx do_dssp"
    )
    parser.add_argument("filename", help=".dat file output of gmx do_dssp")
    parser.add_argument("-o", "--output", help="name and extension for output image")
    return parser.parse_args()


def secondary_structure_percent(secstr: Sequence[str], indices: Sequence[int]) -> NDArray:
    """Count the secondary structure exhibited by each residue throughout the trajectory."""
    code = [
        "H",  # α-helix
        "G",  # 3-helix (310 helix)
        "I",  # 5 helix (π-helix)
        "B",  # residue in isolated β-bridge
        "E",  # extended strand, participates in β ladder
        "T",  # hydrogen bonded turn
        "S",  # bend
        "~",
    ]  # loop or irregular
    len_code = 8  # len(code)
    frames = len(secstr)
    output = np.zeros((len_code, len(indices)))
    for i, ind in enumerate(indices):
        counts = {"H": 0, "B": 0, "E": 0, "G": 0, "I": 0, "T": 0, "S": 0, "~": 0}
        for f in range(frames):
            counts[secstr[f][ind]] = counts[secstr[f][ind]] + 1
        for j in range(len_code):
            output[j][i] = counts[code[j]]
    return np.divide(output, frames / 100)


def do_dssp_count(
    filename: str | PathLike[str], start: int = 0, end: int | None = None
) -> list[NDArray]:
    """Parse the DSSP data"""
    with open(filename) as dat_file:
        lines = dat_file.readlines()
        lines.pop(0)
        num_lines = len(lines)
    num_chains = len(lines[0].split("="))
    secstr = [[""] * num_lines for _ in range(num_chains)]
    for line_index, line in enumerate(lines):
        split_lines = line.strip().split("=")
        for chain in range(num_chains):
            secstr[chain][line_index] = split_lines[chain][start:end]
    output = []
    for chain_secstr in secstr:
        output.append(secondary_structure_percent(chain_secstr, range(len(chain_secstr[0]))))
    return output


def main() -> None:
    """Parse and plot the DSSP data"""
    args = get_arguments()
    counts = do_dssp_count(args.filename)

    plot_secondary_structure_counts(counts, output=args.output)


if __name__ == "__main__":
    main()
