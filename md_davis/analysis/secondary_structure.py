#!/usr/bin/env python

"""
This module parses the strings in .dat output of gmx do_dssp and
plots the secondary structure exhibited by each residue.

File name: secStr_per_residue.py
Author: Dibyajyoti Maity
Date created: 21 February 2019
Python Version: 3.7
"""

from __future__ import annotations

import argparse
import pickle
import statistics
from collections import Counter
from os import PathLike

SECSTR_CODES = {
    "H": "α-helix",
    "G": "3_10-helix",
    "I": "π-helix",
    "B": "β-bridge",
    "E": "β strand",
    "T": "Turn",
    "S": "Bend",
    "~": "Loop",
}


def parse_dat(filename: str | PathLike[str]) -> list[list[str]]:
    """Parse the DSSP data file"""
    secondary_structure = []
    with open(filename) as dat_file:
        next(dat_file)
        for line in dat_file:
            secondary_structure.append(list(line.strip()))
    return secondary_structure


def get_counts_per_residue(ss_strings: list[list[str]]) -> list[Counter[str]]:
    output = []
    for ss in map(list, zip(*ss_strings, strict=False)):
        output.append(Counter(ss))
    return output


def get_ss_percentage(ss_strings: list[list[str]], start: int = 0, end: int | None = None) -> str:
    output = {
        "Helix": [],
        "Sheet": [],
        "Turn + Bend": [],
        "Loop": [],
    }
    for ss in map(list, ss_strings):
        counter = Counter(ss[start:end])
        total = sum(counter.values())
        output["Helix"].append((counter["H"] + counter["G"] + counter["I"]) * 100 / total)
        output["Sheet"].append((counter["E"] + counter["B"]) * 100 / total)
        output["Turn + Bend"].append((counter["T"] + counter["S"]) * 100 / total)
        output["Loop"].append((counter["~"]) * 100 / total)
    output_string = ""
    for structure in output.keys():
        mean = statistics.mean(output[structure])
        std = statistics.stdev(output[structure])
        output_string += f"{structure:<12}: {mean:.2f} % +/- {std:.2f} %\n"
    return output_string


def get_arguments() -> argparse.Namespace:
    """Get input and output filename from the commandline"""
    parser = argparse.ArgumentParser(
        description="Calculate and plot the secondary structure exhibited by each residue in the output of gmx do_dssp"
    )
    parser.add_argument("filename", help=".dat file output of gmx do_dssp")
    parser.add_argument("-o", "--output", help="name and extension for output pickle file")
    return parser.parse_args()


def main() -> None:
    """Parse and plot the DSSP data"""
    args = get_arguments()
    data = parse_dat(args.filename)

    ss_percentage = get_ss_percentage(data)
    print(ss_percentage)

    if args.output:
        ss_counts = get_counts_per_residue(data)
        pickle.dump(ss_counts, open(args.output, "wb"))


if __name__ == "__main__":
    main()
