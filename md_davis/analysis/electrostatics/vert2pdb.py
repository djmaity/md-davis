"""Convert MSMS surface vertices to PDB records for Delphi."""

from __future__ import annotations

import argparse
import sys
from os import PathLike
from typing import TextIO

import numpy
from biopandas.pdb import PandasPdb
from scipy.spatial import cKDTree


def vertices_to_pdb(
    vertex_filename: str | PathLike[str], pdb_filename: str | PathLike[str], output_stream: TextIO
) -> None:
    """Write each surface vertex using metadata from its closest PDB atom."""
    with open(vertex_filename, errors="replace") as vertex_file:
        vertices = numpy.loadtxt(vertex_file, usecols=(0, 1, 2), skiprows=3)

    atom_data = PandasPdb().read_pdb(pdb_filename).df["ATOM"]
    coordinates = atom_data[["x_coord", "y_coord", "z_coord"]]
    _, atom_indices = cKDTree(coordinates).query(vertices)

    for vertex, index in zip(vertices, atom_indices, strict=True):
        residue = atom_data.at[int(index), "residue_name"]
        chain = atom_data.at[int(index), "chain_id"]
        atom_name = atom_data.at[int(index), "atom_name"]
        residue_number = atom_data.at[int(index), "residue_number"]
        x, y, z = vertex
        print(
            f"ATOM        {atom_name:^4} {residue:3} {chain:1}{residue_number:4}    {x:8.3f}{y:8.3f}{z:8.3f}",
            file=output_stream,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert an MSMS surface to PDB records for Delphi")
    parser.add_argument("vert", metavar="filename.vert", help="MSMS vertices file")
    parser.add_argument("pdb", metavar="filename.pdb", help="Source PDB file")
    parser.add_argument(
        "-o",
        "--output",
        metavar="surface.pdb",
        type=argparse.FileType("w"),
        help="Output PDB file; defaults to standard output",
    )
    args = parser.parse_args()
    vertices_to_pdb(args.vert, args.pdb, args.output or sys.stdout)


if __name__ == "__main__":
    main()
