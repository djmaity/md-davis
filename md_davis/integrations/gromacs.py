from __future__ import annotations

import pathlib
import re
import subprocess
from os import PathLike

def first_frame(trajectory: str | PathLike[str], gmx_path: str | PathLike[str] = "gmx") -> float:
    """Get the time of the first frame of a GROMACS trajectory file"""
    process = subprocess.run([gmx_path, "dump", "-f", trajectory], capture_output=True, encoding="utf-8")
    m = re.search(r"time=([\d.e+-]+)", process.stdout)
    time = float(m.group(1))
    return time


def process_gromacs_trajectory(
    trajectory: str,
    tpr: str,
    output: str = "./",
    index: str | None = None,
    prefix: str = "",
    selection: str = "non-Water_&_!Ion",
    gmx_path: str = "gmx",
) -> None:
    """Process a GROMACS trajecotry to correct for periodic boundary condition."""

    pathlib.Path(output).mkdir(parents=True, exist_ok=True)
    if index is None:
        index = f"{output}/{prefix}index.ndx"
        if selection == "non-Water_&_!Ion":
            subprocess.run(
                [gmx_path, "make_ndx", "-f", tpr, "-o", index], input='"non-Water"&!"Ion"\nq\n', encoding="ascii"
            )
        if selection == "Protein":
            subprocess.run([gmx_path, "make_ndx", "-f", tpr, "-o", index], input='"Protein"\nq\n', encoding="ascii")

    begin = str(first_frame(trajectory, gmx_path=gmx_path))

    subprocess.run(
        [
            gmx_path,
            "trjconv",
            "-f",
            trajectory,
            "-s",
            tpr,
            "-o",
            f"{output}/{prefix}whole.pdb",
            "-n",
            index,
            "-pbc",
            "whole",
            "-dump",
            begin,
        ],
        input=selection,
        encoding="ascii",
    )

    subprocess.run(
        [
            gmx_path,
            "trjconv",
            "-f",
            trajectory,
            "-s",
            tpr,
            "-o",
            f"{output}/{prefix}whole.xtc",
            "-n",
            index,
            "-pbc",
            "whole",
        ],
        input=selection,
        encoding="ascii",
    )

    subprocess.run(
        [
            gmx_path,
            "trjconv",
            "-f",
            f"{output}/{prefix}whole.xtc",
            "-s",
            f"{output}/{prefix}whole.pdb",
            "-o",
            f"{output}/{prefix}nojump.pdb",
            "-pbc",
            "nojump",
            "-dump",
            begin,
        ],
        input="System",
        encoding="ascii",
    )

    subprocess.run(
        [
            gmx_path,
            "trjconv",
            "-f",
            f"{output}/{prefix}whole.xtc",
            "-s",
            f"{output}/{prefix}whole.pdb",
            "-o",
            f"{output}/{prefix}nojump.xtc",
            "-pbc",
            "nojump",
        ],
        input="System",
        encoding="ascii",
    )

    subprocess.run(
        [
            gmx_path,
            "trjconv",
            "-f",
            f"{output}/{prefix}nojump.xtc",
            "-s",
            tpr,
            "-o",
            f"{output}/{prefix}center_mol.pdb",
            "-n",
            index,
            "-center",
            "-pbc",
            "mol",
            "-ur",
            "compact",
            "-dump",
            begin,
        ],
        input=f"{selection}\n{selection}",
        encoding="ascii",
    )

    subprocess.run(
        [
            gmx_path,
            "trjconv",
            "-f",
            f"{output}/{prefix}nojump.xtc",
            "-s",
            tpr,
            "-o",
            f"{output}/{prefix}center_mol.xtc",
            "-n",
            index,
            "-center",
            "-pbc",
            "mol",
            "-ur",
            "compact",
        ],
        input=f"{selection}\n{selection}",
        encoding="ascii",
    )

    subprocess.run(
        [
            gmx_path,
            "trjconv",
            "-f",
            f"{output}/{prefix}center_mol.xtc",
            "-s",
            f"{output}/{prefix}center_mol.pdb",
            "-o",
            f"{output}/{prefix}processed.pdb",
            "-fit",
            "progressive",
            "-dump",
            begin,
        ],
        input="System  System",
        encoding="ascii",
    )

    subprocess.run(
        [
            gmx_path,
            "trjconv",
            "-f",
            f"{output}/{prefix}center_mol.xtc",
            "-s",
            f"{output}/{prefix}center_mol.pdb",
            "-o",
            f"{output}/{prefix}processed.xtc",
            "-fit",
            "progressive",
        ],
        input="System  System",
        encoding="ascii",
    )


if __name__ == "__main__":
    import sys

    from md_davis.cli import main as cli_main

    raise SystemExit(cli_main(["process", *sys.argv[1:]]))
