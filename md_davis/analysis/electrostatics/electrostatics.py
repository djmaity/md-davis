"""Calculate and aggregate electrostatic potentials with MSMS and Delphi."""

from __future__ import annotations

import subprocess
import warnings
from collections.abc import Sequence
from os import PathLike
from pathlib import Path

import pandas as pd

from md_davis.analysis.electrostatics.vert2pdb import vertices_to_pdb
from md_davis.io.delphi import read_delphi_potential


def run_msms(
    pdb_file: str | PathLike[str],
    output_directory: str | PathLike[str] = ".",
    msms: str | PathLike[str] | None = None,
) -> str:
    """Calculate an MSMS triangulated surface and return its vertex filename."""
    pdb_path = Path(pdb_file).resolve()
    output_path = Path(output_directory).resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    msms_path = Path(msms).resolve()
    xyz_path = output_path / f"{pdb_path.stem}.xyz"

    with xyz_path.open("w", encoding="utf-8") as xyz_file:
        subprocess.run(
            [str(msms_path.with_name("pdb_to_xyzrn")), str(pdb_path)],
            stdout=xyz_file,
            cwd=msms_path.parent,
            check=True,
        )

    subprocess.run(
        [str(msms_path), "-if", xyz_path.name, "-of", pdb_path.stem, "-probe_radius", "1.4"],
        cwd=output_path,
        check=True,
    )
    return str(output_path / f"{pdb_path.stem}.vert")


def run_delphi(
    pdb_file: str | PathLike[str],
    output_directory: str | PathLike[str],
    output_filename: str,
    delphi_path: str | PathLike[str],
    radius_file: str | PathLike[str] | None = None,
    charge_file: str | PathLike[str] | None = None,
    pqr: bool = False,
    grid_size: int = 101,
    surface: str | PathLike[str] | None = None,
    center: bool = False,
) -> None:
    """Run Delphi and optionally calculate potential at supplied surface points."""
    output_path = Path(output_directory).resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    pdb_path = Path(pdb_file).resolve()
    cube_path = output_path / f"{output_filename}.cube"

    parameters = [
        f"gsize={grid_size}",
        "salt=0.15",
        "exdi=80",
        "linit=2000",
        "maxc=0.0000000001",
        f'out(phi, file="{cube_path}", format="cube")',
        f'in(pdb,file="{pdb_path}")',
    ]

    if not pqr:
        data_directory = Path(__file__).parent
        radius_path = Path(radius_file).resolve() if radius_file else data_directory / "charmm.siz"
        charge_path = Path(charge_file).resolve() if charge_file else data_directory / "charmm.crg"
        parameters.extend(
            [
                f'in(siz,file="{radius_path}")',
                f'in(crg,file="{charge_path}")',
            ]
        )

    if center:
        parameters.append("acenter(0,0,0)")
    if surface:
        surface_path = Path(surface).resolve()
        potential_path = output_path / f"{output_filename}.pot"
        parameters.extend(
            [
                f'in(frc,file="{surface_path}")',
                f'out(frc, file="{potential_path}")',
                "site(Atom, Potential, Reaction, Coulomb, Field)",
            ]
        )

    parameter_file = output_path / f"{output_filename}_tmp.prm"
    parameter_file.write_text("\n".join(parameters) + "\n", encoding="utf-8")
    try:
        subprocess.run(
            [str(Path(delphi_path).resolve()), parameter_file.name],
            cwd=output_path,
            check=True,
        )
    finally:
        parameter_file.unlink(missing_ok=True)


def parse_electrostatic_potential(
    potential_file: str | PathLike[str], atomic_potentials: bool = False
) -> dict[str, pd.DataFrame]:
    """Aggregate Delphi atom potentials into per-chain residue dataframes."""
    dataframe = read_delphi_potential(potential_file)
    if dataframe["chainID"].isnull().values.any():
        warnings.warn(
            "Chain ID missing; multichain results may be incorrect",
            stacklevel=2,
        )

    dataframe["chainID"] = dataframe["chainID"].fillna("A")
    output = {}
    for chain_index, (_, chain_data) in enumerate(dataframe.sort_values(by="chainID").groupby("chainID")):
        group_columns = ["resSeq", "name"] if atomic_potentials else ["resSeq"]
        grouped = chain_data.groupby(group_columns, as_index=False)["potential"]
        potential = grouped.sum().rename(columns={"potential": "total"})
        potential["mean"] = grouped.mean()["potential"]
        output[f"chain {chain_index}"] = potential
    return output


def main(
    pdb_files: Sequence[str],
    output: str,
    msms: str,
    delphi: str,
    radius: str | None,
    charge: str | None,
    grid_size: int,
    surface: bool,
    center: bool,
) -> None:
    """Calculate electrostatic potential maps for PDB structures."""
    output_path = Path(output)
    for pdb_file in pdb_files:
        output_filename = Path(pdb_file).stem
        surface_file = None
        if surface:
            vertex_file = run_msms(pdb_file, output_directory=output, msms=msms)
            surface_file = output_path / f"{output_filename}_surf.pdb"
            with surface_file.open("w", encoding="utf-8") as output_stream:
                vertices_to_pdb(vertex_file, pdb_file, output_stream)

        run_delphi(
            pdb_file=pdb_file,
            output_directory=output,
            output_filename=output_filename,
            surface=surface_file,
            delphi_path=delphi,
            radius_file=radius,
            charge_file=charge,
            grid_size=grid_size,
            center=center,
        )


if __name__ == "__main__":
    import sys

    from md_davis.cli import main as cli_main

    raise SystemExit(cli_main(["electrostatics", *sys.argv[1:]]))
