"""Parsers for DelPhi site-potential output."""

from __future__ import annotations

from os import PathLike

import pandas as pd

POTENTIAL_COLUMNS = [
    "name",
    "resName",
    "chainID",
    "resSeq",
    "potential",
    "reaction",
    "coulomb",
    "Ex",
    "Ey",
    "Ez",
]
POTENTIAL_WIDTHS = [5, 3, 3, 9, 10, 10, 10, 10, 10, 10]


def read_delphi_potential(filename: str | PathLike[str]) -> pd.DataFrame:
    """Read the atomic rows of a DelPhi ``.pot`` file."""
    return pd.read_fwf(
        filename,
        skiprows=12,
        skipfooter=2,
        dtype={"resSeq": int},
        engine="python",
        names=POTENTIAL_COLUMNS,
        widths=POTENTIAL_WIDTHS,
    )


def residue_potential_series(filename: str | PathLike[str]) -> pd.Series:
    """Sum atomic potentials for each residue index."""
    data = read_delphi_potential(filename)
    return data.groupby("resSeq")["potential"].sum()


def residue_potential_plot_data(
    filename: str | PathLike[str],
) -> tuple[list[int], list[str], list[float]]:
    """Return residue indices, labels, and summed potentials for plotting."""
    data = read_delphi_potential(filename)
    grouped = data.groupby(["resSeq", "resName"])["potential"].sum()
    residues = [int(residue) for residue, _ in grouped.index]
    labels = [str(index) for index in grouped.index]
    return residues, labels, grouped.astype(float).tolist()
