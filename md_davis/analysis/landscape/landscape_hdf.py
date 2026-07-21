"""Create or plot free-energy landscapes from MD DaVis HDF5 files."""

from __future__ import annotations

from collections.abc import Sequence
from os import PathLike

import h5py
import numpy as np
from numpy.typing import NDArray

from md_davis.analysis.landscape.landscape import Landscape
from md_davis.io.config import parse_literal


def get_sorted_dataset_keys(group: h5py.Group) -> list[str]:
    """Return sorted keys for datasets directly contained in an HDF5 group."""
    return sorted(key for key, value in group.items() if isinstance(value, h5py.Dataset))


def read_timeseries(
    filename: str | PathLike[str], begin: int = 0, end: int | None = None
) -> tuple[str, str, NDArray, NDArray, NDArray]:
    """Read name, label, time, RMSD, and Rg arrays from a collated HDF5 file."""
    with h5py.File(filename, "r") as hdf_file:
        name = hdf_file.attrs["name"]
        label = hdf_file.attrs["label"]
        rmsd_rg = hdf_file["/timeseries/rmsd_rg"]

        if isinstance(rmsd_rg, h5py.Dataset):
            time = rmsd_rg["time"][begin:end]
            rmsd = rmsd_rg["rmsd"][begin:end]
            rg = rmsd_rg["rg"][begin:end]
        else:
            keys = get_sorted_dataset_keys(rmsd_rg)
            time = np.hstack([rmsd_rg[key]["time"] for key in keys])[begin:end]
            rmsd = np.hstack([rmsd_rg[key]["rmsd"] for key in keys])[begin:end]
            rg = np.hstack([rmsd_rg[key]["rg"] for key in keys])[begin:end]

    if not len(time) or not len(rmsd) or not len(rg):
        raise ValueError(f"No timeseries data selected from {filename}")
    return name, label, time, rmsd, rg


def order_landscapes(landscapes: Sequence[Landscape], order: str | None) -> list[Landscape]:
    """Apply an optional name ordering while retaining unspecified landscapes."""
    if not order:
        return landscapes
    requested_names = parse_literal(order, (list, tuple))
    remaining = {landscape.name: landscape for landscape in landscapes}
    ordered = [remaining.pop(name) for name in requested_names if name in remaining]
    return ordered + list(remaining.values())


def main(
    hdf_files: Sequence[str],
    common: bool,
    temperature: float | None,
    shape: tuple[int, int],
    output: str,
    save: str | None,
    begin: int,
    end: int | None,
    limits: str | None,
    precomputed: bool,
    order: str | None,
    width: int | None,
    height: int | None,
    layout: str | None,
    title: str | None,
    axis_labels: str | None,
    orthographic: bool,
    font: str | None,
    font_size: int | None,
    dtick: str | None,
) -> None:
    """Plot free-energy landscapes from collated HDF5 data."""
    parsed_limits = parse_literal(limits, dict)

    if precomputed:
        landscapes = []
        for filename in hdf_files:
            landscapes.extend(Landscape.open(filename).values())
    else:
        input_data = {}
        landscapes = []
        for filename in hdf_files:
            name, label, time, rmsd, rg = read_timeseries(filename, begin=begin, end=end)
            if common:
                input_data[name] = [time, rmsd, rg, label]
            else:
                landscapes.append(
                    Landscape.landscape(
                        name=name,
                        time=time,
                        x_data=rmsd,
                        y_data=rg,
                        shape=shape,
                        label=label,
                        temperature=temperature,
                        limits=parsed_limits,
                    )
                )
        if common:
            landscapes = Landscape.common_landscapes(
                data=input_data,
                shape=shape,
                temperature=temperature,
                limits=parsed_limits,
            )

    landscapes = order_landscapes(landscapes, order)
    if save:
        for landscape in landscapes:
            landscape.save(save, landscape.name, xlabel="RMSD (in nm)", ylabel="Radius of Gyration (in nm)")

    Landscape.plot_landscapes(
        landscapes=landscapes,
        title=title,
        filename=output,
        axis_labels=parse_literal(axis_labels, dict),
        width=width,
        height=height,
        othrographic=orthographic,
        dtick=parse_literal(dtick, dict),
        layout=parse_literal(layout, (list, tuple)),
        font_family=font,
        font_size=font_size,
    )


if __name__ == "__main__":
    import sys

    from md_davis.cli import main as cli_main

    raise SystemExit(cli_main(["landscape", *sys.argv[1:]]))
