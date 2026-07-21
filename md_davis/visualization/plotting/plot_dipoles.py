#! /usr/env/python

"""
Plot dipole moment from hdf5 files

Usage:
  md_davis plot dipoles [options] [(--3d | --both)] FILES...

Options:
  -t, --title <string>      Title for the plot
  -o, --output FILENAME     Output filename
  -s, --size <int>          Number of points to plot
                            [default: 1000]
  --3d                      Make 3D plot
  -b, --both                Make both a 3D plot and a 2D polar plot
  -c, --center              Show the center in the 3D plot
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from os import PathLike

import h5py
import numpy
import plotly
import plotly.graph_objs as go
from numpy.typing import NDArray
from plotly.colors import DEFAULT_PLOTLY_COLORS
from plotly.offline import plot

# Local imports
from md_davis.analysis import polar


def sample[T](array: Sequence[T] | NDArray, size: int | None = None) -> Sequence[T] | NDArray:
    if not isinstance(array, (list, tuple, numpy.ndarray)):
        raise ValueError
    if size:
        step = int(len(array) / size)
    else:
        step = 1
    return array[::step]


def plot_dipoles(
    dictionary: Mapping[str, NDArray],
    filename: str | PathLike[str] | None = None,
    title: str | None = None,
) -> None:
    fig = plotly.subplots.make_subplots(rows=3, cols=1, vertical_spacing=0.03, shared_xaxes=True)
    for i, (label, value) in enumerate(dictionary.items()):
        magnitude, azimuth, inclination = value
        x_data = list(range(len(magnitude)))

        chosen_color = DEFAULT_PLOTLY_COLORS[i % len(DEFAULT_PLOTLY_COLORS)]
        trace = go.Scatter(
            x=x_data,
            y=magnitude,
            legendgroup=label,
            name=label,
            mode="lines",
            marker=dict(
                color=chosen_color,
            ),
        )
        fig.append_trace(trace, 1, 1)

        trace = go.Scatter(
            x=x_data,
            y=numpy.degrees(azimuth),
            legendgroup=label,
            name=label,
            mode="markers",
            marker=dict(
                size=3,
                color=chosen_color,
            ),
            showlegend=False,
        )
        fig.append_trace(trace, 2, 1)

        trace = go.Scatter(
            x=x_data,
            y=numpy.degrees(inclination),
            legendgroup=label,
            name=label,
            mode="lines",
            marker=dict(
                color=chosen_color,
            ),
            showlegend=False,
        )
        fig.append_trace(trace, 3, 1)

    if title:
        fig["layout"].update(title=title)

    fig["layout"]["xaxis3"].update(range=[0, 1000], title="Time (ns)")
    fig["layout"]["yaxis1"].update(title="Dipole Moment Magnitude (Debye)")
    fig["layout"]["yaxis2"].update(title="Azimuth (degrees)", range=[-181, 181], tick0=-180, dtick=60)
    fig["layout"]["yaxis3"].update(title="Inclination (degrees)", range=[0, 181], tick0=0, dtick=30)
    if filename:
        plot(fig, filename=filename)
    else:
        fig.show()


def plot_dipoles3d(
    dictionary: Mapping[str, NDArray],
    filename: str | PathLike[str] | None = None,
    title: str | None = None,
    centroid: bool = False,
) -> None:
    data = []
    for i, (label, value) in enumerate(dictionary.items()):
        x_data, y_data, z_data = value
        chosen_color = DEFAULT_PLOTLY_COLORS[i % len(DEFAULT_PLOTLY_COLORS)]
        trace = go.Scatter3d(
            x=x_data,
            y=y_data,
            z=z_data,
            name=label,
            legendgroup=label,
            marker=dict(
                color=chosen_color,
                size=2,
                opacity=0.2,
            ),
            mode="markers",
        )
        data.append(trace)
        if centroid:
            center = numpy.mean(value, axis=1)
            centroid_trace = go.Scatter3d(
                x=[center[0]],
                y=[center[1]],
                z=[center[2]],
                name=label,
                legendgroup=label,
                marker=dict(
                    color=chosen_color,
                ),
                mode="markers",
            )
            data.append(centroid_trace)

    fig = go.Figure(data=data)
    if title:
        fig["layout"].update(title=title)
    fig.update_layout(
        scene=dict(
            xaxis=dict(
                gridcolor="white",
                showbackground=True,
                zerolinecolor="black",
            ),
            yaxis=dict(gridcolor="white", showbackground=True, zerolinecolor="black"),
            zaxis=dict(
                gridcolor="white",
                showbackground=True,
                zerolinecolor="black",
            ),
        ),
    )

    if filename:
        plot(fig, filename=filename)
    else:
        fig.show()


def main(
    files: Sequence[str],
    size: int,
    output: str | None,
    title: str | None,
    three_dimensional: bool,
    both: bool,
    center: bool,
) -> None:
    """Plot dipole moment from HDF5 files"""

    polar_dipoles = {}
    dipoles_dict = {}
    include_2d = not three_dimensional or both
    include_3d = three_dimensional or both
    for fname in files:
        with h5py.File(fname, "r") as hdf_file:
            data = numpy.asarray(hdf_file["dipole_moment"])
            dipoles = numpy.vstack(
                [
                    sample(data["mu_x"], size=size),
                    sample(data["mu_y"], size=size),
                    sample(data["mu_z"], size=size),
                ]
            )
            label = hdf_file.attrs["short_html"]
            if include_2d:
                polar_dipoles[label] = polar.spherical_numpy(xyz=dipoles.T).T
            if include_3d:
                dipoles_dict[label] = dipoles

    if include_2d:
        plot_dipoles(polar_dipoles, filename=output, title=title)
    if include_3d:
        three_dimensional_output = f"3D_{output}" if both and output else output
        plot_dipoles3d(dipoles_dict, filename=three_dimensional_output, title=title, centroid=center)


if __name__ == "__main__":
    import sys

    from md_davis.cli import main as cli_main

    raise SystemExit(cli_main(["plot-dipoles", *sys.argv[1:]]))
