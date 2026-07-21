#!/usr/bin/env python3

from __future__ import annotations

import colorsys
from collections.abc import Iterator, Sequence

import h5py
import numpy
import plotly
import plotly.graph_objs as go
from numpy.typing import ArrayLike
from plotly.offline import plot

from md_davis.analysis.statistics import moving_avg, moving_avg_std


def get_colors(
    num_colors: int, lightness: float = 0.5, saturation: float = 0.5, opacity: float = 1.0
) -> Iterator[str]:
    for hue in numpy.linspace(0, 1, num_colors, endpoint=False):
        red, green, blue = colorsys.hls_to_rgb(h=hue, l=lightness, s=saturation)
        yield f"rgba({int(red * 255)}, {int(green * 255)}, {int(blue * 255)}, {opacity})"


def continuous_errorbar(
    x: ArrayLike,
    y: ArrayLike,
    err: ArrayLike,
    name: str,
    hover_text: Sequence[str] | None = None,
    line_color: str | None = None,
    fill_color: str | None = None,
    dash: str | None = None,
    showlegend: bool = True,
) -> list[go.Scatter]:
    """return continuous errorbar plotly trace"""
    if not line_color:
        line_color = "rgb(31,119,180)"
    if not fill_color:
        fill_color = "rgba(31,119,180,0.3)"

    upper_bound = go.Scatter(
        name=name,
        x=x,
        y=y + err,
        mode="lines",
        line=dict(width=0),
        legendgroup=name,
        showlegend=False,
        hoverinfo="none",
        fillcolor=fill_color,
        fill="tonexty",
    )
    trace = go.Scatter(
        name=name,
        x=x,
        y=y,
        mode="lines",
        line=dict(color=line_color, dash=dash),
        showlegend=showlegend,
        legendgroup=name,
        text=hover_text,
        hoverinfo="text+y",
        fillcolor=fill_color,
        fill="tonexty",
    )
    lower_bound = go.Scatter(
        name=name,
        x=x,
        y=y - err,
        mode="lines",
        line=dict(width=0),
        legendgroup=name,
        showlegend=False,
        hoverinfo="none",
    )
    # Trace order can be important
    # with continuous error bars
    return [lower_bound, trace, upper_bound]


def timeseries_plot(
    files: Sequence[str], output: str, title: str, window: int, dataset: str
) -> None:
    fig = plotly.subplots.make_subplots(rows=1, cols=1)
    fig["layout"].update(
        {"yaxis2": dict(anchor="x", overlaying="y", side="right", showgrid=False, title="R<sub>G</sub>")}
    )
    rmsd_traces, rg_traces = [], []

    num_files = len(files)
    rmsd_color = get_colors(num_files, lightness=0.5, saturation=0.8)
    rmsd_fill = get_colors(num_files, lightness=0.5, saturation=0.8, opacity=0.3)
    rg_color = get_colors(num_files, lightness=0.6, saturation=0.3)
    rg_fill = get_colors(num_files, lightness=0.6, saturation=0.3, opacity=0.3)

    for fname in files:
        print(fname)
        with h5py.File(fname, "r") as datafile:
            dset = datafile[dataset][1:]

            time = moving_avg(dset["time"], window=window, shift=window)
            rmsd, rmsd_err = moving_avg_std(dset["rmsd"], window=window)
            rg, rg_err = moving_avg_std(dset["rg"], window=window)

            rmsd_traces += continuous_errorbar(
                x=time,
                y=rmsd,
                err=rmsd_err,
                name="RMSD " + datafile.attrs["short_html"],
                line_color=next(rmsd_color),
                fill_color=next(rmsd_fill),
            )
            rg_traces += continuous_errorbar(
                x=time,
                y=rg,
                err=rg_err,
                name="R<sub>G</sub> " + datafile.attrs["short_html"],
                line_color=next(rg_color),
                fill_color=next(rg_fill),
            )
    for trace in rmsd_traces:
        fig.append_trace(trace, 1, 1)
    for trace in rg_traces:
        fig.append_trace(trace, 1, 1)
        fig["data"][-1].update(yaxis="y2")

    fig["layout"]["legend"].update(x=1.05, y=1)
    fig["layout"]["xaxis1"].update(title="Time (in ns)")
    fig["layout"]["yaxis1"].update(title="RMSD (in Å)")
    fig["layout"].update(title=title)
    plot(fig, filename=output)


if __name__ == "__main__":
    import sys

    from md_davis.cli import main as cli_main

    raise SystemExit(cli_main(["timeseries", *sys.argv[1:]]))
