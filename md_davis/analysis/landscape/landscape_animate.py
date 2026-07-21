#! /usr/bin/env python
"""
Create an HDF file for collection analysis data

Usage:
  md_davis landscape animation [options] HDF_FILE
  md_davis landscape animation -h | --help

Options:
  -h, --help                    Show this screen.
  -o, --output <filename.html>  Name for the output HTML file containing
                                the plots [default: landscapes.html]
  -i, --index int               Index of the landscape to plot [default: 0]
  -t, --title <string>          Title for the figure [default: Energy Landscape]
  -b, --begin int               Starting index for the data to include
                                [default: 0]
  -e, --end int                 Last index for the data to include
  -s, --step int                Step size for the animation data [default: 1]
  --hide_surface                Hide the energy landscape surface
  --ortho                       Orthographic projection for 3D plots
  --width <int>                 Width of the plot
  --height <int>                Height of the plot
  --font_size <int>             Size of text and labels in the plot
  --dtick <dict>                Tick interval on each axes
  --hide_labels                 Hide the axes labels
  --static                      Make a static plot with the time points colored
"""

from __future__ import annotations

from collections.abc import Mapping
from os import PathLike
from typing import Any

import numpy
import plotly
import plotly.graph_objects as go
from numpy.typing import NDArray

from md_davis.analysis.landscape.landscape import Landscape
from md_davis.io.config import parse_literal


def get_animation_data(landscape: Landscape) -> list[tuple[float, float, float, float]]:
    output = []
    for i, row in enumerate(landscape.binned_data):
        for j, time in enumerate(row):
            for t in time:
                output.append((t, landscape.xbins[i], landscape.ybins[j], landscape.zValues[i, j]))
    return sorted(output)


def frame_args(duration: int) -> dict[str, Any]:
    return {
        "frame": {"duration": duration},
        "mode": "immediate",
        "fromcurrent": True,
        "transition": {"duration": duration, "easing": "linear"},
    }


def landscape_trajectory(
    landscape: Landscape,
    data: NDArray,
    filename: str | PathLike[str] = "landscape.html",
    save_image: str | PathLike[str] | None = None,
    xlabel: str = "x",
    ylabel: str = "y",
    zlabel: str = "z",
    width: int | None = None,
    height: int | None = None,
    hide_surface: bool = False,
    marker_size: int | None = None,
    title: str | None = None,
    camera: Mapping[str, Any] | None = None,
    othrographic: bool = False,
    font_family: str | None = None,
    font_size: int | None = None,
    dtick: Mapping[str, float] | None = None,
    hide_ticks: bool = False,
) -> None:
    """Show the X and Y values from the trajectory on the landscape"""
    x, y = numpy.meshgrid(landscape.xbins, landscape.ybins, indexing="ij")
    z = landscape.zValues
    time, x_val, y_val, z_val = data.T
    data = []
    if not hide_surface:
        trace = go.Surface(
            x=x,
            y=y,
            z=z,
            colorscale="Cividis",
            cmin=landscape.limits["z"][0],
            cmax=landscape.limits["z"][1],
            contours_z=dict(show=True, usecolormap=True, highlightcolor="limegreen", project_z=True),
        )
        data.append(trace)
    trace = go.Scatter3d(
        x=x_val,
        name="Trajectory",
        y=y_val,
        z=z_val,
        mode="markers",
        marker=dict(
            color=time,
            reversescale=True,
            colorbar=dict(x=0.9, title="Time"),
            size=marker_size,
        ),
    )
    data.append(trace)
    fig = go.Figure(data=data)
    for axis in ["x", "y", "z"]:
        if axis in landscape.limits:
            fig["layout"]["scene"][axis + "axis"].update(range=landscape.limits[axis])

    if dtick:
        for axis in ["x", "y", "z"]:
            fig["layout"]["scene"][axis + "axis"].update(dtick=dtick[axis])

    if hide_ticks:
        for axis in ["x", "y", "z"]:
            fig["layout"]["scene"][axis + "axis"].update(showticklabels=False)

    fig["layout"]["scene"].update(
        xaxis_title=xlabel,
        yaxis_title=ylabel,
        zaxis_title=zlabel,
    )

    if camera is None:
        camera = dict(eye=dict(x=1.25, y=1.25, z=1.25))
    fig.update_layout(
        title=title if title else landscape.label,
        showlegend=True,
        legend_orientation="h",
        scene_camera=camera,
    )
    if othrographic:
        fig.layout.scene.camera.projection.update(type="orthographic")
    if width:
        fig.layout.update(width=width)
    if height:
        fig.layout.update(height=height)
    if font_size:
        fig.layout.font.update(
            family=font_family,
            size=font_size - 10 if font_size > 20 else 12,
            # Fonts size for ticks do not match other fonts
        )
    if save_image:
        fig.write_image(
            save_image,
            engine="kaleido",
            width=width,
            height=height,
        )
    fig.write_html(filename)


def landscape_animation(
    landscape: Landscape,
    data: NDArray,
    filename: str | PathLike[str] = "landscape.html",
    xlabel: str = "x",
    ylabel: str = "y",
    zlabel: str = "z",
    title: str | None = None,
    width: int | None = None,
    height: int | None = None,
    othrographic: bool = False,
    font_size: int | None = None,
    dtick: Mapping[str, float] | None = None,
) -> None:
    x, y = numpy.meshgrid(landscape.xbins, landscape.ybins, indexing="ij")
    z = landscape.zValues
    time, x_val, y_val, z_val = data.T
    fig = go.Figure(
        data=[
            go.Scatter3d(
                x=x_val,
                y=y_val,
                z=z_val,
            ),
            go.Surface(
                x=x,
                y=y,
                z=z,
                colorscale="Cividis",
                cmin=landscape.limits["z"][0],
                cmax=landscape.limits["z"][1],
                contours_z=dict(show=True, usecolormap=True, highlightcolor="limegreen", project_z=True),
            ),
        ],
        frames=[
            go.Frame(
                data=[go.Scatter3d(x=[x], y=[y], z=[z], mode="markers", marker=dict(color="red", size=10))],
                name=t,
            )
            for t, x, y, z in data
        ],
    )
    for axis in ["x", "y", "z"]:
        if axis in landscape.limits:
            fig["layout"]["scene"][axis + "axis"].update(range=landscape.limits[axis])

    sliders = [
        {
            "active": 0,
            "yanchor": "top",
            "xanchor": "left",
            "currentvalue": {"font": {"size": 20}, "prefix": "Time:", "visible": True, "xanchor": "right"},
            "transition": {"duration": 300, "easing": "cubic-in-out"},
            "pad": {"b": 10, "t": 50},
            "len": 0.9,
            "x": 0.1,
            "y": 0,
            "steps": [
                {
                    "args": [[f.name], frame_args(50)],
                    "label": str(k),
                    "method": "animate",
                }
                for k, f in enumerate(fig.frames)
            ],
        }
    ]
    fig.update_layout(
        title=title if title else landscape.label,
        updatemenus=[
            {
                "buttons": [
                    {
                        "args": [None, frame_args(50)],
                        "label": "&#9654;",  # play symbol
                        "method": "animate",
                    },
                    {
                        "args": [[None], frame_args(0)],
                        "label": "&#9724;",  # pause symbol
                        "method": "animate",
                    },
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 70},
                "type": "buttons",
                "x": 0.1,
                "y": 0,
            }
        ],
        sliders=sliders,
    )
    plotly.offline.plot(fig, filename=filename)


def main(
    hdf_file: str,
    output: str,
    static: bool = True,
    select: str | None = None,
    begin: int = 0,
    end: int | None = None,
    step: int = 1,
    title: str | None = None,
    orthographic: bool = False,
    axis_labels: Mapping[str, str] | str | None = None,
    width: int | None = None,
    height: int | None = None,
    font: str | None = None,
    marker_size: int | None = None,
    font_size: int | None = None,
    dtick: Mapping[str, float] | str | None = None,
    hide_surface: bool = False,
    camera: Mapping[str, Any] | str | None = None,
) -> None:
    """ """
    landscapes = Landscape.open(hdf_file, select=select)

    if select is not None:
        if select not in landscapes:
            print("selection not found in landscapes")
            return
        else:
            landscape = landscapes[select]
    else:
        landscape = next(iter(landscapes.values()))

    data = get_animation_data(landscape)[begin:end:step]
    data = numpy.array(data)

    if axis_labels is None:
        axis_labels = dict(x=" <br>RMSD (in nm)", y=" <br>Rg (in nm)", z="Free Energy (kJ mol<sup>-1</sup>)<br> ")
    else:
        if isinstance(axis_labels, str):
            axis_labels = parse_literal(axis_labels, dict)

    if static:
        landscape_trajectory(
            landscape=landscape,
            data=data,
            hide_surface=hide_surface,
            title=title,
            filename=output,
            xlabel=axis_labels["x"],
            ylabel=axis_labels["y"],
            zlabel=axis_labels["z"],
            width=width,
            height=height,
            font_size=font_size,
            dtick=dtick,
            othrographic=orthographic,
            camera=parse_literal(camera, dict),
            marker_size=marker_size,
        )
    else:
        landscape_animation(
            landscape=landscape,
            data=data,
            title=title,
            filename=output,
            xlabel=axis_labels["x"],
            ylabel=axis_labels["y"],
            zlabel=axis_labels["z"],
            width=width,
            height=height,
            font_size=font_size,
            dtick=dtick,
            othrographic=orthographic,
        )


if __name__ == "__main__":
    import sys

    from md_davis.cli import main as cli_main

    raise SystemExit(cli_main(["landscape_animate", *sys.argv[1:]]))
