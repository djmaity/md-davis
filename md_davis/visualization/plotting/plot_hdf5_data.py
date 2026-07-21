from __future__ import annotations

from typing import Any, NamedTuple

import h5py
import pandas
import plotly
import plotly.graph_objs as go


class StructureInfo(NamedTuple):
    label: str
    html: str
    color: str


SECONDARY_STRUCTURE = {
    "H": StructureInfo(label="α-helix", html="α helix", color="rgb(255, 0, 0)"),
    "G": StructureInfo(label="3_10-helix", html="3<sub>10</sub> helix", color="rgb(120, 0, 100)"),
    "I": StructureInfo(label="π-helix", html="π helix", color="rgb(255, 15, 139)"),
    "E": StructureInfo(label="β strand", html="β strand", color="rgb(0, 0, 255)"),
    "B": StructureInfo(label="β-bridge", html="β bridge", color="rgb(15, 125, 255)"),
    "T": StructureInfo(label="Turn", html="Turn", color="rgb(0, 200, 0)"),
    "S": StructureInfo(label="Bend", html="Bend", color="rgb(120, 200, 0)"),
    "~": StructureInfo(label="Loop", html="Loop", color="rgb(159, 255, 15)"),
}


def plot_secondary_structure(
    data: pandas.DataFrame | h5py.Dataset,
    name: str,
    figure: go.Figure,
    axis: str | None = None,
    row: int = 1,
    column: int = 1,
    height: float = 100,
    showlegend: bool = True,
) -> None:
    """Plot secondary structure stacked bar graph"""
    if isinstance(data, pandas.DataFrame):
        frames = data.sum(axis=1).max()
    elif isinstance(data, h5py.Dataset):
        frames = sum(data[0])
    else:
        raise TypeError

    for code, structure in SECONDARY_STRUCTURE.items():
        if isinstance(data, pandas.DataFrame):
            x_values = data.index + 1
        if isinstance(data, h5py.Dataset):
            x_values = list(range(1, len(data) + 1))
        y_values = data[code] * (height / frames)
        trace = go.Bar(
            x=x_values,
            y=y_values,
            name=name,
            legendgroup=name,
            showlegend=False,
            marker={
                "color": structure.color,
                "opacity": 0.3,
            },
            hoverinfo="none",
        )
        figure.append_trace(trace, row=row, col=column)
        if axis:
            figure["data"][-1].update(yaxis=axis)
    if showlegend:
        figure["data"][-1]["showlegend"] = True
    figure["layout"]["barmode"] = "stack"


def add_secondary_structure_legend(
    figure: go.Figure,
    xloc: float = 0.8,
    yloc: float = -0.2,
    xspace: float = 0,
    yspace: float = -0.06,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    dssp_legend = []
    dssp_shapes = []
    for index, structure in enumerate(SECONDARY_STRUCTURE.values()):
        dssp_legend.append(
            dict(
                x=xloc + index * xspace + 0.03,
                y=yloc + index * yspace,
                showarrow=False,
                text=structure.html,
                xref="paper",
                yref="paper",
                font=dict(size=16),
                xanchor="left",
                ax=0,
                ay=0,
            )
        )

        dssp_shapes.append(
            dict(
                type="rect",
                x0=xloc + index * xspace,
                y0=yloc + index * yspace + 0.01,
                x1=xloc + index * xspace + 0.02,
                y1=yloc + index * yspace + 0.05,
                xref="paper",
                yref="paper",
                fillcolor=structure.color,
                opacity=0.3,
            )
        )
    return dssp_legend, dssp_shapes


def plot_dihedral_sd(
    dataset: h5py.Dataset,
    name: str,
    figure: go.Figure,
    twin_axis: str,
    row: int = 1,
    column: int = 1,
    height: float = 100,
) -> None:
    """Add dihedral standard deviation into a figure"""
    if twin_axis != "":
        figure["layout"].update({f"yaxis{twin_axis}": dict(anchor="x", overlaying="y", side="left", showgrid=False)})
    for index, structure in enumerate(SECONDARY_STRUCTURE):
        trace = go.Bar(
            x=list(range(len(dataset))),
            y=dataset[structure] / sum(dataset[0]) * height,
            name=name,
            legendgroup=name,
            showlegend=False,
            marker={"color": plotly.colors.DEFAULT_PLOTLY_COLORS[index], "opacity": 0.6},
            hoverinfo="none",
        )
        figure.append_trace(trace, row=row, col=column)
        figure["data"][-1].update(yaxis=f"y{twin_axis}")
    figure["data"][-1]["showlegend"] = True
    figure["layout"]["barmode"] = "stack"
