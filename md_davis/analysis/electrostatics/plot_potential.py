#! /usr/bin/env python
"""Plot the site potentials file from Delphi"""

from __future__ import annotations

import argparse
from os import PathLike

import plotly.graph_objs as go
from plotly.offline import plot

from md_davis.io.delphi import residue_potential_plot_data


def parse_potential(potential_file: str | PathLike[str]) -> tuple[list[int], list[str], list[float]]:
    return residue_potential_plot_data(potential_file)


def main() -> None:
    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument("pot", help="Input frc file (.pot) from Delphi Calculation", metavar="FILENAME.pot")
    parser.add_argument("-o", "--output", help="Output filename (.html)", default="potential.html", metavar="plot.html")
    args = parser.parse_args()

    residues, labels, potentials = parse_potential(args.pot)

    trace = go.Bar(
        x=residues,
        y=potentials,
        text=labels,
        hoverinfo="text+y",
    )

    layout = go.Layout(
        title="Total Electrostatic Potential on the Protein Surface per Residue",
        xaxis=dict(
            title="Residue Sequence Number",
            zeroline=False,
            gridwidth=2,
        ),
        yaxis=dict(
            title="Electrostatic Potential (kT/e)",
            gridwidth=2,
        ),
        showlegend=False,
    )

    fig = go.Figure(data=[trace], layout=layout)
    plot(fig, filename=args.output)


if __name__ == "__main__":
    main()
