from __future__ import annotations

import operator

import numpy
import pandas
import plotly.graph_objects as go

from md_davis.visualization.contacts import contact_labels

text = contact_labels


def main(
    hbond_file: str,
    total_frames: int | None = None,
    cutoff: float | None = None,
    output: str = "hbond_matrix.html",
    percent: bool = False,
    title: str | None = None,
) -> None:
    """H-bond or contact dataframe obtained from md_davis hbond"""
    df = pandas.read_pickle(hbond_file)

    zmin = 0
    zmax = total_frames
    if percent:
        df["Count"] = df["Count"] * 100 / total_frames  # Convert to percentage of frames
        if cutoff:
            df = df[df["Count"] > cutoff]
            zmin = cutoff
        zmax = 100
    else:
        if cutoff:
            df = df[df["Count"] > cutoff]
            zmin = cutoff

    atoms1 = list(df.groupby(["Chain1", "Residue1", "ResSeq1", "Atom1"]).groups.keys())
    atoms1 = sorted(atoms1, key=operator.itemgetter(0, 2, 1))
    atom_group2 = df.groupby(["Chain2", "Residue2", "ResSeq2", "Atom2"])
    atoms2 = list(atom_group2.groups.keys())
    atoms2 = sorted(atoms2, key=operator.itemgetter(0, 2, 1))

    matrix = numpy.empty(shape=(len(atoms1), len(atoms2)))
    matrix[:] = numpy.nan
    for group, atoms in atom_group2:
        j = atoms2.index(group)
        for _index, row in atoms.iterrows():
            i = atoms1.index((row["Chain1"], row["Residue1"], row["ResSeq1"], row["Atom1"]))
            matrix[i, j] = row["Count"]

    fig = go.Figure(
        data=go.Heatmap(
            z=matrix,
            x=text(atoms1),
            y=text(atoms2),
            colorscale="Greys",
            zmin=zmin,
            zmax=zmax,
            transpose=True,
        )
    )

    if title is not None:
        fig.update_layout(
            title=title,
        )
    fig.write_html(file=output, auto_open=True)


if __name__ == "__main__":
    import sys

    from md_davis.cli import main as cli_main

    raise SystemExit(cli_main(["plot_hbond", *sys.argv[1:]]))
