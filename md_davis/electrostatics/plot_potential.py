#! /usr/bin/env python
""" Plot the site potentials file from Delphi """

import argparse
import pandas
from plotly.offline import plot
import plotly.graph_objs as go


def parse_potential(potential_file):
    df = pandas.read_csv(potential_file, skiprows=12, delim_whitespace=True, skipfooter=2,
        dtype={'resSeq': int}, engine='python',
        names=['name', 'resName', 'chainID', 'resSeq', 'potential',
            'reaction', 'coulomb', 'Ex', 'Ey', 'Ez'],
    )
    potential_dict =  df.groupby(['resSeq', 'resName'])['potential'].sum().to_frame().to_dict()

    resi, labels, pot = [], [], []
    for res, potential in potential_dict['potential'].items():
        resSeq, _ = res
        labels.append(str(res))
        resi.append(resSeq)
        pot.append(potential)
    return resi, labels, pot


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('pot', help='Input frc file (.pot) from Delphi Calculation',
                        metavar='FILENAME.pot')
    parser.add_argument('-o', '--output', help='Output filename (.html)',
                        default='potential.html', metavar='plot.html')
    args = parser.parse_args()

    residues, labels, potentials = parse_potential(args.pot)

    trace = go.Bar(
        x = residues,
        y = potentials,
        text=labels,
        hoverinfo='text+y',
    )

    layout= go.Layout(
        title= 'Total Electrostatic Potential on the Protein Surface per Residue',
        xaxis= dict(
            title= 'Residue Sequence Number',
            zeroline= False,
            gridwidth= 2,
        ),
        yaxis=dict(
            title= 'Electrostatic Potential (kT/e)',
            gridwidth= 2,
        ),
        showlegend= False
    )


    fig = go.Figure(data=[trace], layout=layout)
    plot(fig, filename=args.output)


if __name__ == "__main__":
    main()