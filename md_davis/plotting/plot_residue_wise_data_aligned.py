#! /usr/bin/env python
""" Plot the site potentials file from Delphi """

import argparse
from plotly.offline import plot
import plotly.graph_objs as go
from plotly import tools
import pickle
# import sys
# import os
# import numpy
# import re
# from Bio.PDB.PDBParser import PDBParser
# import json

# Imports from within MD Davis package
# import plot_rmsd_rg


def main():
    parser = argparse.ArgumentParser(description='Collect residue wise data into pandas dataframes')
    parser.add_argument('file', help='Pickle file containing data')
    parser.add_argument('-o', '--output', default='output.html', help='Output file name')
    parser.add_argument('-t', '--title', default='', help='Title')
    parser.add_argument('--height', default=1000, action="store",
                        type=int, help='Plot height')
    args = parser.parse_args()

    fig = tools.make_subplots(rows=1, cols=1)
    pickled_data = pickle.load(open(args.file, 'rb'))
    ax_start = len(pickled_data) + 1
    fig['layout'].update({
        f'yaxis{ax_start + 1}': dict(anchor='free', overlaying='y', position=0,
                       side='left', showgrid=False, title='Surface Potential (in kT/e)'),
        f'yaxis{ax_start + 2}': dict(anchor='x', overlaying='y',
                       side='right', showgrid=False, title='Dihedral SD (in degrees)'),
        f'yaxis{ax_start + 3}': dict(anchor='x', overlaying='y',
                       side='left', title='RMSF (in Å)'),
        })

    lcolor = iter(line_color)
    fcolor = iter(fill_color)
    ss_axes = iter(['' if _ < 2 else str(_) for _ in range(1, ax_start)])

    for prefix, residue_data in pickled_data.items():
        for chain, res_data in residue_data['data'].items():
            if 'annotations' in residue_data:
                my_annotations = residue_data['annotations'][chain]
            else:
                my_annotations = None
            pot_traces, dih_traces, rmsf_traces, site_annotations = \
                plot_residue_dataframe.residue_data_trace(
                    figure=fig, data=res_data, prefix=prefix + ' ',
                    annotation=my_annotations,
                    line_color=next(lcolor), fill_color=next(fcolor),
                    ss_axis=next(ss_axes)
                )

        for trace in pot_traces:
            fig.append_trace(trace, 1, 1)
            fig['data'][-1].update(yaxis=f'y{ax_start + 1}')
        for trace in dih_traces:
            fig.append_trace(trace, 1, 1)
            fig['data'][-1].update(yaxis=f'y{ax_start + 2}')
        for trace in rmsf_traces:
            fig.append_trace(trace, 1, 1)
            fig['data'][-1].update(yaxis=f'y{ax_start + 3}')
        if site_annotations:
            for trace in site_annotations:
                fig.append_trace(trace, 1, 1)
                fig['data'][-1].update(yaxis=f'y{ax_start + 3}')

    fig['layout']['font'].update(family='Courier New, monospace', size=16, color='black')

    annotations = plot_hdf5_data.add_secondary_structure_legend(figure=fig, xloc=0, yloc=-0.2, spacing=0.08)
    annotations += [
        dict(x=0.5, y=0.02, showarrow=True, text='Residue Index',
             arrowcolor='rgba(0,0,0,0)',
             xref='paper', yref='paper', font=dict(size=16), ax=0, ay=50),
    ]
    fig['layout']['annotations'] = annotations
    fig['layout'].update({'yaxis': dict(showgrid=False,
                                        ticks='',
                                        showticklabels=False),
                          'xaxis':dict(domain=[0.1, 1])
    })

    fig['layout']['legend'].update(x=1.2, y=0, xanchor='left', yanchor='top')
    # fig['layout']['legend'].update(x=0.5, y=-0.2, xanchor='left', yanchor='top')
    fig['layout'].update(title=args.title,
                         height=args.height,
        margin=go.layout.Margin(l=120, r=100, b=(args.height - 580), t=80, pad=2)
    )

    # print(fig)
    plot(fig, filename=args.output)


if __name__ == "__main__":
    main()
