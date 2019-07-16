#! /usr/bin/env python
""" Plot the site potentials file from Delphi """

import argparse
import pandas
from plotly.offline import plot
import plotly.graph_objs as go
from plotly import tools
import sys
import os
import numpy
import re
from Bio.PDB.PDBParser import PDBParser
import pickle
import json

# Imports from within MD Davis package
import plot_hdf5_data
import plot_rmsd_rg

def marker_type(trace_color):
    marker_dict = {
        "Active Site": dict(
            color='rgba(0, 0, 0, 0)',
            size=10,
            line=dict(color=trace_color, width=3)
        ),
        "Nucleotide Binding Regions": dict(
            symbol='square',
            color=trace_color,
            size=6,
        ),
        "NADP Binding Site": dict(
            symbol='square',
            color='rgba(0, 0, 0, 0)',
            size=10,
            line=dict(color=trace_color, width=3)
        ),
        "Substrate Binding Site": dict(
            color='rgba(0, 0, 0, 0)',
            size=10,
            line=dict(color=trace_color, width=3)
        ),
        "Metal Binding Site": dict(
            symbol='circle-dot',
            color=trace_color,
            size=10,
        ),
        "Cofactor Binding Site": dict(
            symbol='diamond',
            color='rgba(0, 0, 0, 0)',
            line=dict(color=trace_color, width=3),
            size=10,
        ),
    }
    return marker_dict


def alternate_join(list1, list2):
    """ Combine two lists in with consecutive elements from alternate lists """
    result = [None]*(len(list1)+len(list2))
    result[::2] = list1
    result[1::2] = list2
    return result


def annotation_traces(data, annotations, prefix, color='rgb(0,0,0)'):
    sites = []
    for text, res_list in annotations.items():
        if text in marker_type(color):
            use_marker = marker_type(color)[text]
        else:
            use_marker = dict(
                              symbol='cross',
                              color=color,
                              size=10,
                              )

        residues = []
        for res in res_list:
            if isinstance(res, list):
                assert len(res) == 2
                residues += list(range(res[0], res[1] + 1))
            else:
                residues.append(res)

        site = data[ data.sequence.resi.isin(residues) ]
        sites += [go.Scatter(
            name=prefix + f'{text}',
            x=site.index + 1,
            y=site.rmsf.values.flatten() * 10,
            mode='markers',
            hoverinfo='none',
            marker = use_marker
            ),
        ]
    return sites


def residue_data_trace(figure, data, prefix,
                       line_color='rgb(0, 0, 0)',
                       fill_color='rgba(0, 0, 0, 0.4)',
                       annotation=None, ss_axis=''):
    """ Create traces for the residue wise data """
    potential_traces, dihedral_traces, rmsf_traces, = [], [], []
    seq = data['sequence']
    hover_text = seq.resn + ' ' + seq.resi.map(str)

    # Secondary Structure
    secstr = data['secondary_structure']
    plot_hdf5_data.plot_secondary_structure(data=secstr,
        name=prefix + 'Secondary Structure',
        figure=figure,
        twin_axis=ss_axis,
    )
    potential = data['surface_potential']
    potential_traces += plot_rmsd_rg.continuous_errorbar(
        x=data.index + 1,
        y=potential['mean_total'],
        err=potential['std_total'],
        name=prefix + 'Total Surface Potential',
        hover_text=hover_text,
        line_color=line_color,
        fill_color=fill_color,
        dash='dot',
    )
    potential_traces += plot_rmsd_rg.continuous_errorbar(
        x=data.index + 1,
        y=potential['mean_mean'],
        err=potential['std_mean'],
        name=prefix + 'Mean Surface Potential',
        hover_text=hover_text,
        line_color=line_color,
        fill_color=fill_color,
        dash='dot',
    )
    # PDB Surface Potential Trace
    potential_traces += [go.Scatter(
        name=prefix + 'PDB Total Surface Potential',
        x=data.index + 1,
        y=potential['pdb_total'],
        mode='markers',
        line=dict(color=line_color),
        text=hover_text,
        hoverinfo='text+y',
        # legendgroup=prefix,
    )]
    potential_traces += [go.Scatter(
        name=prefix + 'PDB Mean Surface Potential',
        x=data.index + 1,
        y=potential['pdb_mean'],
        mode='markers',
        line=dict(color=line_color),
        text=hover_text,
        hoverinfo='text+y',
        # legendgroup=prefix,
    )]
    # RMSF
    rmsf_traces += [go.Scatter(
        name=prefix + 'RMSF',
        x=data.index + 1,
        y=data.rmsf.values.flatten() * 10,
        line=dict(color=line_color),
        text=hover_text,
        hoverinfo='text+y',
        mode='lines',
        # legendgroup=prefix,
    )]
    # Dihedral SD
    x_values = alternate_join(numpy.array(data.index) + 0.75,
                    numpy.array(data.index) + 1.25)
    dih_sd = alternate_join(data['dihedral_standard_deviation'].phi,
                    data['dihedral_standard_deviation'].psi)
    dih_hover_text = alternate_join(seq.resn + ' ϕ ' + seq.resi.map(str), 
                                    seq.resn + ' ψ ' + seq.resi.map(str))
    dihedral_traces += [go.Scatter(
        name=prefix + 'Dihedral Standard Deviation',
        x=x_values[1:-1],   # Omit the first and last angle which do not exist
        y=numpy.degrees(dih_sd[1:-1]),
        line=dict(color=line_color, dash='dash'),
        mode='lines',
        text=dih_hover_text[1:-1],
        hoverinfo='text+y',
        # legendgroup=prefix,
    )]
    if annotation:
        sites = annotation_traces(data=data, annotations=annotation,
                                  prefix=prefix, color=line_color)
    else:
        sites = None

    return potential_traces, dihedral_traces, rmsf_traces, sites


def main():
    parser = argparse.ArgumentParser(description='Collect residue wise data into pandas dataframes')
    parser.add_argument('file', help='Pickle file containing data')
    parser.add_argument('-o', '--output', default='output.html', help='Output file name')
    parser.add_argument('-t', '--title', default='', help='Title')
    args = parser.parse_args()

    residue_data = pickle.load(open(args.file, 'rb'))
    prefix = residue_data['prefix']

    fig = tools.make_subplots(rows=1, cols=1)
    fig['layout'].update({
        'yaxis2': dict(anchor='free', overlaying='y', position=0,
                       side='left', showgrid=False, title='Surface Potential (in kT/e)'),
        'yaxis3': dict(anchor='x', overlaying='y',
                       side='right', showgrid=False, title='Dihedral SD (in degrees)'),
        'yaxis4': dict(anchor='x', overlaying='y',
                       side='left', title='RMSF (in Å)'),
        })

    for chain, res_data in residue_data['data'].items():
        # print(chain)
        if 'annotations' in residue_data:
            pot_traces, dih_traces, rmsf_traces, site_annotations = residue_data_trace(
                figure=fig, data=res_data, prefix=prefix + ' ',
                annotation=residue_data['annotations'][chain])
        else:
            pot_traces, dih_traces, rmsf_traces, site_annotations = residue_data_trace(
                figure=fig, data=res_data, prefix=prefix + ' ')
            
    for trace in pot_traces:
        fig.append_trace(trace, 1, 1)
        fig['data'][-1].update(yaxis=f'y2')
    for trace in dih_traces:
        fig.append_trace(trace, 1, 1)
        fig['data'][-1].update(yaxis=f'y3')
    for trace in rmsf_traces:
        fig.append_trace(trace, 1, 1)
        fig['data'][-1].update(yaxis=f'y4')

    if site_annotations:
        for trace in site_annotations:
            fig.append_trace(trace, 1, 1)
            fig['data'][-1].update(yaxis=f'y4')

    fig['layout']['font'].update(family='Courier New, monospace', size=24, color='black')

    annotations = plot_hdf5_data.add_secondary_structure_legend(figure=fig, spacing=0.08)
    annotations += [
        dict(x=0.5, y=-0.05, showarrow=True, text='Residue Index',
             arrowcolor='rgba(0,0,0,0)',
             xref='paper', yref='paper', font=dict(size=24), ax=0, ay=50),
    ]
    fig['layout']['annotations'] = annotations
    fig['layout'].update({'yaxis': dict(showgrid=False,
                                        ticks='',
                                        showticklabels=False),
                          'xaxis':dict(domain=[0.1, 1])
    })
    fig['layout']['legend'].update(x=0, y=-0.2, xanchor='left', yanchor='top')
    fig['layout'].update(title=args.title,
                         height=1000,
        margin=go.layout.Margin(l=120, r=100, b=350, t=80, pad=2)
    )
    plot(fig, filename=args.output)


if __name__ == "__main__":
    main()
