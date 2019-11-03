#! /usr/bin/env python
""" Plot residue wise data from a pandas dataframe containing one or
    more of the following information:
        1. RMSF
        2. Dihedral circular standard deviation
        3. Surface electrostatic potential
        4. Secondary structure

Usage:  md_davis plot residue [options] FILES...
        md_davis plot residue -h | --help

Options:
  -t, --title <string>              title for the plot
  -h, --help                        display this help and exit
  -o, --output <filename.html>      output HTML file [default: output.html]
"""

import docopt
import pandas
from plotly.offline import plot
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import sys
import os
import numpy
import re
from Bio.PDB.PDBParser import PDBParser
import pickle
import json
import itertools

# Imports from within MD Davis package
from md_davis.plotting import plot_hdf5_data
from md_davis.plotting import plot_rmsd_rg
from md_davis.plotting import colors

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
        "Mutation": dict(
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


def residue_data_trace(figure, data, prefix, ss_axis=None,
                       line_color='rgb(0, 0, 0)',
                       fill_color='rgba(0, 0, 0, 0.4)',
                       annotation=None, row=1, column=1, showlegend=True):
    """ Create traces for the residue wise data """
    mean_potential_traces = []
    total_potential_traces = []
    dihedral_traces = []
    rmsf_traces = []
    seq = data['sequence']
    hover_text = seq.resn + ' ' + seq.resi.map(str)
    # Secondary Structure
    if 'secondary_structure' in data:
        secstr = data['secondary_structure']
        plot_hdf5_data.plot_secondary_structure(
            figure=figure,
            data=secstr,
            name=prefix + 'Secondary Structure',
            axis=ss_axis,
            row=row,
            column=column,
            showlegend=showlegend,
        )
    # Surface Electrostatic Potential
    if 'surface_potential' in data:
        potential = data['surface_potential']
        total_potential_traces += plot_rmsd_rg.continuous_errorbar(
            x=data.index + 1,
            y=potential['mean_total'],
            err=potential['std_total'],
            name=prefix + 'Total Surface Potential',
            hover_text=hover_text,
            line_color=line_color,
            fill_color=fill_color,
            dash='dash',
            showlegend=showlegend,
        )
        mean_potential_traces += plot_rmsd_rg.continuous_errorbar(
            x=data.index + 1,
            y=potential['mean_mean'],
            err=potential['std_mean'],
            name=prefix + 'Mean Surface Potential',
            hover_text=hover_text,
            line_color=line_color,
            fill_color=fill_color,
            dash='dashdot',
            showlegend=showlegend,
        )
        if 'pdb_total' in potential:
            # PDB Surface Potential Trace
            total_potential_traces += [go.Scatter(
                name=prefix + 'PDB Total Surface Potential',
                x=data.index + 1,
                y=potential['pdb_total'],
                mode='markers',
                line=dict(color=line_color),
                text=hover_text,
                hoverinfo='text+y',
                showlegend=False,
            )]
        if 'pdb_mean' in potential:
            mean_potential_traces += [go.Scatter(
                name=prefix + 'PDB Mean Surface Potential',
                x=data.index + 1,
                y=potential['pdb_mean'],
                mode='markers',
                line=dict(color=line_color),
                text=hover_text,
                hoverinfo='text+y',
                showlegend=False,
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
        legendgroup=prefix + ' RMSF',
        showlegend=False,
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
        line=dict(color=line_color, dash='dot'),
        mode='lines',
        text=dih_hover_text[1:-1],
        hoverinfo='text+y',
        legendgroup=prefix + ' Dihedral SD',
        showlegend=False,
    )]
    if annotation:
        sites = annotation_traces(data=data, annotations=annotation,
                                  prefix=prefix, color=line_color)
    else:
        sites = None
    traces = (mean_potential_traces,
              total_potential_traces,
              dihedral_traces,
              rmsf_traces,
              sites,
    )
    return traces


def main(args):
    input_files = args['FILES']
    title = args['--title']
    output_filename = args['--output']

    pickled_data = []
    for input_file in input_files:
        pickled_data.append( pickle.load( open(input_file, 'rb') ) )

    # TODO: Rename 'data' to chains in the dataframe pickle file.
    num_chains = [len(_['data']) for _ in pickled_data]
    max_rows = max(num_chains)
    fig = make_subplots(rows=max_rows, cols=1)

    y_axis = max_rows + 1
    # Axes for Secondary Structure
    if len(num_chains) > 1:
        for _ in range(1, len(num_chains) + 1):
            for ax in range(1, max_rows + 1):
                fig['layout'].update({
                    f'yaxis{y_axis}': dict(anchor=f'x{ax}',
                                           overlaying=f'y{ax}',
                                           side='left',
                                           showgrid=False,
                                           domain=[0.0, 1.0],
                                           showticklabels=False,
                                           ticks='',
                                           range=[0, 100]
                    ),
                })
                y_axis += 1

    # Axes for Residue wise Quantities
    ax_start = y_axis
    my_layout = [('right', 0.9),
                 ('left', 0),
                 ('right', 0.8),
                 ('left', 0.09),
    ]
    for side, pos in my_layout:
        for ax in range(1, max_rows + 1):
            fig['layout'].update({
                f'yaxis{y_axis}': dict(anchor=f'free', overlaying=f'y{ax}',
                            side=side, showgrid=False, position=pos),
                })
            y_axis += 1

    lcolor = itertools.cycle(colors.line_color)
    fcolor = itertools.cycle(colors.fill_color)
    for j, residue_data in enumerate(pickled_data):
        prefix = residue_data['prefix']
        rows = len(residue_data['data'])

        line_color = next(lcolor)
        fill_color = next(fcolor)

        first_chain = True
        for i, (chain, res_data) in enumerate(residue_data['data'].items(), start=1):
            ss_ax = '' if j * max_rows + i < 2 else j * max_rows + i
            axis_number = '' if i < 2 else i

            my_annotations = None
            if 'annotations' in residue_data:
                if chain in residue_data['annotations']:
                    my_annotations = residue_data['annotations'][chain]

            traces = residue_data_trace(
                figure=fig, data=res_data, prefix=prefix + ' ', row=i,
                line_color=line_color, fill_color=fill_color,
                showlegend=first_chain, annotation=my_annotations,
                ss_axis=f'y{ss_ax}'
            )

            mean_pot_traces = traces[0]
            total_pot_traces = traces[1]
            dih_traces = traces[2]
            rmsf_traces = traces[3]
            site_annotations = traces[4]

            if total_pot_traces:
                # if first_chain:
                #     total_pot_traces[0]['showlegend'] = True
                for trace in total_pot_traces:
                    fig.append_trace(trace, i, 1)
                    fig['data'][-1].update(yaxis=f'y{ax_start + i-1}')

            if mean_pot_traces:
                # if first_chain:
                #     mean_pot_traces[0]['showlegend'] = True
                for trace in mean_pot_traces:
                    fig.append_trace(trace, i, 1)
                    fig['data'][-1].update(yaxis=f'y{ax_start + max_rows + i-1}')

            if dih_traces:
                if first_chain:
                    dih_traces[0]['showlegend'] = True
                for trace in dih_traces:
                    fig.append_trace(trace, i, 1)
                    fig['data'][-1].update(yaxis=f'y{ax_start + 2*max_rows + i-1}')

            if rmsf_traces:
                if first_chain:
                    rmsf_traces[0]['showlegend'] = True
                for trace in rmsf_traces:
                    fig.append_trace(trace, i, 1)
                    fig['data'][-1].update(yaxis=f'y{ax_start + 3*max_rows + i-1}')

            if site_annotations:
                for trace in site_annotations:
                    fig.append_trace(trace, i, 1)
                    fig['data'][-1].update(yaxis=f'y{ax_start + 3*max_rows + i-1}')

            fig['layout'].update({
                f'yaxis{axis_number}': dict(showgrid=False,
                                            ticks='',
                                            showticklabels=False,
                                            range=[0, 100],
                                            ),
                f'xaxis{axis_number}':dict(domain=[0.1, 0.8]),
            })
            first_chain = False


    fig['layout']['font'].update(family='Courier New, monospace', size=24, color='black')

    annotations = plot_hdf5_data.add_secondary_structure_legend(figure=fig, spacing=0.06, xloc=0.1 )
    annotations += [
        dict(x=0.0, y=0.5, showarrow=True, text='RMSF (in Å)',
             arrowcolor='rgba(0,0,0,0)', textangle=-90,
             xref='paper', yref='paper', font=dict(size=24), ax=50, ay=0),

        dict(x=-0.1, y=0.5, showarrow=True, text='Mean Surface Potential (in kT/e)',
            arrowcolor='rgba(0,0,0,0)', textangle=-90,
            xref='paper', yref='paper', font=dict(size=24), ax=50, ay=0),

        dict(x=0.81, y=0.5, showarrow=True, text='Dihedral SD (in degrees)',
             arrowcolor='rgba(0,0,0,0)', textangle=-90,
             xref='paper', yref='paper', font=dict(size=24), ax=50, ay=0),

        dict(x=0.95, y=0.5, showarrow=True, text='Total Surface Potential (in kT/e)',
            arrowcolor='rgba(0,0,0,0)', textangle=-90,
            xref='paper', yref='paper', font=dict(size=24), ax=50, ay=0),

        dict(x=0.5, y=-0.03, showarrow=True, text='Residue Index',
             arrowcolor='rgba(0,0,0,0)',
             xref='paper', yref='paper', font=dict(size=24), ax=0, ay=50),
    ]
    fig['layout']['annotations'] = annotations
    fig['layout'].update(title=title,
                         height=1200,
                         width=2000,
        margin=go.layout.Margin(l=100, r=200, b=450, t=80, pad=2)
    )
    plot(fig, filename=output_filename)


if __name__ == "__main__":
    arguments = docopt.docopt(__doc__)
    main(arguments)
