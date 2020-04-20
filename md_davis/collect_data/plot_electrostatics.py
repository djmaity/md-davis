#! /usr/bin/env python
""" Plot the site potentials file from Delphi """

import argparse
import pandas
from plotly.offline import plot
import plotly.graph_objs as go
import sys
import os
import numpy
import re
from Bio.PDB.PDBParser import PDBParser

line_color = [
    'rgb(31,119,180)',  # muted blue
    'rgb(255,127,14)',  # safety orange
    'rgb(44,160,44)',  # cooked asparagus green
    'rgb(214,39,40)',  # brick red
    'rgb(148,103,189)',  # muted purple
    'rgb(140,86,75)',  # chestnut brown
    'rgb(227,119,194)',  # raspberry yogurt pink
    'rgb(127,127,127)',  # middle gray
    'rgb(188,189,34)',  # curry yellow-green
    'rgb(23,190,207)',   # blue-teal
]

fill_color = [
    'rgba(31,119,180,0.3)',  # muted blue
    'rgba(255,127,14,0.3)',  # safety orange
    'rgba(44,160,44,0.3)',  # cooked asparagus green
    'rgba(214,39,40,0.3)',  # brick red
    'rgba(148,103,189,0.3)',  # muted purple
    'rgba(140,86,75,0.3)',  # chestnut brown
    'rgba(227,119,194,0.3)',  # raspberry yogurt pink
    'rgba(127,127,127,0.3)',  # middle gray
    'rgba(188,189,34,0.3)',  # curry yellow-green
    'rgba(23,190,207,0.3)',   # blue-teal
]


def parse_alignment(alignment_file):
    seq = {}
    with open(alignment_file) as alnFile:
        for line in alnFile:
            m = re.match('(\S{4})\s+(\S+)', line)
            if m:
                if m.group(1) in seq:
                    seq[m.group(1)] = seq[m.group(1)] + m.group(2)
                else:
                    seq[m.group(1)] = m.group(2)
    return seq


def parse_potential(potential_file):
    df = pandas.read_csv(potential_file, skiprows=12, delim_whitespace=True, skipfooter=2,
        dtype={'resSeq': int}, engine='python',
        names=['name', 'resName', 'chainID', 'resSeq', 'potential',
            'reaction', 'coulomb', 'Ex', 'Ey', 'Ez'],
    )
    potentials =  df.groupby('resSeq')['potential'].sum()
    return potentials


def main(sim_set):
    if sim_set == 'acylphosphatase':
        active = ((23, 41), (16, 34), (21, 39), (23, 41),
                (20, 38), (18, 36), (18, 36), (24, 42))
        alignment_file = 'acylphosphatase_alignment.clustal_num'
    elif sim_set == 'dhfr':
        active = ((27, 32, 60, 80, 100, 113), (27, 49, 57, 92, 100, 121),
                (9, 64, 70), (21, 26, 31, 57, 79, 117), (12, 75),
                (5, 7, 27, 52, 57, 76, 113), (8, 63, 69), (9, 64, 70))
        alignment_file = 'dhfr_alignment.clustal_num'
    else:
        print('Error')

    alignment = parse_alignment(alignment_file)

    pdb_electrostatics_dir = f'/home/djmaity/electrostatics/{sim_set}'

    parser = PDBParser()

    data = []
    for i, prefix in enumerate(my_simulations_dict[f'{sim_set}']):
        structure = parser.get_structure(prefix,
            f'/home/djmaity/all_simulated_structures_gro/{sim_set}/{prefix}_processed.pdb')
        resSeq = [res.id[1] for res in structure[0]['A']]

        aligned_resSeq = numpy.zeros( numpy.shape(list( alignment[prefix]) ), dtype=int )
        j = 0
        for r, res_code in enumerate(alignment[prefix]):
            if res_code != '-':
                aligned_resSeq[r] = resSeq[j]
                j += 1
        directory = f'/home/djmaity/uniform_sample_electrostatics/{prefix}_uniform_sample'

        aligned_df = pandas.DataFrame(aligned_resSeq, columns=['resSeq'] )
        pdb_pot = parse_potential(f'{pdb_electrostatics_dir}/{prefix}_processed.pot')
        aligned_df['pdb'] = numpy.zeros( numpy.shape(aligned_resSeq) )
        for idx in pdb_pot.index:
            aligned_df.loc[aligned_df.resSeq == idx, 'pdb'] = pdb_pot[idx]

        pdb_trace = go.Scatter(
            name=labels[sim_set][prefix] + ' [PDB]',
            x=aligned_df['pdb'].index.values,
            y=aligned_df['pdb'],
            mode='markers',
            line=dict(color=line_color[i]),
            # text=hovertext,
            # hoverinfo='text+y',
            # legendgroup=prefix,
        )
        del aligned_df['pdb']

        column = 2
        for file in os.listdir(directory):
            if file.endswith('.pot'):
                print('Parsing: ', file)
                pot_file = os.path.join(directory, file)
                pot_series = parse_potential(pot_file)
                aligned_df[str(column)] = numpy.zeros( numpy.shape(aligned_resSeq) )
                for idx in pot_series.index:
                    aligned_df.loc[aligned_df.resSeq == idx, str(column)] = pot_series[idx]
                column += 1

        del aligned_df['resSeq']

        mean = aligned_df.mean(axis=1)
        std = aligned_df.std(axis=1)



        hover_text = []
        for n, r in zip(aligned_resSeq, alignment[prefix]):
            if r == '-':
                hover_text.append('')
            else:
                hover_text.append(f'{r} {n}')

        upper_bound = go.Scatter(
            name=labels[sim_set][prefix],
            x=mean.index.values,
            y=mean + std,
            mode='lines',
            line=dict(width=0),
            legendgroup=prefix,
            showlegend=False,
            hoverinfo='none',
            fillcolor=fill_color[i],
            fill='tonexty',
        )
        trace = go.Scatter(
            name=labels[sim_set][prefix],
            x=mean.index.values,
            y=mean,
            mode='lines',
            line=dict(color=line_color[i]),
            legendgroup=prefix,
            text=hover_text,
            hoverinfo='text+y',
            fillcolor=fill_color[i],
            fill='tonexty',
        )
        lower_bound = go.Scatter(
            name=labels[sim_set][prefix],
            x=mean.index.values,
            y=mean - std,
            mode='lines',
            line=dict(width=0),
            legendgroup=prefix,
            showlegend=False,
            hoverinfo='none',
        )
        # Trace order can be important
        # with continuous error bars
        data = data + [lower_bound, trace, upper_bound, pdb_trace]

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
        showlegend= True
    )
    fig = go.Figure(data=data, layout=layout)
    plot(fig, filename=f'{sim_set}_surface_potentials.html')


if __name__ == "__main__":
    main('acylphosphatase')
    # main('dhfr')
