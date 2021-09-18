import collections
from plotly import tools
from plotly.offline import plot
import numpy as np
import pickle
import plotly
import plotly.graph_objs as go
import sys
import string
import h5py
import os
import pandas

structure_info = collections.namedtuple('structure_info', ['label', 'html', 'color'])

secondary_structure = collections.OrderedDict([
    ('H', structure_info(label='α-helix', html='α helix', color='rgb(255, 0, 0)') ),
    ('G', structure_info(label='3_10-helix', html='3<sub>10</sub> helix', color='rgb(120, 0, 100)') ),
    ('I', structure_info(label='π-helix', html='π helix', color='rgb(255, 15, 139)') ),
    ('E', structure_info(label='β strand', html='β strand', color='rgb(0, 0, 255)') ),
    ('B', structure_info(label='β-bridge', html='β bridge', color='rgb(15, 125, 255)') ),
    ('T', structure_info(label='Turn', html='Turn', color='rgb(0, 200, 0)') ),
    ('S', structure_info(label='Bend', html='Bend', color='rgb(120, 200, 0)') ),
    ('~', structure_info(label='Loop', html='Loop', color='rgb(159, 255, 15)') ),
])

# def plot_secondary_structure(data, name, figure, twin_axis, row=1, column=1, height=100):
def plot_secondary_structure(data, name, figure, axis=None, row=1, column=1, height=100, showlegend=True):
    """ Plot secondary structure stacked bar graph """
    if isinstance(data, pandas.DataFrame):
        frames = data.sum(axis=1).max()
    elif isinstance(data, h5py.Dataset):
        frames = sum(data[0])
    else:
        raise TypeError

    for code, structure in secondary_structure.items():
        if isinstance(data, pandas.DataFrame):
            x_values = data.index + 1
        if isinstance(data, h5py.Dataset):
            x_values = list(range(1, data.len() + 1 )) + 1
        y_values = data[code] * (height / frames)
        trace = go.Bar(
            x=x_values,
            y=y_values,
            name=name,
            legendgroup=name,
            showlegend=False,
            marker={'color': structure.color,
                    'opacity': 0.3,
            },
            hoverinfo='none',
        )
        figure.append_trace(trace, row=row, col=column)
        if axis:
            figure['data'][-1].update(yaxis=axis)
    if showlegend:
        figure['data'][-1]['showlegend']=True
    figure['layout']['barmode'] = 'stack'


def add_secondary_structure_legend(figure, xloc=0.8, yloc=-0.2, xspace=0, yspace=-0.06):
    dssp_legend = []
    for index, structure in enumerate(secondary_structure.values()):
        dssp_legend.append(
            dict(x=xloc + index * xspace,
                 y=yloc + index * yspace,
                 showarrow=False,
                 text=structure.html,
                 xref='paper',
                 yref='paper',
                 font=dict(size=24, color=structure.color),
                 xanchor='left', ax=0, ay=0,
            )
        )
    return dssp_legend


def plot_dihedral_sd(dataset, name, figure, twin_axis, row=1, column=1, height=100):
    """ Add dihedral standard deviation into a figure """
    if twin_axis != '':
        figure['layout'].update({f'yaxis{twin_axis}': dict(
            anchor='x', overlaying='y', side='left', showgrid=False) })
    for index, structure in enumerate(SECSTR_CODES):
        trace = go.Bar(
            x=list(range(dataset.len())),
            y=dataset[structure] / sum(dataset[0]) * height,
            name=name,
            legendgroup=name,
            showlegend=False,
            marker={'color': plotly.colors.DEFAULT_PLOTLY_COLORS[index],
                    'opacity': 0.6
            },
            hoverinfo='none',
        )
        figure.append_trace(trace, row=row, col=column)
        figure['data'][-1].update(yaxis=f'y{twin_axis}')
    figure['data'][-1]['showlegend']=True
    figure['layout']['barmode'] = 'stack'


def main():
    directory = '/home/djmaity/Lab_Home/my_research/research_sorted/acylphosphatase/acylphosphatase_hdf5_data/'
    fig = tools.make_subplots(rows=1, cols=1)
    count = 1
    for file in os.listdir(directory):
        if file == '2VH7_data_bak.h5':
            continue
        print(file)
        with h5py.File(directory + file, 'r') as datafile:
            for chain, dataset in datafile['secondary_structure/counts_per_residue'].items():
                plot_secondary_structure(data=dataset,
                                         name=file[:-8] + ' Secondary Structure',
                                         figure=fig,
                                         twin_axis='' if count <= 1 else count
                )
        count += 1

    fig['layout']['font'].update(family='Courier New, monospace', size=24, color='black')

    annotations = add_secondary_structure_legend(figure=fig)
    annotations += [
        dict(x=0, y=0.5, showarrow=True, text='Root Mean Square Fluctutation (in Å)',
             arrowcolor='rgba(0,0,0,0)',
             textangle=-90, xref='paper', yref='paper', font=dict(size=24), ax=-60, ay=0),
        dict(x=1, y=0.5, showarrow=True, text='Dihedral Standard Deviation (in degrees)',
             arrowcolor='rgba(0,0,0,0)',
             textangle=-90, xref='paper', yref='paper', font=dict(size=24), ax=50, ay=0),
        dict(x=0.5, y=0, showarrow=True, text='Residue Index',
             arrowcolor='rgba(0,0,0,0)',
             xref='paper', yref='paper', font=dict(size=24), ax=0, ay=50),
    ]
    fig['layout']['annotations'] = annotations

    fig['layout']['legend'].update(x=0, y=-0.1, xanchor='left', yanchor='top')
    fig['layout'].update(title='Comparison of HbA and HbS in MD simulation of crystal',
        autosize=False, width=1200, height=1200,
        margin=go.layout.Margin(l=80, r=100, b=450, t=80, pad=2)
    )
    plot(fig, filename='test.html')


if __name__ == "__main__":
    main()
