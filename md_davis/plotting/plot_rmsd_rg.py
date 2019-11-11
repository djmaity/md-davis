import collections
from plotly import tools
from plotly.offline import plot
import numpy
import plotly.graph_objs as go
import h5py
import os
import argparse
import colorsys


def get_colors(num_colors, lightness=0.5, saturation=0.5, opacity=1.0):
    for hue in numpy.linspace(0, 1, num_colors, endpoint=False):
        red, green, blue = colorsys.hls_to_rgb(h=hue, l=lightness, s=saturation)
        yield f'rgba({int(red*255)}, {int(green*255)}, {int(blue*255)}, {opacity})'


def window_mean(data, window=100):
    rows = int(numpy.ceil(len(data) / window))
    means = numpy.empty((rows,))
    index = 0
    start, stop = 0, window
    while stop < len(data):
        means[index] = numpy.mean(data[start:stop])
        start = stop
        stop = stop + window
        index += 1
    means[index] = numpy.mean(data[start:])
    return means


def window_std(data, window=100):
    rows = int(numpy.ceil(len(data) / window))
    std = numpy.empty((rows,))
    index = 0
    start, stop = 0, window
    while stop < len(data):
        std[index] = numpy.std(data[start:stop])
        start = stop
        stop = stop + window
        index += 1
    std[index] = numpy.std(data[start:])
    return std


def continuous_errorbar(x, y, err, name, hover_text=None,
        line_color=None, fill_color=None, dash=None, showlegend=True):
    """ return continuous errorbar plotly trace """
    if not line_color:
        line_color = 'rgb(31,119,180)'
    if not fill_color:
        fill_color = 'rgba(31,119,180,0.3)'

    upper_bound = go.Scatter(name=name, x=x, y=y + err, mode='lines',
        line=dict(width=0), legendgroup=name, showlegend=False, hoverinfo='none',
        fillcolor=fill_color, fill='tonexty')
    trace = go.Scatter(name=name, x=x, y=y, mode='lines',
        line=dict(color=line_color, dash=dash), showlegend=showlegend,
        legendgroup=name, text=hover_text,
        hoverinfo='text+y', fillcolor=fill_color, fill='tonexty')
    lower_bound = go.Scatter(name=name, x=x, y=y - err, mode='lines',
        line=dict(width=0), legendgroup=name, showlegend=False, hoverinfo='none')
    # Trace order can be important
    # with continuous error bars
    return [lower_bound, trace, upper_bound]


def main():
    fig = tools.make_subplots(rows=1, cols=1)
    fig['layout'].update({f'yaxis2': dict(anchor='x', overlaying='y',
        side='right', showgrid=False, title='R<sub>G</sub>') })
    # TODO: Have a button to toggle between twin axes and subplots

    parser = argparse.ArgumentParser(description='Plot RMSD and Radius of gyration from HDF5 data files.')
    parser.add_argument('files', nargs='+', help='Input HDF5 data files')
    parser.add_argument('-t', '--title', default='', help='Title')
    parser.add_argument('-o', '--output', default='output.html', help='Output file name')
    parser.add_argument('-w', '--window', default=200, help='Window size for averaging')
    args = parser.parse_args()

    rmsd_traces, rg_traces = [], []

    num_files = len(args.files)
    rmsd_color = get_colors(num_files, lightness=0.5, saturation=0.8)
    rmsd_fill = get_colors(num_files, lightness=0.5, saturation=0.8, opacity=0.3)
    rg_color = get_colors(num_files, lightness=0.6, saturation=0.3)
    rg_fill = get_colors(num_files, lightness=0.6, saturation=0.3, opacity=0.3)

    for fname in args.files:
        print(fname)
        with h5py.File(fname, 'r') as datafile:
            # dset = datafile['rmsd_rg/backbone'][1:]
            dset = datafile['rmsd_rg/all_atom'][1:]
            time = window_mean(dset['time'], window=args.window) / 1000
            rmsd_traces += continuous_errorbar(
                x=time,
                y=window_mean(dset['rmsd'], window=args.window) * 10,  # Convet nm to Å
                err=window_std(dset['rmsd'], window=args.window) * 10,
                name='RMSD ' + datafile.attrs['short_html'],
                line_color=next(rmsd_color),
                fill_color=next(rmsd_fill),
            )
            rg_traces += continuous_errorbar(
                x=time,
                y=window_mean(dset['rg'], window=args.window) * 10,  # Convet nm to Å
                err=window_std(dset['rg'], window=args.window) * 10,
                name='R<sub>G</sub> ' + datafile.attrs['short_html'],
                line_color=next(rg_color),
                fill_color=next(rg_fill),
            )
    for trace in rmsd_traces:
        fig.append_trace(trace, 1, 1)
    for trace in rg_traces:
        fig.append_trace(trace, 1, 1)
        fig['data'][-1].update(yaxis=f'y2')

    fig['layout']['legend'].update(x=1.05, y=1)

    fig['layout']['xaxis1'].update(title='Time (in ns)')
    fig['layout']['yaxis1'].update(title='RMSD (in Å)')
    fig['layout'].update(title=args.title)
    plot(fig, filename=args.output)


if __name__ == "__main__":
    main()
