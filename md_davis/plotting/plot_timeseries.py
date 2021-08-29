#!/usr/bin/env python3

from warnings import simplefilter
import numpy
import plotly
from plotly.offline import plot
import plotly.graph_objs as go
import h5py
import click
import colorsys

from md_davis.common.stats import *


class TimeSeries:
    """Create a time series object from data like RMSD and radius of gyration """

    def __init__(self, time, data, unit=None, label=None) -> None:
        self.time = time
        self.data = data
        self.unit = unit
        self.label = label
        super().__init__()


    def from_xvg(filename):
        """Load data from .xvg file """
        pass


    def from_hdf5(filename):
        """Load data from HDF5 file """
        pass



def get_colors(num_colors, lightness=0.5, saturation=0.5, opacity=1.0):
    for hue in numpy.linspace(0, 1, num_colors, endpoint=False):
        red, green, blue = colorsys.hls_to_rgb(h=hue, l=lightness, s=saturation)
        yield f'rgba({int(red*255)}, {int(green*255)}, {int(blue*255)}, {opacity})'


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


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(name='timeseries', context_settings=CONTEXT_SETTINGS)
@click.option('-o', '--output', default='output.html', help='Output file name')
@click.option('-t', '--title', default='', help='Title')
@click.option('-w', '--window', default=100, help='Window size for averaging')
@click.argument('files')
def timeseries_plot(args):
    fig = plotly.subplots.make_subplots(rows=1, cols=1)
    fig['layout'].update({f'yaxis2': dict(anchor='x', overlaying='y',
        side='right', showgrid=False, title='R<sub>G</sub>') })
    # TODO: Have a button to toggle between twin axes and subplots
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
            dset = datafile['rmsd_rg/all_atom'][1:]  # TODO: remove these hard coded locations

            time = moving_avg(dset['time'], window=args.window, shift=args.window)
            rmsd, rmsd_err = moving_avg_std(dset['rmsd'], window=args.window)
            rg, rg_err = moving_avg_std(dset['rg'], window=args.window)

            rmsd_traces += continuous_errorbar(
                x=time,
                y=rmsd,
                err=rmsd_err,
                name='RMSD ' + datafile.attrs['short_html'],
                line_color=next(rmsd_color),
                fill_color=next(rmsd_fill),
            )
            rg_traces += continuous_errorbar(
                x=time,
                y=window_mean(dset['rg'], window=args.window),
                err=window_std(dset['rg'], window=args.window),
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
    timeseries_plot()
