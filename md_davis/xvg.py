# -*- coding: utf-8 -*-
"""
This module parses and plots the data from xmgrace (.xvg) files
generated by GROMACS. Although xmgrace can show such plots; the look is
outdated and difficult for quick comparison of multiple files.
"""

import numpy
import os
import re

import click
import matplotlib.pyplot as plt
from plotly.offline import plot
import plotly.graph_objs as go


class Xvg:
    def __init__(self, xvg_file) -> None:
        """ Parse header of .xvg files and create an Xvg object """

        title = re.compile(r'\s+title\s+\"(.+)\"')
        subtitle = re.compile(r'\s+subtitle\s+\"(.+)\"')
        xlabel = re.compile(r'\s+xaxis  label\s+\"(.+)\"')
        ylabel = re.compile(r'\s+yaxis  label\s+\"(.+)\"')
        legend = re.compile(r'\s+s(\d+)\s+legend\s+\"(.+)\"')

        data = []
        legend_dict = {}
        for current_line in xvg_file:

            # Parse header lines
            if current_line.startswith("@"):
                title_match = title.search(current_line)
                subtitle_match = subtitle.search(current_line)
                xlabel_match = xlabel.search(current_line)
                ylabel_match = ylabel.search(current_line)
                legend_match = legend.search(current_line)
                if title_match:
                    self.title = title_match.group(1)
                if subtitle_match:
                    self.subtitle = subtitle_match.group(1)
                if xlabel_match:
                    self.xlabel = xlabel_match.group(1)
                    legend_dict[0] = self.xlabel
                if ylabel_match:
                    self.ylabel = ylabel_match.group(1)
                if legend_match:
                    label = legend_match.group(2).replace('\\s', ' ').replace('\\N', '')
                    legend_dict[int(legend_match.group(1)) + 1] = label

            # Ignore comments
            elif current_line.startswith("#") or current_line.startswith("&"):
                continue

            # Parse data
            else:
                row = current_line.strip().split()
                data.append([float(_) for _ in row])

        self.mode = 'lines'
        self.name = os.path.splitext(os.path.basename(xvg_file.name))[0]
        self.data = numpy.array(data)
        self.legend = legend_dict

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return (
            '\ntitle = ' + self.title + '\n' +
            'subtitle = ' + self.subtitle + '\n' +
            'xlabel = ' + self.xlabel + '\n' +
            'ylabel = ' + self.ylabel + '\n' +
            'legend = ' + self.legend.__repr__() + '\n' +
            'data = ' + self.data.__repr__() + '\n'
        )


def plot_xvg(data, filename=None, xcol=0):
    """Plot multiple Xvg objects using matplotlib."""
    axis = plt.subplot(1, 1, 1)
    for xvg_data in data:
        if xvg_data.data.ndim == 1:
            if xvg_data.legend:
                label = xvg_data.name + ' ' + xvg_data.legend[1]
            else:
                label = xvg_data.name
            axis.plot(xvg_data, label=label)
        else:
            x_val = xvg_data.data[:, xcol]
            _, cols = numpy.shape(xvg_data.data)
            for index in range(1, cols):
                y_val = xvg_data.data[:, index]
                if len(xvg_data.legend) > 1:
                    label = xvg_data.name + ' ' + xvg_data.legend[index]
                else:
                    label = xvg_data.name
                axis.plot(x_val, y_val, label=label)
        if xvg_data.title:
            plt.suptitle(xvg_data.title)
        if xvg_data.xlabel:
            axis.set_xlabel(xvg_data.xlabel)
        if xvg_data.ylabel:
            axis.set_ylabel(xvg_data.ylabel)
    plt.legend()
    if filename is not None:
        plt.savefig(fname=filename)
    else:
        plt.show()


def plotly_xvg(data, filename=None, layout=None):
    """Plot multiple Xvg objects using plotly.

    Arguments:
        data (array): Input array of Xvg objects to plot
        filename (string): output filename or path ending in .html
        layout (dict): plotly layout dict

    """

    if filename is None:
        filename = 'plot_xvg.html'

    if layout is None:
        layout = dict(titlefont=dict(size=18),
                      xaxis=dict(titlefont=dict(size=14),
                                 tickfont=dict(size=12),
                                 ),
                      yaxis=dict(titlefont=dict(size=14),
                                 tickfont=dict(size=12),
                                 #  tickformat='.2f',
                                 ),
                      )
    else:
        layout = eval(layout)

    trace_data = []
    for xvg_data in data:
        if xvg_data.data.ndim == 1:
            if xvg_data.legend:
                label = xvg_data.name + ' ' + xvg_data.legend[1]
            else:
                label = xvg_data.name
            trace = go.Scatter(
                x=list(range(len(xvg_data.data))),
                y=xvg_data,
                mode=xvg_data.mode,
                name=label
            )
            trace_data.append(trace)
        else:
            _, cols = numpy.shape(xvg_data.data)
            for index in range(1, cols):
                if len(xvg_data.legend) > 1:
                    label = xvg_data.name + ' ' + xvg_data.legend[index]
                else:
                    label = xvg_data.name
                trace = go.Scatter(
                    x=xvg_data.data[:, 0],
                    y=xvg_data.data[:, index],
                    mode=xvg_data.mode,
                    name=label,
                )
                trace_data.append(trace)

        if xvg_data.title:
            layout['title'] = xvg_data.title
        if xvg_data.xlabel:
            layout['xaxis']['title'] = xvg_data.xlabel
        if xvg_data.ylabel:
            layout['yaxis']['title'] = xvg_data.ylabel
        # if xvg_data.legend:
        #     for
        #     layout['yaxis']['title'] = xvg_data.legend

    figure = dict(data=trace_data, layout=layout)
    plot(figure, filename=filename)


# TODO: Allow columns to be selected by fixing '-c' option
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(name='xvg', context_settings=CONTEXT_SETTINGS)
@click.option('-o', '--output', default=None, help='Output filename')
@click.option('--matplotlib/--plotly', 'use_matplotlib',
              help='Use plotly for plotting')
# @click.option('-c', '--columns',
#               help='Choose columns to plot')
@click.argument('xvg_files', nargs=-1, type=click.File('r'))
def main(xvg_files, output, columns=None, use_matplotlib=False):
    """Plot xmgrace (.xvg) files generated by GROMACS.

    Arguments:\n
        xvg_files (Files): Input xmgrace (.xvg) files'
    """

    data = []
    for xvg_file in xvg_files:
        xvg_data = Xvg(xvg_file)
        data.append(xvg_data)

    if use_matplotlib:
        plot_xvg(data=data, filename=output)
    else:
        plotly_xvg(data=data, filename=output)


if __name__ == '__main__':
    main()
