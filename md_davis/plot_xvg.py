#! /usr/bin/env python

"""
This module parses and plots the data from a .xvg file generated by
Gromacs. Although xmgrace can show .xvg plots; it is outdated, the
plots are not pretty and are difficult to customize

Usage:
  md_davis plot_xvg [options] XVG_FILE
  md_davis plot_xvg -h | --help
  md_davis plot_xvg -v | --version

Options:
  -h --help                     Show this screen
  -v --version                  Show version
  -o, --output FILENAME         Output filename
  -c, --columns <columns>       Columns to use (start from second column = 1)')
  -p, --plotly HTML_FILENAME    Use plotly for plotting. 
                                By default plots are created using matplotlib.

"""

import docopt
# import argparse
import numpy
import re


# def get_arguments():
#     """ Get filename and arguments from the commandline """
#     parser = argparse.ArgumentParser(description='Parse .xvg files created by Gromacs')
#     parser.add_argument('file', help='Input file',
#                         metavar='filename.xvg')
#     parser.add_argument('-o', '--output', dest='output',
#                         help='Output filename', metavar='output.npy')
#     parser.add_argument('-c', '--columns', dest='columns', type=int,
#                         nargs='+', metavar='int',
#                         help='Columns to use (start from second column = 1)')
#     parser.add_argument('-p', dest='plotly', help='Use plotly for plotting.'
#                         ' By default plots are created using matplotlib')
#     return parser.parse_args()


def parse_header(in_file):
    """ Parse header of .xvg files """
    with open(in_file, "r") as inp_file:
        header = []
        for current_line in inp_file:
            if current_line.startswith("@"):
                header.append(current_line[1:])
    output = {}
    title = re.compile(r'\s+title\s+\"(.+)\"')
    xlabel = re.compile(r'\s+xaxis  label\s+\"(.+)\"')
    ylabel = re.compile(r'\s+yaxis  label\s+\"(.+)\"')
    legend = re.compile(r'\s+s(\d+)\s+legend\s+\"(.+)\"')
    legend_dict = {}
    for header_line in header:
        title_match = title.search(header_line)
        xlabel_match = xlabel.search(header_line)
        ylabel_match = ylabel.search(header_line)
        legend_match = legend.search(header_line)
        if title_match:
            output['title'] = title_match.group(1)
        if xlabel_match:
            output['xlabel'] = xlabel_match.group(1)
        if ylabel_match:
            output['ylabel'] = ylabel_match.group(1)
        if legend_match:
            label = legend_match.group(2).replace('\\s', ' ').replace('\\N', '')
            legend_dict[int(legend_match.group(1))] = label
    if legend_dict:
        output['legend'] = legend_dict
    return output


def parse_xvg(input_file, columns=None, output_file=None):
    """ Parse the contents of .xvg files generated by Gromacs. """
    if columns:
        data = numpy.loadtxt(input_file, comments=('#', '@', '&'), unpack=True,
                             usecols=[0] + columns)
    else:
        data = numpy.loadtxt(input_file, comments=('#', '@', '&'), unpack=True)
        columns = range(1, len(data))
    if output_file:
        numpy.save(output_file, data)
    header = parse_header(input_file)
    if 'legend' in header:
        legend_dict = header['legend']
        labels = []
        for index in columns:
            index = index - 1  # The first coulmn is counted as 0
            if index in legend_dict:
                labels.append(legend_dict[index])
            else:
                labels.append("")
        return data, header, labels
    return data, header


def plot_xvg(data_arrays, header=None, labels=None):
    """ Plot the data parsed by parse_xvg """
    import matplotlib.pyplot as plt

    axis = plt.subplot(1, 1, 1)
    if data_arrays.ndim == 1:
        plot_array = [axis.plot(data_arrays)]
    else:
        x_val = data_arrays[0]  # First coulmn of .xvg file is time
        plot_array = []
        for y_val in data_arrays[1:]:
            current_plot, = axis.plot(x_val, y_val)
            plot_array.append(current_plot)
        if header:
            if 'title' in header:
                plt.suptitle(header['title'])
            if 'xlabel' in header:
                axis.set_xlabel(header['xlabel'])
            if 'ylabel' in header:
                axis.set_ylabel(header['ylabel'])
            if labels:
                plt.legend(plot_array, labels)
    plt.show()


def plotly_xvg(filename, data_arrays, header=None, labels=None):
    """ Plot the data parsed by parse_xvg using plotly"""
    from plotly.offline import plot
    import plotly.graph_objs as go

    layout = dict(titlefont=dict(size=30),
                  xaxis = dict(titlefont=dict(size=24),
                               tickfont=dict(size=20),
                  ),
                  yaxis = dict(titlefont=dict(size=24),
                               tickfont=dict(size=20),
                               tickformat='.2f',
                  ),
    )
    if data_arrays.ndim == 1:
        trace = go.Scatter(
            x=list(range(len(data_arrays))),
            y=data_arrays,
            mode = 'lines+markers',
        )
        data = [trace]
    else:
        x_val = data_arrays[0]  # First coulmn of .xvg file is time
        data = []
        for y_val in data_arrays[1:]:
            trace = go.Scatter(
                x=x_val,
                y=y_val,
                mode = 'lines+markers',
                name='',
            )
            data.append(trace)

        if header:
            if 'title' in header:
                layout['title'] = header['title']
            if 'xlabel' in header:
                layout['xaxis']['title'] = header['xlabel']
            if 'ylabel' in header:
                layout['yaxis']['title'] = header['ylabel']
            if labels:
                for i, label in enumerate(labels):
                    data[i].name = label
    figure = dict(data=data, layout=layout)
    plot(figure, filename=filename)


def main(args):
    """ Plot .xvg files genrated by Gromacs """
    if args['--columns']:
        columns = args['--columns'].strip().split()
    else:
        columns = None
    parsed_data = parse_xvg(input_file=args['XVG_FILE'],
                            columns=columns,
                            output_file=args['--output'],
    )
    if args['--plotly']:
        plotly_xvg(args['--plotly'], *parsed_data)
    else:
        plot_xvg(*parsed_data)


if __name__ == '__main__':
    arguments = docopt.docopt(__doc__)
    main(arguments)
