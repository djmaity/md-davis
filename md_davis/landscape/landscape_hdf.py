#!/usr/bin/env python3
"""
Plot free energy landscapes using RMSD and radius of gyration stored in HDF5 files.

Usage:
  command [options] HDF_FILES...
  command -h | --help

Options:
  -h, --help                    Show this screen.
  -c, --common                  Create landscapes with same ranges in
                                all three dimensions
  -s, --save <filename.h5>      Name for HDF5 file to save the landscapes
  -T, --temperature <int>       Temperature of system. If this option
                                is provided the energy landscape is
                                calculated using Boltzmann inversion,
                                else only the histogram is evaluated
  -o, --output <filename.html>  Name for the output HTML file containing
                                the plots [default: landscapes.html]
  -t, --title <string>          Title for the figure
  -x, --x_bins int              Number of bins in the X direction
                                [default: 100]
  -y, --y_bins int              Number of bins in the Y direction
                                [default: 100]
  -b, --begin int               Starting index for the data to include
                                [default: 0]
  -e, --end int                 Last index for the data to include
  -p, --plot                    Plot the precomputed energy landscape
  --ortho                       Orthographic projection for 3D plots
  -d, --dict                    Dictionary containing input data
  --order <string>              Specify the order of landscapes in the plot
  --width <int>                 Width of the plot
  --height <int>                Height of the plot
  --font_family <string>        Font style
  --font_size <int>             Size of text and labels in the plot
  --dtick <dict>                Tick interval on each axes
  --select <string>             Select: all_atom, c-alpha or backbone
                                [default: all_atom]
  --layout <string>             Layout of subplots [default: None]
"""

import numpy
import itertools
import click
import h5py
import numpy as np

from md_davis.landscape.landscape import Landscape


import click

from md_davis.landscape.landscape import Landscape
from md_davis.xvg import Xvg


def get_sorted_dset_keys(group):
    dataset_keys = []
    for key in group.keys():
        if isinstance(group[key], h5py.Dataset):
            dataset_keys.append(key)
    return sorted(dataset_keys)


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(name='landscape', context_settings=CONTEXT_SETTINGS)
# Options required for plotting .xvg files
@click.option('-x', type=click.File('r'), multiple=True,
              help='Data to plot on x-axis')
@click.option('-y', type=click.File('r'), multiple=True,
              help='Data to plot on y-axis')
@click.option('-n', '--name', type=str, multiple=True,
              help='Names of each landscape object')
@click.option('-l', '--label', type=str, multiple=True,
              help='Label to show in plots')
# Common options
@click.option('-c', '--common', is_flag=True,
              help='Use common ranges for all the landscapes')
@click.option('-T', '--temperature', type=float,
              help='Temperature of the system. If this option is provided the '
              'energy landscape is calculated using Boltzmann inversion, '
              'else only the histogram is evaluated')
@click.option('--shape', nargs=2, default=[100, 100], type=int,
              metavar="('X-bins', 'Y-bins')",
              help='Number of bins in the X and Y direction')
@click.option('-o', '--output', type=str, default='landscapes.html',
              help='Name for the output HTML file containing the plots')
@click.option('-s', '--save', metavar='FILENAME',
              help='Name for HDF5 file to save the landscapes')
@click.option('-b', '--begin', type=int, metavar='<int>', default=0,
              help='Starting index for the data to include')
@click.option('-e', '--end', type=int, metavar='<int>',
              help='Last index for the data to include')
@click.option('--limits', type=dict,
              help='A dictionary containing the limits for x, y, and x axes')
# Plotting options
@click.option('--width', type=int, metavar='<int>', help='Width of the plot')
@click.option('--height', type=int, metavar='<int>', help='Height of the plot')
@click.option('--layout', help='Layout of subplots')
@click.option('--title', type=str, metavar='<str>', help='Title for the figure')
@click.option('--axis_labels', type=dict,
              default=dict(x=' <br>RMSD (in  nm)',
                           y=' <br>Rg (in  nm)',
                           z='Free Energy (kJ mol<sup>-1</sup>)<br> '),
              help='A dictionary of strings specifying the labels for '
                   'the x, y and z-axis')
@click.option('--orthographic/--perspective', default=False,
              help='Orthographic projection for 3D plots')
@click.option('--font', type=str, metavar='<str>', help='Font style')
@click.option('--font_size', type=int, metavar='<int>',
              help='Size of text and labels in the plot')
@click.option('--dtick', help='Tick interval on each axes')
# TODO: Options
@click.option('-p', '--plot', is_flag=True, help='Plot the precomputed energy landscape')
@click.option('-d', '--dict', help='Dictionary containing input data')
@click.option('--order', type=list, help='Array specifying the order in '
              'which to plot the landscapes')
@click.option('--columns', help='Columns to use (start from second column = 1)')
# HDF option
@click.argument('hdf_files', nargs=-1)
def main(hdf_files, x, y, name, label, common, temperature, shape, output,
         save, begin, end, limits, width, height, layout, title, axis_labels,
         orthographic, font, font_size, dtick, plot, dict, order, columns):
    """Plot free energy landscapes from .xvg files generated by GROMACS
    or HDF files created by md_davis collate
    """

    # TODO: Add checks to ensure the correct number of arguments are passed

    if len(hdf_files) > 0:
        landscapes = []
        if plot:
            print(hdf_files)
            for filename in hdf_files:
                landscapes += Landscape.open(filename)
            if order:
                landscape_dict = {}
                for landscape in landscapes:
                    landscape_dict[landscape.name] = landscape
                landscapes = []
                for key in eval(order):
                    if key in landscape_dict:
                        landscapes.append(landscape_dict[key])
                        del landscape_dict[key]
                landscapes += landscape_dict.values()
        else:
            if dict:
                landscapes = Landscape.common_landscapes(data=dict,
                    shape=shape, temperature=temperature)
            else:
                input_data = {}
                for filename in hdf_files:
                    with h5py.File(filename, 'r') as hdf_file:
                        name = hdf_file.attrs['name']
                        label = hdf_file.attrs['label']
                        rmsd_rg = hdf_file['/timeseries/rmsd_rg']
                        if isinstance(rmsd_rg, h5py.Dataset):
                            time = rmsd_rg['time'][begin:end]
                            rmsd = rmsd_rg['rmsd'][begin:end]
                            rg = rmsd_rg['rg'][begin:end]
                        else:
                            time, rmsd, rg = [], [], []
                            for key in get_sorted_dset_keys(rmsd_rg):
                                time.append(rmsd_rg[key]['time'])
                                rmsd.append(rmsd_rg[key]['rmsd'])
                                rg.append(rmsd_rg[key]['rg'])
                            time = np.hstack(time)
                            rmsd = np.hstack(rmsd)
                            rg = np.hstack(rg)
                        if len(time) < 1 or len(rmsd) < 1 or len(rg) < 1:
                            raise ValueError('Invalid value for begin or end')
                        if common:
                            input_data[name] = [time, rmsd, rg, label]
                        else:
                            print(f'Generating Landscape for {name}')
                            landscape = Landscape.landscape(
                                name=name,
                                time=time,
                                x_data=rmsd,
                                y_data=rg,
                                shape=shape,
                                label=label,
                                temperature=temperature,
                            )
                            landscapes.append(landscape)
    else:
        input_data = {}
        landscapes = []
        for f1, f2, name, label in zip(x, y, name, label):
            data1 = Xvg(f1)
            data2 = Xvg(f2)

            time = data1.data[begin:end, 0]
            x_data = data1.data[begin:end, 1]
            y_data = data2.data[begin:end, 1]

            if len(time) < 1 or len(x_data) < 1 or len(y_data) < 1:
                raise ValueError('Invalid value for begin or end')

            if common:
                input_data[name] = [time, x_data, y_data, label]
            else:
                print(f'Generating Landscape for {name}')
                landscape = Landscape.landscape(
                    name=name,
                    time=time,
                    x_data=x_data,
                    y_data=y_data,
                    shape=shape,
                    label=label,
                    temperature=temperature,
                    limits=limits,
                )
                landscapes.append(landscape)

    if common and len(input_data) > 0:
        landscapes = Landscape.common_landscapes(
            data=input_data,
            shape=shape, temperature=temperature,
            limits=limits,
        )

        if save:
            for ls in landscapes:
                ls.save(filename=save,
                        name=ls.name,
                        xlabel='RMSD (in nm)',
                        ylabel='Radius of Gyration (in nm)',
                )

    Landscape.plot_landscapes(
        landscapes=landscapes,
        title=title,
        filename=output,
        axis_labels=axis_labels,
        width=width,
        height=height,
        othrographic=orthographic,
        dtick=dtick,
        layout=layout,
        font_family=font,
        font_size=font_size,
    )

if __name__ == '__main__':
    main()
