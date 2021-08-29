#! /usr/bin/env python
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

from md_davis.landscape.landscape import Landscape


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(name='landscape', context_settings=CONTEXT_SETTINGS)
@click.option('--layout', nargs=2, help='Layout of subplots')
@click.option('--select', default='all_atom',
              type=click.Choice(['all_atom', 'backbone', 'c-alpha'],
                                 case_sensitive=False),
              help='Select which data to plot: all_atom, c-alpha or backbone')
# @click.option('--axis_labels',
#               default=dict(x=' <br>RMSD (in  Å)',
#                            y=' <br>Rg (in  Å)',
#                            z='Free Energy (kJ mol<sup>-1</sup>)<br> '),
#               help='A dictionary of strings specifying the labels for '
#                    'the x, y and z-axis')
# @click.option('-l', '--labels', nargs=-1, type=str,
#               help='Label to show in plots')
# @click.option('-n', '--name', action='append', type=str,
#               help='Name for the landscape object plots')
# @click.option('-c', '--common', action='store_true',
#               help='Use the same ranges for all landscapes')
# @click.option('-T', '--temperature', type=float, metavar='<float>',
#               help='Temperature of the system. If this option is provided '
#                    'the energy landscape is calculated using Boltzmann '
#                    'inversion, else only the histogram is evaluated')
# @click.option('--shape', nargs=2, default=[100, 100], type=int,
#               help='Number of bins in the X and Y direction')
# @click.option('-s', '--save',
#               help='Name for HDF5 file to save the landscapes')
# @click.option('-o', '--output', type=str, default='landscapes.html',
#               help='Name for the output HTML file containing the plots')
# @click.option('-t', '--title', type=str, help='Title for the figure')
# @click.option('-b', '--begin', type=int, default=0,
#               help='Starting index for the data to include')
# @click.option('-e', '--end', type=int, help='Last index for the data to include')
# @click.option('-p', '--plot', help='Plot the precomputed energy landscape')
# @click.option('--ortho', is_flag=True,
#               help='Orthographic projection for 3D plots')
# @click.option('-d', '--dict', help='Dictionary containing input data')
# @click.option('--width', type=int, help='Width of the plot')
# @click.option('--height', type=int, help='Height of the plot')
# @click.option('--font', type=str, help='Font style')
# @click.option('--font_size',
#     help='Size of text and labels in the plot')
# @click.option('--dtick', help='Tick interval on each axes')
@click.argument('hdf_files', nargs=-1)
def main(hdf_files=None,
         layout=None,
         select=None,
         axis_labels=None,
         label=None,
         name=None,
         common=None,
         temperature=None,
         shape=None,
         dict=None,
         save=None,
         output=None,
         title=None,
         plot=None,
         ortho=False,
         width=None,
         height=None,
         start=0,
         end=None,
         font=None,
         font_size=None,
         dtick=None
    ):

    print(hdf_files,
         layout,
         select,
         axis_labels,
         label,
         name,
         common,
         temperature,
         shape,
         dict,
         save,
         output,
         title,
         plot,
         ortho,
         width,
         height,
         start,
         end,
         font,
         font_size,
         dtick)
    # shape = (int(x_bins), int(y_bins))
    # if temperature:
    #     temperature = float(temperature)
    # else:
    #     temperature = None
    #
    # num_files = len(HDF_FILES)
    # landscapes = []
    # if plot:
    #     for filename in HDF_FILES:
    #         landscapes += Landscape.open(filename)
    #     if order:
    #         landscape_dict = {}
    #         for landscape in landscapes:
    #             landscape_dict[landscape.name] = landscape
    #         landscapes = []
    #         for key in eval(order):
    #             if key in landscape_dict:
    #                 landscapes.append(landscape_dict[key])
    #                 del landscape_dict[key]
    #         landscapes += landscape_dict.values()
    # else:
    #     if dict:
    #         landscapes = Landscape.common_landscapes(data=dict,
    #             shape=shape, temperature=temperature)
    #     else:
    #         input_data = {}
    #         for filename in HDF_FILES:
    #             with h5py.File(filename, 'r') as hdf_file:
    #                 name = hdf_file.attrs['short_label
    #                 label = hdf_file.attrs['short_html
    #                 group = hdf_file['/rmsd_rg/' + select]
    #                 time = group['time[start:end]
    #                 rmsd = group['rmsd[start:end] * 10  # Convert from nm to Å
    #                 rg = group['rg[start:end] * 10    # Convert from nm to Å
    #                 if len(time) < 1 or len(rmsd) < 1 or len(rg) < 1:
    #                     raise ValueError('Invalid value for begin or end')
    #                 if common:
    #                     input_data[name] = [time, rmsd, rg, label]
    #                 else:
    #                     print(f'Generating Landscape for {name}')
    #                     landscape = Landscape.landscape(
    #                         name=name,
    #                         time=time,
    #                         x_data=rmsd,
    #                         y_data=rg,
    #                         shape=shape,
    #                         label=label,
    #                         temperature=temperature,
    #                     )
    #                     landscapes.append(landscape)
    #
    #     if common and len(input_data) > 0:
    #         landscapes = Landscape.common_landscapes(data=input_data,
    #             shape=shape, temperature=temperature,
    #             # dimensions={'x': (0.5, 5.5), 'y': (23.2, 25), 'z': (-20, 0)}  # Hardcoded values remember to remove
    #             )
    #
    #     if save:
    #         for ls in landscapes:
    #             ls.save(filename=save,
    #                     name=ls.name,
    #                     xlabel='RMSD (in Å)',
    #                     ylabel='Radius of Gyration (in Å)',
    #             )
    #
    # Landscape.plot_landscapes(
    #     landscapes=landscapes,
    #     title=title,
    #     filename=output,
    #     axis_labels=eval(axis_labels),
    #     width=int(width) if width else None,
    #     height=int(height) if height else None,
    #     font_size=int(font_size) if font_size else None,
    #     othrographic=ortho,
    #     dtick=eval(dtick) if dtick else None,
    #     layout=eval(layout) if layout else None,
    #     font_family=font_family)


if __name__ == '__main__':
    main()
