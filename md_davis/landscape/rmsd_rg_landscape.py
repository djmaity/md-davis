#! /usr/bin/env python
"""
Create an HDF file for collection analysis data

Usage:
  md_davis landscape rmsd_rg [options] HDF_FILES...
  md_davis landscape rmsd_rg -h | --help

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
  --font_size <int>             Size of text and labels in the plot
  --dtick <dict>                Tick interval on each axes
  --hide_labels                 Hide the axes labels
"""

import numpy
import itertools
import docopt
import h5py

from .landscape import Landscape


def main(argv=None):
    if argv:
        args = docopt.docopt(__doc__, argv=argv)
    else:
        args = docopt.docopt(__doc__)

    shape = (int(args['--x_bins']), int(args['--y_bins']))
    start = int(args['--begin'])
    end = None if args['--end'] == None else int(args['--end'])

    if args['--temperature']:
        temperature = float(args['--temperature'])
    else:
        temperature = None

    num_files = len(args['HDF_FILES'])
    landscapes = []
    if args['--plot']:
        for filename in args['HDF_FILES']:
            landscapes += Landscape.open(filename)
        if args['--order']:
            landscape_dict = {}
            for landscape in landscapes:
                landscape_dict[landscape.name] = landscape
            landscapes = []
            for key in eval(args['--order']):
                if key in landscape_dict:
                    landscapes.append(landscape_dict[key])
                    del landscape_dict[key]
            landscapes += landscape_dict.values()
    else:
        if args['--dict']:
            landscapes = Landscape.common_landscapes(data=args['--dict'],
                shape=shape, temperature=temperature)
        else:
            input_data = {}
            for filename in args['HDF_FILES']:
                with h5py.File(filename, 'r') as hdf_file:
                    name = hdf_file.attrs['short_label']
                    label = hdf_file.attrs['short_html']
                    group = hdf_file['/rmsd_rg/all_atom']
                    time = group['time'][start:end]
                    rmsd = group['rmsd'][start:end] * 10  # Convert from nm to Å
                    rg = group['rg'][start:end] * 10    # Convert from nm to Å
                    if len(time) < 1 or len(rmsd) < 1 or len(rg) < 1:
                        raise ValueError('Invalid value for --begin or --end')
                    if args['--common']:
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

        if args['--common'] and len(input_data) > 0:
            landscapes = Landscape.common_landscapes(data=input_data,
                shape=shape, temperature=temperature)

        if args['--save']:
            for ls in landscapes:
                ls.save(filename=args['--save'],
                        name=ls.name,
                        xlabel='RMSD (in Å)',
                        ylabel='Radius of Gyration (in Å)',
                )
    if args['--hide_labels']:
        xlabel, ylabel, zlabel = '', '', ''
    else:
        xlabel=' <br>RMSD (in  Å)'
        ylabel=' <br>Rg (in  Å)'
        zlabel='Energy (kJ mol<sup>-1</sup>)<br> '

    Landscape.plot_landscapes(
        landscapes=landscapes,
        title=args['--title'],
        filename=args['--output'],
        xlabel=xlabel,
        ylabel=ylabel,
        zlabel=zlabel,
        width=int(args['--width']) if args['--width'] else None,
        height=int(args['--height']) if args['--height'] else None,
        font_size=int(args['--font_size']) if args['--font_size'] else None,
        othrographic=args['--ortho'],
        dtick=eval(args['--dtick']),
    )


if __name__ == '__main__':
    main()
