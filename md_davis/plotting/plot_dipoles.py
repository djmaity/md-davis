#! /usr/env/python

"""
Plot dipole moment from hdf5 files

Usage:
  md_davis plot dipoles [options] [(--3d | --both)] FILES...

Options:
  -t, --title <string>      Title for the plot
  -o, --output FILENAME     Output filename
  -s, --size <int>          Number of points to plot
                            [default: 1000]
  --3d                      Make 3D plot
  -b, --both                Make both a 3D plot and a 2D polar plot
  -c, --center              Show the center in the 3D plot
"""

import numpy
import docopt
import h5py
import plotly.graph_objs as go
from plotly.offline import plot
import plotly
from plotly.colors import DEFAULT_PLOTLY_COLORS

# Local imports
from md_davis.common import polar


def sample(array, size=None):
    if not isinstance(array, (list, tuple, numpy.ndarray) ):
        raise ValueError
    if size:
        step = int(len(array) / size)
    else:
        step = 1
    return array[::step]


def plot_dipoles(dictionary, filename=None, title=None):
    fig = plotly.subplots.make_subplots(rows=3, cols=1, vertical_spacing=0.03, shared_xaxes=True)
    for i, (label, value) in enumerate(dictionary.items()):
        magnitude, azimuth, inclination = value
        x_data = list(range(len(magnitude)))

        chosen_color = DEFAULT_PLOTLY_COLORS[i % len(DEFAULT_PLOTLY_COLORS)]
        trace = go.Scatter(
            x=x_data,
            y=magnitude,
            legendgroup=label,
            name=label,
            mode='lines',
            marker=dict(
                color=chosen_color,
            ),
        )
        fig.append_trace(trace, 1, 1)

        trace = go.Scatter(
            x=x_data,
            y=numpy.degrees(azimuth),
            legendgroup=label,
            name=label,
            mode='markers',
            marker=dict(
                size=3,
                color=chosen_color,
            ),
            showlegend=False,
        )
        fig.append_trace(trace, 2, 1)

        trace = go.Scatter(
            x=x_data,
            y=numpy.degrees(inclination),
            legendgroup=label,
            name=label,
            mode='lines',
            marker=dict(
                color=chosen_color,
            ),
            showlegend=False,
        )
        fig.append_trace(trace, 3, 1)

    if title:
        fig['layout'].update(title=title)

    fig['layout']['xaxis3'].update(range=[0, 1000], title='Time (ns)')
    fig['layout']['yaxis1'].update(title='Dipole Moment Magnitude (Debye)')
    fig['layout']['yaxis2'].update(title='Azimuth (degrees)', range=[-181, 181],
                                   tick0=-180, dtick=60)
    fig['layout']['yaxis3'].update(title='Inclination (degrees)', range=[0, 181],
                                   tick0=0, dtick=30)
    if filename:
        plot(fig, filename=filename)
    else:
        fig.show()


def plot_dipoles3d(dictionary, filename=None, title=None, centroid=False):
    data = []
    for i, (label, value) in enumerate(dictionary.items()):
        x_data, y_data, z_data = value
        chosen_color = DEFAULT_PLOTLY_COLORS[i % len(DEFAULT_PLOTLY_COLORS)]
        trace = go.Scatter3d(
            x=x_data,
            y=y_data,
            z=z_data,
            name=label,
            legendgroup=label,
            marker=dict(
                color=chosen_color,
                size=2,
                opacity=0.2,
            ),
            mode='markers'
        )
        data.append(trace)
        if centroid:
            center = numpy.mean(value, axis=1)
            centroid_trace = go.Scatter3d(
                x=[center[0]],
                y=[center[1]],
                z=[center[2]],
                name=label,
                legendgroup=label,
                marker=dict(
                    color=chosen_color,
                ),
                mode='markers'
            )
            data.append(centroid_trace)

    fig = go.Figure(data=data)
    if title:
        fig['layout'].update(title=title)
    fig.update_layout(scene = dict(
                        xaxis = dict(
                            gridcolor="white",
                            showbackground=True,
                            zerolinecolor="black",),
                        yaxis = dict(
                            gridcolor="white",
                            showbackground=True,
                            zerolinecolor="black"),
                        zaxis = dict(
                            gridcolor="white",
                            showbackground=True,
                            zerolinecolor="black",),),
    )

    if filename:
        plot(fig, filename=filename)
    else:
        fig.show()



def main(argv):
    """ Plot dipole moment from HDF5 files """
    if argv:
        args = docopt.docopt(__doc__, argv=argv)
    else:
        args = docopt.docopt(__doc__)
    size = int(args['--size'])
    output = args['--output']
    title = args['--title']
    polar_dipoles = {}
    dipoles_dict = {}
    for fname in args['FILES']:
        hdf_file = h5py.File(fname, 'r')
        data = numpy.asarray(hdf_file['dipole_moment'])
        dipoles = numpy.vstack([sample(data['mu_x'], size=size),
                                sample(data['mu_y'], size=size),
                                sample(data['mu_z'], size=size),
                                ])
        if args['--3d']:
            dipoles_dict[ hdf_file.attrs['short_html'] ] = dipoles
        else:
            polar_dipoles[ hdf_file.attrs['short_html'] ] = polar.spherical_numpy(xyz=dipoles.T).T

        if args['--both']:
            dipoles_dict[ hdf_file.attrs['short_html'] ] = dipoles

    if args['--3d']:
        plot_dipoles3d(dipoles_dict, filename=output, title=title, centroid=args['--center'])
    else:
        plot_dipoles(polar_dipoles, filename=output, title=title)

    if args['--both']:
        plot_dipoles3d(dipoles_dict, filename='3D_' + output, title=title, centroid=args['--center'])


if __name__ == '__main__':
    main()
