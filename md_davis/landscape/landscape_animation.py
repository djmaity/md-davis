#! /usr/bin/env python
"""
Create an HDF file for collection analysis data

Usage:
  md_davis landscape animation [options] HDF_FILE
  md_davis landscape animation -h | --help

Options:
  -h, --help                    Show this screen.
  -o, --output <filename.html>  Name for the output HTML file containing
                                the plots [default: landscapes.html]
  -i, --index int               Index of the landscape to plot [default: 0]
  -t, --title <string>          Title for the figure [default: Energy Landscape]
  -b, --begin int               Starting index for the data to include
                                [default: 0]
  -e, --end int                 Last index for the data to include
  -s, --step int                Step size for the animation data [default: 1]
  --hide_surface                Hide the energy landscape surface
  --ortho                       Orthographic projection for 3D plots
  --width <int>                 Width of the plot
  --height <int>                Height of the plot
  --font_size <int>             Size of text and labels in the plot
  --dtick <dict>                Tick interval on each axes
  --hide_labels                 Hide the axes labels
  --static                      Make a static plot with the time points colored
"""

import docopt
import plotly
import numpy
import plotly.graph_objects as go

from .landscape import Landscape


def get_animation_data(landscape):
    output = []
    for i, row in enumerate(landscape.time_bins):
        for j, time in enumerate(row):
            for t in time:
                output.append( (t, landscape.xBins[i], landscape.yBins[j], landscape.zValues[i, j]) )
    return sorted(output)


def frame_args(duration):
    return {
        "frame": {"duration": duration},
        "mode": "immediate",
        "fromcurrent": True,
        "transition": {"duration": duration, "easing": "linear"},
    }


def landscape_trajectory(landscape, data, title='Energy Landscape',
                         filename='landscape.html',
                         xlabel='x', ylabel='y', zlabel='z',
                         width=None, height=None, hide_surface=False, marker_size=8,
                         othrographic=False, font_size=None, dtick=None):
    """ Show the X and Y quiatity values from the trajectory on the landscape """
    x, y = numpy.meshgrid(landscape.xBins, landscape.yBins, indexing='ij')
    z = landscape.zValues
    time, x_val, y_val, z_val = data.T
    data = []
    if not hide_surface:
        trace = go.Surface(x=x, y=y, z=z,
            colorscale='Cividis',
            cmin=landscape.dims['z'][0],
            cmax=landscape.dims['z'][1],
            contours_z=dict(show=True, usecolormap=True, highlightcolor="limegreen", project_z=True)
        )
        data.append(trace)
    trace = go.Scatter3d(x=x_val,
                         y=y_val,
                         z=z_val,
                         mode='markers',
                         marker=dict(color=time,
                                     reversescale=True,
                                     colorbar=dict(x=0.9, title='Time'),
                                     size=marker_size,
                         ),
    )
    data.append(trace)
    fig = go.Figure(data=data)
    for axis in ['x', 'y', 'z']:
        if axis in landscape.dims:
            fig['layout'][f'scene'][axis + 'axis'].update(range=landscape.dims[axis])
    plotly.offline.plot(fig, filename=filename)


def landscape_animation(landscape, data, title='Energy Landscape',
                      filename='landscape.html',
                      xlabel='x', ylabel='y', zlabel='z',
                      width=None, height=None,
                      othrographic=False, font_size=None, dtick=None):
    x, y = numpy.meshgrid(landscape.xBins, landscape.yBins, indexing='ij')
    z = landscape.zValues
    time, x_val, y_val, z_val = data.T
    fig = go.Figure(
        data=[
            go.Scatter3d(x=x_val,
                         y=y_val,
                         z=z_val,
            ),
            go.Surface(x=x, y=y, z=z,
              colorscale='Cividis',
              cmin=landscape.dims['z'][0],
              cmax=landscape.dims['z'][1],
              contours_z=dict(show=True, usecolormap=True, highlightcolor="limegreen", project_z=True)
            )
        ],
        frames=[go.Frame(data=[go.Scatter3d(x=[x],
                                      y=[y],
                                      z=[z],
                                      mode="markers",
                                      marker=dict(color="red", size=10))],
                        name=t,
                        )
        for t, x, y, z in data]
    )
    for axis in ['x', 'y', 'z']:
        if axis in landscape.dims:
            fig['layout'][f'scene'][axis + 'axis'].update(range=landscape.dims[axis])

    sliders = [{
        "active": 0,
        "yanchor": "top",
        "xanchor": "left",
        "currentvalue": {
            "font": {"size": 20},
            "prefix": "Time:",
            "visible": True,
            "xanchor": "right"
        },
        "transition": {"duration": 300, "easing": "cubic-in-out"},
        "pad": {"b": 10, "t": 50},
        "len": 0.9,
        "x": 0.1,
        "y": 0,
        "steps": [{"args": [[f.name], frame_args(50)],
                   "label": str(k),
                   "method": "animate",
                  }
        for k, f in enumerate(fig.frames)]
    }]
    fig.update_layout(
         title=title,
         updatemenus = [
            {
                "buttons": [
                    {
                        "args": [None, frame_args(50)],
                        "label": "&#9654;", # play symbol
                        "method": "animate",
                    },
                    {
                        "args": [[None], frame_args(0)],
                        "label": "&#9724;", # pause symbol
                        "method": "animate",
                    },
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 70},
                "type": "buttons",
                "x": 0.1,
                "y": 0,
            }
         ],
         sliders=sliders
    )
    plotly.offline.plot(fig, filename=filename)


def main(argv=None):
    if argv:
        args = docopt.docopt(__doc__, argv=argv)
    else:
        args = docopt.docopt(__doc__)

    start = int(args['--begin'])
    stop = int(args['--end']) if args['--end'] is not None else None
    step = int(args['--step'])

    landscapes = Landscape.open(args['HDF_FILE'])
    landscape = landscapes[int(args['--index'])]
    data = get_animation_data(landscape)[start:stop:step]
    data = numpy.array(data)

    if args['--static']:
        landscape_trajectory(landscape=landscape, data=data, hide_surface=args['--hide_surface'])
    else:
        landscape_animation(landscape=landscape, data=data)


if __name__ == "__main__":
    main()
