import plotly
import numpy
import plotly.graph_objects as go

from md_davis.landscape.landscape import Landscape

def frame_args(duration):
    return {
            "frame": {"duration": duration},
            "mode": "immediate",
            "fromcurrent": True,
            "transition": {"duration": duration, "easing": "linear"},
        }


def get_animation_data(landscape):
    output = []
    for i, row in enumerate(landscape.time_bins):
        for j, time in enumerate(row):
            for t in time:
                output.append( (t, landscape.xBins[i], landscape.yBins[j], landscape.zValues[i, j]) )
    return sorted(output)


def plot_landscapes(landscape, title='Energy Landscape', filename='landscape.html',
                    xlabel='x', ylabel='y', zlabel='z'):
    x, y = numpy.meshgrid(landscape.xBins, landscape.yBins, indexing='ij')
    z = landscape.zValues
    animation_data = get_animation_data(landscape)[::10]
    dat = list(zip(*animation_data))
    fig = go.Figure(
        data=[
            go.Scatter3d(x=dat[1],
                         y=dat[2],
                         z=dat[3],
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
                                      marker=dict(color="red", size=10))])
        for _, x, y, z in animation_data]
    )
    for axis in ['x', 'y', 'z']:
        if axis in landscape.dims:
            fig['layout'][f'scene'][axis + 'axis'].update(range=landscape.dims[axis])

    fig.layout.updatemenus = [
        {
            "buttons": [
                {
                    "args": [None],
                    "label": "Play",
                    "method": "animate"
                },
                {
                    "args": [[None], {"frame": {"duration": 0, "redraw": False},
                                    "mode": "immediate",
                                    "transition": {"duration": 0}}],
                    "label": "Pause",
                    "method": "animate"
                }
            ],
            "direction": "left",
            "pad": {"r": 10, "t": 87},
            "showactive": False,
            "type": "buttons",
            "x": 0.1,
            "xanchor": "right",
            "y": 0,
            "yanchor": "top"
        }
    ]

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
        "steps": [{"args": [[f.name], frame_args(0)],
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
    fig.show()

def main():
    landscapes = Landscape.open('/home/djmaity/Lab_Home/my_research/research_sorted/hemoglobin/15_hemoglobin_crystal_hdf5/Hb_Crystal-landscape.h5')
    plot_landscapes(landscapes[0])

if __name__ == "__main__":
    main()
