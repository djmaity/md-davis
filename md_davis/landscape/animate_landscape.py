import plotly
import numpy
import plotly.graph_objects as go

def get_animation_data(landscape):
    output = []
    for i, row in enumerate(landscape.time_bins):
        for j, time in enumerate(row):
            for t in time:
                output.append( (t, landscape.xBins[i], landscape.yBins[j], landscape.zValues[i, j]) )
    return sorted(output)

def plot_landscapes(landscape, title='Energy Landscape', filename='landscape.html',
                    xlabel='x', ylabel='y', zlabel='z'):
    x, y = numpy.meshgrid(landscape.xBins, landscape.yBins)
    z = landscape.zValues
    animation_data = get_animation_data(landscape)[:100]
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
    fig.update_layout(title=title)

    # for axis in ['x', 'y', 'z']:
    #     if axis in landscape.dims:
    #         fig['layout'][f'scene'][axis + 'axis'].update(range=landscape.dims[axis])

    # fig.layout.updatemenus=[dict(type="buttons",
    #                              buttons=[dict(label="Play",
    #                              method="animate",
    #                              args=[None])])]

    fig.layout.updatemenus = [
        {
            "buttons": [
                {
                    "args": [None, {"frame": {"duration": 500, "redraw": False},
                                    "fromcurrent": True, "transition": {"duration": 300,
                                                                        "easing": "quadratic-in-out"}}],
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

    sliders_dict = {
        "active": 0,
        "yanchor": "top",
        "xanchor": "left",
        "currentvalue": {
            "font": {"size": 20},
            "prefix": "Year:",
            "visible": True,
            "xanchor": "right"
        },
        "transition": {"duration": 300, "easing": "cubic-in-out"},
        "pad": {"b": 10, "t": 50},
        "len": 0.9,
        "x": 0.1,
        "y": 0,
        "steps": []
    }
    fig.show()

def main():
    pass

if __name__ == "__main__":
    main()
