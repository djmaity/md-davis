"""Free energy landscape from two properties

Description

    Typical usage example:

    fel = Landscape('my_landscapes',
                              [1.0, 2.0, 3.0, 4.0, 5.0],
                              [10.5, 12.5, 14.5, 16.5, 18.5])

"""

import math
import numpy
import h5py
import itertools

import plotly.offline as py
import plotly.subplots

# Exact value after redifinition of SI units in May 2019
R = 0.00831446261815324  # Energy in KJ / mol


class Landscape:
    """A class for Free energy landscapes (FEL).

    Create a FEL object and populate it with data from two xvg files or
    data from MD DaVis HDF5 file.

    Parameters
    ----------
        name : string without spaces
            A name for the object; crucial when saving the FEL in HDF5 format.
        xbins : array_like
            The bins in the x-dimension
        ybins : array_like
            The bins in the y-dimension
        label : string, optional
            Label to be displayed in plots. Same as name if not provided.

    Attributes
    ----------
        binned_data : list of lists of lists
            The binned data, each bin of binned_data[i][j] stores the time
            of each frame in the bin corresponding to xbins[i] and ybins[j]
        zValues : matrix of floats
            The two-dimensional histogram of the data or free energy of
            each bin.
        limits : dict
            A dict with the keys 'x', 'y' and 'z' each specifying the
            ranges along x, y, and z dimension
    """

    def __init__(self, name, xbins, ybins, label=None):
        """Constructor for Landscape"""
        self.name = name
        self.label = label if label else name
        self.xbins = numpy.array(xbins)
        self.ybins = numpy.array(ybins)
        self.binned_data = [[[] for _ in range(len(ybins))] for _ in range(len(xbins))]
        self.zValues = numpy.zeros((len(xbins), len(ybins)))
        self.limits = {'x': None, 'y': None, 'z': None}

    def __repr__(self):
        return 'Name:  ' + self.name + '\n' \
                         + 'Label: ' + self.label + '\n' \
                         + self.limits.__repr__() + '\n\n' \
                         + 'Z values:\n' + self.zValues.__repr__() + '\n\n' \
                         + 'X bins:\n' + self.xbins.__repr__() + '\n\n' \
                         + 'Y bins:\n' + self.ybins.__repr__() + '\n\n'

    def __str__(self):
        return 'Name:  ' + self.name + '\n' \
                         + 'Label: ' + self.label + '\n' \
                         + self.limits.__str__() + '\n\n' \
                         + 'Z values:\n' + self.zValues.__str__() + '\n\n' \
                         + 'X bins:\n' + self.xbins.__str__() + '\n\n' \
                         + 'Y bins:\n' + self.ybins.__str__() + '\n\n'

    def __iter__(self):
        for i, rows in enumerate(self.binned_data):
            for j, time in enumerate(rows):
                yield i, j, time

    def __len__(self):
        return len(self.xbins)

    def add_data(self, time, x_data, y_data):
        """ Create a 2D histogram and add times to the corresponding bins """
        dx = numpy.digitize(x_data, self.xbins)
        dy = numpy.digitize(y_data, self.ybins)
        for t in range(len(time)):
            self.zValues[dx[t] - 1][dy[t] - 1] += 1
            self.binned_data[dx[t] - 1][dy[t] - 1].append(time[t])
        return self

    def save(self, filename, name, xlabel='', ylabel=''):
        print(f'Saving {name} ...')
        with h5py.File(filename, 'a') as hdf_file:
            landscape_grp = hdf_file.require_group('landscapes')
            grp = landscape_grp.require_group(name)
            dset = grp.create_dataset('zValues', data=self.zValues)
            dset.attrs['name'] = self.name
            dset.attrs['label'] = self.label

            if isinstance(self.limits, dict):
                dims_grp = grp.require_group('ranges')
                for axis in ['x', 'y', 'z']:
                    if axis in self.limits:
                        dims_grp.create_dataset(axis, data=self.limits[axis])

            grp.create_dataset('xbins', data=self.xbins)
            grp['xbins'].make_scale(xlabel)
            dset.dims[0].attach_scale(grp['xbins'])
            dset.dims[0].label = xlabel

            grp.create_dataset('ybins', data=self.ybins)
            grp['ybins'].make_scale(ylabel)
            dset.dims[1].attach_scale(grp['ybins'])
            dset.dims[1].label = ylabel

            ref_dtype = h5py.special_dtype(ref=h5py.Reference)
            t_dset = grp.create_dataset('data', dset.shape, dtype=ref_dtype)
            time_grp = grp.require_group('time_group')
            for x, rows in enumerate(self.binned_data):
                for y, time in enumerate(rows):
                    if len(time) > 0:
                        t = time_grp.create_dataset(
                            f'({self.xbins[x]}, {self.ybins[y]})',
                            data=numpy.array(time),
                        )
                        t_dset[x, y] = t.ref
        return

    @classmethod
    def open(cls, filename):
        landscapes = {}
        with h5py.File(filename, 'r') as hdf_file:
            for key, group in hdf_file['/landscapes/'].items():
                print(f'Loading landscape for {key} ...')
                dset = group['zValues']
                name = dset.attrs['name']
                label = dset.attrs['label']
                xbins = group['xbins']
                ybins = group['ybins']
                t_dset = group['data']
                landscape = cls(name=name,
                                xbins=xbins,
                                ybins=ybins,
                                label=label)
                landscape.zValues = dset[...]
                for i, j in itertools.product(range(len(xbins)), range(len(ybins))):
                    if t_dset[i, j]:
                        landscape.binned_data[i][j] = hdf_file[t_dset[i, j]][...]
                if group['ranges']:
                    for axis in ['x', 'y', 'z']:
                        landscape.limits[axis] = group['ranges'][axis][...]
                landscapes[landscape.name] = landscape
            return landscapes

    @staticmethod
    def minmax(array):
        min_value, max_value = numpy.Inf, -numpy.Inf
        for item in array:
            if item != numpy.NaN:
                if item < min_value:
                    min_value = item
                if item > max_value:
                    max_value = item
        return min_value, max_value

    def energy_landscape(self, temperature=298):
        """ Perform Boltzmann inversion to get the energy landscape """
        def boltzmann_inversion(p, p_max, T=298):
            return numpy.NaN if p <= 0 else -R*T*numpy.log(p/p_max)
        boltzmann_inversion = numpy.vectorize(boltzmann_inversion)

        z_max = numpy.nanmax(self.zValues)
        energies = boltzmann_inversion(self.zValues, z_max, temperature)
        self.zValues = energies - numpy.nanmax(energies)
        return self

    @classmethod
    def landscape(cls, name, time, x_data, y_data, shape=(100, 100),
                  temperature=None, limits=None, label=None):
        if len(time) < 1 or len(x_data) < 1 or len(y_data) < 1:
            raise ValueError('Empty array provided in input')
        if len(time) != len(x_data) != len(y_data):
            raise ValueError(
                'Unequal length of arrays: time, x_data and/or y_data')
        if limits is None:
            limits = dict(x=cls.minmax(x_data),
                          y=cls.minmax(y_data))
        x_bins = numpy.linspace(start=limits['x'][0],
                                stop=limits['x'][1],
                                num=shape[0])
        y_bins = numpy.linspace(start=limits['y'][0],
                                stop=limits['y'][1],
                                num=shape[1])
        landscape = cls(name=name, xbins=x_bins, ybins=y_bins, label=label)
        landscape.add_data(time=time, x_data=x_data, y_data=y_data)
        if temperature:
            landscape.energy_landscape(temperature=temperature)
        if 'z' not in limits:
            limits['z'] = (numpy.nanmin(landscape.zValues),
                           numpy.nanmax(landscape.zValues))
        landscape.limits = limits
        return landscape

    @classmethod
    def common_landscapes(cls, data, shape=(100, 100),
                          temperature=None, limits=None):
        """ Create energy landscapes on axes with identical ranges """
        # Find the common range for all landscapes
        if not limits:
            x_range, y_range = [], []
            for _, x_data, y_data, _ in data.values():
                if len(x_data) == len(y_data):
                    x_range.append(cls.minmax(x_data))
                    y_range.append(cls.minmax(y_data))
                else:
                    raise RuntimeError("sizes of x and y data do not match")
            limits = dict(x=cls.minmax(numpy.array(x_range).flatten()),
                          y=cls.minmax(numpy.array(y_range).flatten()))
        # Create a common set of bins
        x_bins = numpy.linspace(*limits['x'], shape[0])
        y_bins = numpy.linspace(*limits['y'], shape[1])
        # Create landscapes with common limits
        landscapes = []
        z_range = []
        for name, (time, x_data, y_data, label) in data.items():
            print(f'Generating Landscape for {name}')
            landscape = cls(name=name, xbins=x_bins, ybins=y_bins, label=label)
            landscape.add_data(time=time, x_data=x_data, y_data=y_data)
            if temperature:
                landscape.energy_landscape(temperature=temperature)
            z_range.append(numpy.nanmin(landscape.zValues))
            z_range.append(numpy.nanmax(landscape.zValues))
            landscapes.append(landscape)

        if 'z' not in limits:
            limits['z'] = cls.minmax(numpy.array(z_range).flatten())
        # Update the range for the 3 directions for each landscape
        for ls in landscapes:
            ls.limits = limits
        return landscapes

    @staticmethod
    def get_layout(num_subplots):
        """ Create a aesthetic layout for large number of subplots """
        if num_subplots < 1:
            raise ValueError('No subplots')
        elif num_subplots < 4:
            return num_subplots, 1
        elif 4 < num_subplots < 9:
            return math.ceil(num_subplots/2.0), 2
        elif num_subplots == 4:
            return 2, 2
        elif num_subplots == 9:
            return 3, 3
        else:
            return 4, math.ceil(num_subplots/4.0)

    @classmethod
    def plot_landscapes(cls, landscapes,
                        title='Landscapes',
                        filename='landscape.html',
                        xlabel='x', ylabel='y', zlabel='z',
                        width=None, height=None, layout=None,
                        othrographic=False, dtick=None,
                        font_family='Courier New, monospace', font_size=None):
        """ Make 2 subplots """
        subtitles = [ls.label for ls in landscapes]
        if layout:
            columns, rows = layout
        else:
            columns, rows = cls.get_layout(len(landscapes))
        fig = plotly.subplots.make_subplots(
            rows=rows, cols=columns, shared_xaxes=True, shared_yaxes=True,
            subplot_titles=subtitles,
            specs=[[{'is_3d': True}] * columns] * rows,
            horizontal_spacing=0.02, vertical_spacing=0.05)
        # adding surfaces to subplots.
        axes_indices = itertools.product(range(1, rows + 1), range(1, columns + 1))
        for landscape, (current_row, current_column) in zip(landscapes, axes_indices):
            x, y = numpy.meshgrid(landscape.xbins,
                                  landscape.ybins,
                                  indexing='ij')
            z = landscape.zValues
            current_trace = dict(
                type='surface', x=x, y=y, z=z,
                colorscale='Cividis',
                cmin=landscape.limits['z'][0],
                cmax=landscape.limits['z'][1],
                colorbar=dict(
                    title=dict(
                        text=zlabel,
                        font=dict(size=font_size),
                        side='right',
                    ),
                    len=0.5,
                    thickness=50,
                    tickfont=dict(size=font_size)
                ),
                #  showscale=False,
            )
            fig.append_trace(current_trace, current_row, current_column)
        # Show contour lines
        fig.update_traces(contours_z=dict(show=True, usecolormap=True,
                          highlightcolor="limegreen", project_z=True))
        fig.layout.update(title=title)

        for i, ls in enumerate(landscapes, start=1):

            # Rotate Camera
            fig.layout[f'scene{i}'].camera = dict(
                center=dict(x=0, y=0, z=-0.17),
                eye=dict(x=-1.25, y=1.25, z=1)
            )

            fig['layout'][f'scene{i}'].update(
                xaxis_title=xlabel,
                yaxis_title=ylabel,
                zaxis_title=zlabel,
            )
            # fig['layout'][f'scene{i}'].colorbar.update(len=0.5)
            if othrographic:
                fig.layout[f'scene{i}'].camera.projection.update(
                    type='orthographic')
            for axis in ['x', 'y', 'z']:
                if axis in ls.limits:
                    fig['layout'][f'scene{i}'][axis + 'axis'].update(
                        range=ls.limits[axis])

        if width:
            fig.layout.update(width=width)
        if height:
            fig.layout.update(height=height)

        if dtick:
            for i in range(1, len(landscapes) + 1):
                for axis in ['x', 'y', 'z']:
                    fig['layout'][f'scene{i}'][axis + 'axis'].update(
                        dtick=dtick[axis]
                    )
        if font_size:
            fig.layout.font.update(
                family=font_family,
                # Fonts size for ticks do not match other fonts
                size=font_size - 10 if font_size > 20 else 12,
            )
            for i in range(1, len(landscapes) + 1):
                # fig['layout'][f'scene{i}'].colorbar.tickfont.update(
                #     size=font_size)
                for axis in ['x', 'y', 'z']:
                    fig['layout'][f'scene{i}'][axis + 'axis'].title.font.update(
                        family=font_family,
                        size=font_size,
                    )
            for annotation in fig.layout.annotations:
                annotation.font = dict(family=font_family,
                                       size=font_size)

        div = py.plot(fig, include_plotlyjs=False, output_type='div')
        # retrieve the div id for the figure
        div_id = div.split('=')[1].split()[0].replace("'", "").replace('"', '')

        # Header for the output HTML file
        header = ('<!DOCTYPE html><html>\n<head><metacharset="utf-8"/>\n'
                  '<script src="https://cdn.plot.ly/plotly-latest.min.js">'
                  '</script>\n</head>\n<body>\n')

        # Custom JS code to sync the view of the plots
        if othrographic:
            js = '''
<script>
var gd = document.getElementById('{div_id}');
var isUnderRelayout = false;
var oldCamera, newCamera;
var oldZoom, newZoom;
var scenes = [];
for (propertyName in gd.layout) {{
    if ( propertyName.startsWith('scene') ) {{
        scenes.push(propertyName);
    }}
}}
gd.on('plotly_relayout', () => {{
if ( !isUnderRelayout ) {{
    for (i = 0; i < scenes.length; i++) {{
        if (gd.layout[scenes[i]].camera != oldCamera && !undefined) {{
            newCamera = gd.layout[scenes[i]].camera
        }}
        if (gd.layout[scenes[i]].aspectratio != oldZoom && !undefined) {{
            newZoom = gd.layout[scenes[i]].aspectratio
        }}
    }}
    for (j = 0; j < scenes.length; j++) {{
        Plotly.relayout(gd, scenes[j] + '.camera', newCamera)
        Plotly.relayout(gd, scenes[j] + '.aspectratio', newZoom)
        .then(() => {{ isUnderRelayout = false }}  );
    }}
    oldCamera = newCamera
    oldZoom = newZoom
}}
isUnderRelayout = true;
}})
</script>'''.format(div_id=div_id)
        else:
            js = '''
<script>
var gd = document.getElementById('{div_id}');
var isUnderRelayout = false;
var oldCamera, newCamera;
var scenes = [];
for (propertyName in gd.layout) {{
    if ( propertyName.startsWith('scene') ) {{
        scenes.push(propertyName);
    }}
}}
gd.on('plotly_relayout', () => {{
if ( !isUnderRelayout ) {{
    for (i = 0; i < scenes.length; i++) {{
        if (gd.layout[scenes[i]].camera != oldCamera && !undefined) {{
            newCamera = gd.layout[scenes[i]].camera
        }}
    }}
    for (j = 0; j < scenes.length; j++) {{
        Plotly.relayout(gd, scenes[j] + '.camera', newCamera)
        .then(() => {{ isUnderRelayout = false }}  );
    }}
    oldCamera = newCamera
}}
isUnderRelayout = true;
}})
</script>'''.format(div_id=div_id)
        # merge everything
        div = header + div + js + '\n</body>\n</html>'
        print(div, file=open(filename, "w"))


def main():
    pass


if __name__ == "__main__":
    pass
