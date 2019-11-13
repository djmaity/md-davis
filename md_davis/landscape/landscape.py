import math
import numpy
import h5py
import itertools

import plotly.offline as py
import plotly.graph_objs as go
import plotly.subplots

# Exact value after redifinition of SI units in May 2019
R = 0.00831446261815324  # Energy in KJ / mol


class Landscape(object):

    def __init__(self, name, xBins, yBins, label=None):
        self.name = name
        self.label = label if label else name
        self.time_bins = [[[] for _ in range(len(yBins))] for _ in range(len(xBins))]
        self.zValues = numpy.zeros((len(xBins), len(yBins)))
        self.xBins = numpy.array(xBins)
        self.yBins = numpy.array(yBins)
        self.dims = {'x': None, 'y': None, 'z': None}

    def __repr__(self):
        return 'Name:  ' + self.name + '\n' \
                         + 'Label: ' + self.label + '\n' \
                         + self.dims.__repr__() + '\n\n' \
                         + 'Z values:\n' + self.zValues.__repr__() +'\n\n' \
                         + 'X bins:\n' + self.xBins.__repr__() +'\n\n' \
                         + 'Y bins:\n' + self.yBins.__repr__() +'\n\n'

    def __str__(self):
        return 'Name:  ' + self.name + '\n' \
                         + 'Label: ' + self.label + '\n' \
                         + self.dims.__str__() + '\n\n' \
                         + 'Z values:\n' + self.zValues.__str__() +'\n\n' \
                         + 'X bins:\n' + self.xBins.__str__() +'\n\n' \
                         + 'Y bins:\n' + self.yBins.__str__() +'\n\n'

    def __iter__(self):
        for i, rows in enumerate(self.time_bins):
            for j, time in enumerate(rows):
                yield i, j, time

    def __len__(self):
        return len(self.xBins)

    def add_data(self, time, x_data, y_data):
        """ Create a 2D counts and add times to the corresponding bins """
        dx = numpy.digitize(x_data, self.xBins)
        dy = numpy.digitize(y_data, self.yBins)
        for t in range(len(time)):
            self.zValues[dx[t] - 1][dy[t] - 1] += 1
            self.time_bins[dx[t] - 1][dy[t] - 1].append(time[t])
        return self

    def save(self, filename, name, xlabel='', ylabel=''):
        print(f'Saving {name} ...')
        with h5py.File(filename, 'a') as hdf_file:
            landscape_grp = hdf_file.require_group('landscapes')
            grp = landscape_grp.create_group(name)
            dset = grp.create_dataset('zValues', data=self.zValues)
            dset.attrs['name'] = self.name
            dset.attrs['label'] = self.label

            if isinstance(self.dims, dict):
                dims_grp = grp.create_group('dimensions')
                for axis in ['x', 'y', 'z']:
                    if axis in self.dims:
                        dims_grp.create_dataset(axis, data=self.dims[axis])

            grp.create_dataset('xBins', data=self.xBins)
            grp['xBins'].make_scale(xlabel)
            dset.dims[0].attach_scale(grp['xBins'])
            dset.dims[0].label = xlabel

            grp.create_dataset('yBins', data=self.yBins)
            grp['yBins'].make_scale(ylabel)
            dset.dims[1].attach_scale(grp['yBins'])
            dset.dims[1].label = ylabel

            ref_dtype = h5py.special_dtype(ref=h5py.Reference)
            t_dset = grp.create_dataset('time_bins', dset.shape, dtype=ref_dtype)
            time_grp = grp.create_group('time_group')
            for x, rows in enumerate(self.time_bins):
                for y, time in enumerate(rows):
                    print(x, y, self.xBins[x], self.yBins[y], len(time), self.zValues[x, y] )
                    if len(time) > 0:
                        t = time_grp.create_dataset(f'({self.xBins[x]}, {self.yBins[y]}, {x}, {y})',
                            data=numpy.array(time),
                        )
                        t_dset[x, y] = t.ref
        return

    @classmethod
    def open(cls, filename):
        landscapes = []
        dimensions = None
        with h5py.File(filename, 'r') as hdf_file:
            for key, group in hdf_file['/landscapes/'].items():
                print(f'Loading landscape for {key} ...')
                dset = group['zValues']
                name = dset.attrs['name']
                label = dset.attrs['label']
                xBins = group['xBins']
                yBins = group['yBins']
                t_dset = group['time_bins']
                landscape = cls(name=name, xBins=xBins, yBins=yBins, label=label)
                landscape.zValues = dset[...]
                for i, j in itertools.product(range(len(xBins)), range(len(yBins))):
                    if t_dset[i,j]:
                        landscape.time_bins[i][j] = hdf_file[t_dset[i,j]][...]
                if group['dimensions']:
                    for axis in ['x', 'y', 'z']:
                        landscape.dims[axis] = group['dimensions'][axis][...]
                landscapes.append(landscape)
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
        return (min_value, max_value)

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
    def landscape(cls, name, time, x_data, y_data, shape=(100, 100), temperature=None, label=None):
        if len(time) < 1 or len(x_data) < 1 or len(y_data) < 1:
            raise ValueError('Empty array provided in input')
        if len(time) != len(x_data) != len(y_data):
            raise ValueError('Unequal length of arrays: time, x_data and/or y_data')
        dimensions = dict(x=cls.minmax(x_data),
                          y=cls.minmax(y_data),
        )
        x_bins = numpy.linspace(start = dimensions['x'][0],
                                stop = dimensions['x'][1],
                                num = shape[0],
        )
        y_bins = numpy.linspace(start = dimensions['y'][0],
                                stop = dimensions['y'][1],
                                num = shape[1],
        )
        landscape = cls(name=name, xBins=x_bins, yBins=y_bins, label=label)
        landscape.add_data(time=time, x_data=x_data, y_data=y_data)
        if temperature:
            landscape.energy_landscape(temperature=temperature)
        dimensions['z'] = (numpy.nanmin(landscape.zValues),
                           numpy.nanmax(landscape.zValues),
        )
        landscape.dims = dimensions
        return landscape

    @classmethod
    def common_landscapes(cls, data, shape=(100, 100), temperature=None):
        """ Create enegry landscapes on axes with identical ranges """
        # Find the common range for all landscapes
        x_range, y_range = [], []
        for _, x_data, y_data, _ in data.values():
            if len(x_data) == len(y_data):
                x_range.append(cls.minmax(x_data))
                y_range.append(cls.minmax(y_data))
            else:
                raise RuntimeError
        dimensions = dict(x=cls.minmax(numpy.array(x_range).flatten()),
                          y=cls.minmax(numpy.array(y_range).flatten()),
        )
        # Create a common set of bins
        x_bins = numpy.linspace(*dimensions['x'], shape[0])
        y_bins = numpy.linspace(*dimensions['y'], shape[1])
        # Create landscapes with common dimensions
        landscapes = []
        z_range = []
        for name, (time, x_data, y_data, label) in data.items():
            print(f'Generating Landscape for {name}')
            landscape = cls(name=name, xBins=x_bins, yBins=y_bins, label=label)
            landscape.add_data(time=time, x_data=x_data, y_data=y_data)
            if temperature:
                landscape.energy_landscape(temperature=temperature)
            z_range.append(numpy.nanmin(landscape.zValues))
            z_range.append(numpy.nanmax(landscape.zValues))
            landscapes.append(landscape)
        dimensions['z'] = cls.minmax(numpy.array(z_range).flatten())
        # Update the range for the 3 directions for each landscape
        for ls in landscapes:
            ls.dims = dimensions
        return landscapes

    @staticmethod
    def get_layout(num_subplots):
        """ Create a aethetically pleasing layout for large number of subplots """
        if num_subplots < 1:
            raise ValueError
        elif num_subplots < 4:
            return (num_subplots, 1)
        elif 4 < num_subplots < 9:
            return ( math.ceil(num_subplots/2.0) , 2)
        elif num_subplots == 4:
            return (2, 2)
        elif num_subplots == 9:
            return (3, 3)
        else:
            return (4, math.ceil(num_subplots/4.0) )

    @classmethod
    def plot_landscapes(cls, landscapes,
                        title='Landscapes',
                        filename='landscape.html',
                        xlabel='x', ylabel='y', zlabel='z',
                        width=None, height=None, othrographic=False, font_size=None):
        """ Make 2 subplots """
        subtitles = [ls.label for ls in landscapes]
        columns, rows = cls.get_layout(len(landscapes))
        fig = plotly.subplots.make_subplots(rows=rows,
                                            cols=columns,
                                            shared_xaxes=True,
                                            shared_yaxes=True,
                                            subplot_titles=subtitles,
                                            specs=[ [{'is_3d': True}] * columns ] * rows,
                                            horizontal_spacing = 0.02,
                                            vertical_spacing = 0.05,
                                            )
        # adding surfaces to subplots.
        axes_indices = itertools.product(range(1, rows + 1), range(1, columns + 1))
        for landscape, (current_row, current_column) in zip(landscapes, axes_indices):
            x, y = numpy.meshgrid(landscape.xBins, landscape.yBins)
            print(x, y)
            z = landscape.zValues
            current_trace = dict(type='surface', x=x, y=y, z=z,
                                 colorscale='Cividis',
                                 cmin=landscape.dims['z'][0],
                                 cmax=landscape.dims['z'][1],
                                #  showscale=False,
            )
            fig.append_trace(current_trace, current_row, current_column)
        # Show contour lines
        fig.update_traces(contours_z=dict(show=True, usecolormap=True,
                          highlightcolor="limegreen", project_z=True))
        fig.layout.update(title=title)

        for i, ls in enumerate(landscapes, start=1):
            fig['layout'][f'scene{i}'].update(
                xaxis_title=xlabel,
                yaxis_title=ylabel,
                zaxis_title=zlabel,
            )
            if othrographic:
                fig.layout[f'scene{i}'].camera.projection.update(type='orthographic')
            for axis in ['x', 'y', 'z']:
                if axis in ls.dims:
                    fig['layout'][f'scene{i}'][axis + 'axis'].update(range=ls.dims[axis])
                if font_size:
                    fig['layout'][f'scene{i}'][axis + 'axis'].title.font.update(
                        family='Courier New, monospace',
                        size=font_size,
                    )
        if width:
            fig.layout.update(width=width)
        if height:
            fig.layout.update(height=height)
        if font_size:
            for annotation in fig.layout.annotations:
                annotation.font = dict(family='Courier New, monospace',
                                       size=font_size)
        div = py.plot(fig, include_plotlyjs=False, output_type='div')
        # retrieve the div id for the figure
        div_id = div.split('=')[1].split()[0].replace("'", "").replace('"', '')

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
        div = '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>' + div + js
        print(div, file=open(filename, "w"))


def main():
    pass


if __name__ == "__main__":
    pass
