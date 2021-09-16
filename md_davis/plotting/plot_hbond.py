import operator

import click
import numpy
import pandas
import plotly
import plotly.graph_objects as go


def text(atoms):
    output = []
    for atom in atoms:
        if atom:
            output.append('-'.join([str(_) for _ in atom]))
        else:
            output.append('')
    return output


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(name='plot_hbond', context_settings=CONTEXT_SETTINGS)
@click.option('-t', '--total_frames', type=int, required=True,
              help='Total number of frames in the trajectory used for H-bond/Contact evaluation.')
@click.option('-c', '--cutoff', type=float,
              help='Cutoff for percentage or number of Hbonds above which to include in the plot.')
@click.option('-o', '--output', default='hbond_matrix.html', type=click.Path(),
              help='Output filename')
@click.option('--percent', is_flag=True, help='Plot the percentages')
@click.option('--title', help='Title for the plot')
@click.argument('hbond_file')
def main(hbond_file, total_frames=None, cutoff=None, output='hbond_matrix.html', percent=False, title=None):
    """ H-bond or contact dataframe obtained from md_davis hbond """
    df = pandas.read_pickle(hbond_file)

    zmin = 0
    zmax = total_frames
    if percent:
        df['Count'] = df['Count'] * 100 / total_frames  # Convert to percentage of frames
        if cutoff:
            df = df[df['Count'] > cutoff]
            zmin = cutoff
        zmax = 100
    else:
        if cutoff:
            df = df[df['Count'] > cutoff]
            zmin = cutoff

    atoms1 = list(df.groupby(['Chain1', 'Residue1', 'ResSeq1', 'Atom1']).groups.keys())
    atoms1 = sorted(atoms1, key=operator.itemgetter(0, 2, 1))
    atom_group2 = df.groupby(['Chain2', 'Residue2', 'ResSeq2', 'Atom2'])
    atoms2 = list(atom_group2.groups.keys())
    atoms2 = sorted(atoms2, key=operator.itemgetter(0, 2, 1))

    matrix = numpy.empty(shape=(len(atoms1), len(atoms2)))
    matrix[:] = numpy.NaN
    for group, atoms in atom_group2:
        j = atoms2.index(group)
        for index, row in atoms.iterrows():
            i = atoms1.index((row['Chain1'], row['Residue1'], row['ResSeq1'], row['Atom1']))
            matrix[i, j] = row['Count']

    fig = go.Figure(data=go.Heatmap(
        z=matrix,
        x=text(atoms1),
        y=text(atoms2),
        colorscale='Greys', zmin=zmin, zmax=zmax,
        transpose=True,
    ))

    if title is not None:
        fig.update_layout(
            title=title,
        )
    fig.write_html(file=output, auto_open=True)


if __name__ == '__main__':
    main()
