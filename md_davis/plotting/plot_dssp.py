import argparse
import collections
import h5py
import plotly.express as px
import numpy

structure_info = collections.namedtuple('structure_info', ['label', 'html', 'color'])

secondary_structure = collections.OrderedDict([
    ('H', structure_info(label='α-helix', html='α helix', color='rgb(255, 0, 0)') ),
    ('G', structure_info(label='3_10-helix', html='3<sub>10</sub> helix', color='rgb(120, 0, 100)') ),
    ('I', structure_info(label='π-helix', html='π helix', color='rgb(255, 15, 139)') ),
    ('E', structure_info(label='β strand', html='β strand', color='rgb(0, 0, 255)') ),
    ('B', structure_info(label='β-bridge', html='β bridge', color='rgb(15, 125, 255)') ), 
    ('T', structure_info(label='Turn', html='Turn', color='rgb(0, 200, 0)') ),
    ('S', structure_info(label='Bend', html='Bend', color='rgb(120, 200, 0)') ),
    ('~', structure_info(label='Loop', html='Loop', color='rgb(159, 255, 15)') ),
])


def dssp2color(dssp_char):
    if isinstance(dssp_char, bytes):
        dssp_char = dssp_char.decode('ascii')
    rgb = secondary_structure[dssp_char].color[4:-1].split(',')
    return [int(_) for _ in rgb]


def main():
    parser = argparse.ArgumentParser(description='Plot DSSP data from.')
    parser.add_argument('hdf_file')
    parser.add_argument('-c', '--chain', type=int)
    args = parser.parse_args()

    time_start, time_end, time_step = 0, None, 10
    res_start, res_end, res_step = 0, None, 1
    chain = args.chain

    with h5py.File(args.hdf_file, 'r') as hdf_file:
        dataset = hdf_file[f'secondary_structure/dssp_data/chain {chain}']
        dataset = dataset[time_start : time_end : time_step,
                            res_start :  res_end :  res_step].T
        output = numpy.empty(shape=(*dataset.shape, 3), dtype=numpy.uint8)
        for i, line in enumerate(dataset):
            for j, char in enumerate(line):
                output[i, j] = dssp2color(char)
        fig = px.imshow(output)
        # fig.update_layout(width=1000)
        fig.show()


if __name__ == "__main__":
    main()
