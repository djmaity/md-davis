import argparse
import numpy
import h5py
import scipy.stats as statistics


def split_array(array, lengths):
    assert sum(lengths) == len(array)
    start = 0
    for end in lengths:
        end = start + end
        sub_array = array[start:end]
        start = end
        yield sub_array


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-f', '--file', required=True,
                        help='HDF5 file')
    parser.add_argument('-b', '--begin', default=0, type=float,
                        help='Start frame')
    parser.add_argument('-e', '--end', default=1000, type=float,
                        help='End frame')
    parser.add_argument('-s', '--step', default=200, type=float,
                        help='Number of frames per ns')
    args = parser.parse_args()

    start = int(args.begin * args.step)
    end = int(args.end * args.step)

    with h5py.File(args.file, 'a') as hdf5_file:
        print(hdf5_file)
        # start = int(float(hdf5_file['rmsf'].attrs['begin'][:-3]) * args.step)
        # frames, _ = numpy.shape(hdf5_file[f'/dihedrals/phi'])
        # step = int((frames - 1) / 1000)

        dih_sd_group = hdf5_file.create_group('dihedral_standard_deviation')
        dih_sd_group.attrs['unit'] = 'radians'
        dih_sd_group.attrs['begin'] = f'{args.begin} ns'
        dih_sd_group.attrs['end'] = f'{args.end} ns'

        lengths = [len(_) - 1  for _ in hdf5_file.attrs['sequence'].split('/')]
        data_type =  numpy.dtype([("phi", numpy.float),
                                  ("psi", numpy.float),
                                  ("omega", numpy.float)
                                ])

        for ch, length in enumerate(lengths):
            dih_sd_group.create_dataset(f'chain {ch}', (length + 1,), dtype=data_type)

        for angle in ['phi', 'psi', 'omega']:
            _, num_angles = numpy.shape(hdf5_file[f'/dihedrals/{angle}'])
            assert sum(lengths) == num_angles
            circ_sd = statistics.circstd(hdf5_file[f'/dihedrals/{angle}'][start:end], axis=0, high=numpy.pi)
            for ch, chain_dih_sd in enumerate(split_array(circ_sd, lengths)):
                dset = hdf5_file[f'/dihedral_standard_deviation/chain {ch}']
                if angle == 'phi':
                    dset[f'{angle}', 1:] = chain_dih_sd
                else:
                    dset[f'{angle}', :-1] = chain_dih_sd


if __name__ == "__main__":
    main()
