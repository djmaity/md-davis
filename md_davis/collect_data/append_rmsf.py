import argparse
import numpy
import h5py


def split_increasing(array, use_col=0):
    """ Split an array when it is not increasing """
    array = numpy.array(array)
    if array.ndim == 1:
        column = array
    elif array.ndim == 2:
        column = array[use_col]
    else:
        raise RuntimeError("Can not handle data with dimension other than 1 or 2")
    # Check if following value is greater
    where_greater = numpy.greater(column[:-1], column[1:])
    # Use the indices where the above is True
    split_indices = numpy.flatnonzero(where_greater)
    # to split the array
    split_arrays = numpy.split(array, split_indices + 1, axis=array.ndim - 1)
    return list(map(numpy.transpose, split_arrays))

# # Tests
# for arr in ( [1, 2, 3, 4], [4, 3, 2, 1], [1, 2, 1, 2], [] ):
#     print(arr)
#     for i in split_increasing(arr):
#         print(i)
#     print()
# exit()


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-f', '--file', required=True,
                        help='HDF5 file')
    parser.add_argument('-r', '--rmsf', required=True,
                        help='RMSF file from gmx rmsf')
    parser.add_argument('-b', '--begin', default=0, type=float, help='Begin time in ns')
    parser.add_argument('-e', '--end', default=1000, type=float, help='End time in ns')
    args = parser.parse_args()

    with h5py.File(args.file, 'a') as hdf_file:
        rmsf = numpy.loadtxt(args.rmsf, comments=('#', '@'), dtype=numpy.single)
        rmsf_split_by_chains = split_increasing(rmsf.T)
        rmsf_group = hdf_file.create_group('rmsf')
        for ch, chain_rmsf in enumerate(rmsf_split_by_chains):
            rmsf_group.create_dataset(f'chain {ch}', data=chain_rmsf)
        rmsf_group.attrs['unit'] = 'nanometer'
        rmsf_group.attrs['begin'] = f'{args.begin} ns'
        rmsf_group.attrs['end'] = f'{args.end} ns'


if __name__ == "__main__":
    main()
