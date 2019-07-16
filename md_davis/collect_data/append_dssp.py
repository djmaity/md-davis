import argparse
import numpy
import h5py
import collections
import statistics

SECSTR_CODES = {'H':'α-helix',
                'G':'3_10-helix',
                'I':'π-helix',
                'B':'β-bridge',
                'E':'β strand',
                'T':'Turn',
                'S':'Bend',
                '~':'Loop',
}

def parse_dat(filename):
    """ Parse the DSSP data file """
    secondary_structure = []
    with open(filename, 'r') as dat_file:
        next(dat_file)
        for line in dat_file:
            secondary_structure.append(line.strip().split('='))
    output = []
    for chain in zip(*secondary_structure):
        assert len(set([len(_) for _ in chain])) <= 1 , \
            "The length of the secondary structure strings are not equal. Please check your input .dat file."
        output.append(numpy.array([list(_) for _ in chain], dtype=f'S1'))
    return output


def ss_count_per_residue(chain_dssp_data):
    ss_dtype = [(_, numpy.uint32) for _ in SECSTR_CODES.keys()]
    output = numpy.zeros(len(chain_dssp_data.T), dtype=ss_dtype)
    for resSeq, residue in enumerate(chain_dssp_data.T):
        counts = numpy.unique(residue, return_counts=True)
        for structure, count in zip(*counts):
            output[structure.decode('ascii')][resSeq] = count
    return output


def main():
    """ Parse and plot the DSSP data """
    parser = argparse.ArgumentParser(description='Add the secondary structure'
        'data exhibited by each residue in the output of gmx do_dssp')
    parser.add_argument('-f', '--file', required=True,
                        help='HDF5 file')
    parser.add_argument('-d', '--dssp', required=True,
                        help='DSSP .dat file from gmx do_dssp')
    args = parser.parse_args()

    with h5py.File(args.file, 'a') as datafile:
        dssp_group = datafile.create_group('secondary_structure/dssp_data')
        ss_counts_group = datafile.create_group('secondary_structure/counts_per_residue')
        dssp_data = parse_dat(args.dssp)
        for ch, chain_dssp_data in enumerate(dssp_data):
            dssp_group.create_dataset(f'chain {ch}', data=chain_dssp_data)
            ss_counts_group.create_dataset(f'chain {ch}',
                data=ss_count_per_residue(chain_dssp_data)
            )


if __name__ == "__main__":
    main()

