import argparse
import numpy
import h5py
import mdtraj
import scipy.stats as statistics


def split_chains(array, lengths):
    """ """
    assert sum(lengths) == len(array)
    start = 0
    for end in lengths:
        end = start + end
        sub_array = array[start:end]
        start = end
        yield sub_array


def get_dihedrals(group, trajectory, structure, chunk=1000, atoms=None):
    """ Evaluate the dihedral angles using MDTraj and store it as HDF5 dataset """
    # if atoms:
    #     trj_iterator = mdtraj.iterload(trajectory, top=structure, chunk=chunk, atom_indices=atoms)
    # else:
    #     trj_iterator = mdtraj.iterload(trajectory, top=structure, chunk=chunk)
    trj_iterator = mdtraj.iterload(trajectory, top=structure, chunk=chunk, atom_indices=atoms)

    first_trj_chunk = next(trj_iterator)
    phi_indices, phi = mdtraj.compute_phi(first_trj_chunk)
    psi_indices, psi = mdtraj.compute_psi(first_trj_chunk)
    omega_indices, omega = mdtraj.compute_omega(first_trj_chunk)

    time, phi_array, psi_array, omega_array = [first_trj_chunk.time], [phi], [psi], [omega]
    for trj_chunk in trj_iterator:
        print(f'Calculating dihedrals for trajectory between {trj_chunk.time[0]} and {trj_chunk.time[-1]} ps')
        time.append(trj_chunk.time)
        phi_array.append( mdtraj.compute_phi(trj_chunk)[1] )
        psi_array.append( mdtraj.compute_psi(trj_chunk)[1] )
        omega_array.append( mdtraj.compute_omega(trj_chunk)[1] )

    phi_dset = group.create_dataset( 'phi', data=numpy.vstack( phi_array ) )
    psi_dset = group.create_dataset( 'psi', data=numpy.vstack( psi_array ) )
    omega_dset = group.create_dataset( 'omega', data=numpy.vstack( omega_array ) )

    phi_dset.attrs['indices'] = phi_indices.astype(numpy.int32)
    psi_dset.attrs['indices'] = psi_indices.astype(numpy.int32)
    omega_dset.attrs['indices'] = omega_indices.astype(numpy.int32)
    return numpy.hstack( time )


def get_dihedral_sd(hdf5_file, begin=0, end=1000, step=200):
    """
        hdf5_file must contain dihedrals and sequence

        begin = Start frame in ns
        end = End frame in ns
        step = Number of frames per ns
    """
    start = int(begin * step)
    end = int(end * step)

    dih_sd_group = hdf5_file.create_group('dihedral_standard_deviation')
    dih_sd_group.attrs['unit'] = 'radians'
    dih_sd_group.attrs['begin'] = f'{begin} ns'
    dih_sd_group.attrs['end'] = f'{end} ns'

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
        for ch, chain_dih_sd in enumerate(split_chains(circ_sd, lengths)):
            dset = hdf5_file[f'/dihedral_standard_deviation/chain {ch}']
            if angle == 'phi':
                dset[f'{angle}', 1:] = chain_dih_sd
            else:
                dset[f'{angle}', :-1] = chain_dih_sd


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-t', '--trajectory', required=True,
                        help='Trajectory file')
    parser.add_argument('-s', '--structure', required=True,
                        help='Structure file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output HDF5 file')
    parser.add_argument('-c', '--chunk', default=1000,
                        help='Number of frames to load at once')
    parser.add_argument('-a', '--atoms', default=None,
                        help='Number of frames to load at once')
    args = parser.parse_args()

    with h5py.File(args.output, 'a') as datafile:
        dih_grp = datafile.create_group('dihedrals')
        time = get_dihedrals(dih_grp,
                             trajectory=args.trajectory,
                             structure=args.structure,
                             chunk=args.chunk,
                             atoms=args.atoms,
        )
        datafile.attrs['time'] = time
        datafile.attrs['time_unit'] = 'picosecond'

        get_dihedral_sd(datafile)


if __name__ == "__main__":
    main()
