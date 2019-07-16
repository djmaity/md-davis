import argparse
import numpy
import h5py
import mdtraj


def add_dihedrals(group, trajectory, structure, chunk=1000, atoms=None):
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

    phi_dset.attrs['indices'] = phi_indices
    psi_dset.attrs['indices'] = psi_indices
    omega_dset.attrs['indices'] = omega_indices

    return numpy.hstack( time )


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
        time = add_dihedrals(dih_grp,
                             trajectory=args.trajectory,
                             structure=args.structure,
                             chunk=args.chunk,
                             atoms=args.atoms,
        )
        datafile.attrs['time'] = time
        datafile.attrs['time_unit'] = 'picosecond'


if __name__ == "__main__":
    main()
