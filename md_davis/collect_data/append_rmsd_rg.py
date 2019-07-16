import argparse
import numpy
import h5py
import mdtraj


def new_dataset(group, dataset, rmsd_file, rg_file):
    """ Make RMSD & Rg dataset in HDF5 file """
    print(f'Collecting data from {rmsd_file} and {rg_file} into HDF5 file.')
    rmsd = numpy.loadtxt(rmsd_file, comments=('#', '@'), dtype=numpy.single)
    rg = numpy.loadtxt(rg_file, comments=('#', '@'), dtype=numpy.single)
    assert numpy.array_equal(rmsd[:, 0], rg[:, 0]), \
        "The times in RMSD file do not match those in radius of gyration file"
    data = numpy.core.records.fromarrays(numpy.hstack([rmsd, rg[:, 1:]]).T, 
        names='time, rmsd, rg, rg_x, rg_y, rg_z')
    dset = group.create_dataset(dataset, data=data)
    dset.attrs['time_unit'] = 'picosecond'
    dset.attrs['unit'] = 'nanometer'


def add_rmsd_rg(hdf5file, directory, prefix):
    """ Add RMSD & Rg Data HDF5 file """
     # RMSD w.r.t. last frame
    rmsd_file_all_atom = f'{directory}/rmsd_rg/{prefix}_processed-rmsd.xvg'
    rmsd_file_backbone = f'{directory}/rmsd_rg/{prefix}_backbone-rmsd.xvg'
    rmsd_file_c_alpha = f'{directory}/rmsd_rg/{prefix}_c-alpha-rmsd.xvg'

    rg_file_all_atom = f'{directory}/rmsd_rg/{prefix}_processed-rg.xvg'
    rg_file_backbone = f'{directory}/rmsd_rg/{prefix}_backbone-rg.xvg'
    rg_file_c_alpha = f'{directory}/rmsd_rg/{prefix}_c-alpha-rg.xvg'
    
    with h5py.File(hdf5file, 'a') as datafile:
        rmsd_rg_group = datafile.create_group('rmsd_rg')
        rmsd_rg_group.attrs['comment'] = 'RMSD is calculated with respect to last frame'
        new_dataset(rmsd_rg_group, 'all_atom', rmsd_file_all_atom, rg_file_all_atom)
        new_dataset(rmsd_rg_group, 'backbone', rmsd_file_backbone, rg_file_backbone)
        new_dataset(rmsd_rg_group, 'c-alpha', rmsd_file_c_alpha, rg_file_c_alpha) 


def double_to_single(hdf5file):
    """ Change the precision from double to single """
    with h5py.File(hdf5file, 'a') as datafile:
        datafile['dihedrals/phi'].attrs['indices'] = datafile['dihedrals/phi'].attrs['indices'].astype(numpy.int32)
        datafile['dihedrals/psi'].attrs['indices'] = datafile['dihedrals/psi'].attrs['indices'].astype(numpy.int32)
        datafile['dihedrals/omega'].attrs['indices'] = datafile['dihedrals/omega'].attrs['indices'].astype(numpy.int32)


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-f', '--file', required=True,
                        help='HDF5 file')
    parser.add_argument('-d', '--directory', required=True,
                        help='Directory containing RMSD & Rg files')
    parser.add_argument('-p', '--prefix', required=True,
                        help='Prefix of the files parsed')
    args = parser.parse_args()

    double_to_single(args.file)
    add_rmsd_rg(hdf5file=args.file, directory=args.directory, prefix=args.prefix)

if __name__ == "__main__":
    main()
