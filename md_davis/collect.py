#! /usr/bin/env python
"""
Create an HDF file for collection analysis data

Usage:
  md_davis collect [(--all_atom_rmsd FILE --all_atom_rg FILE)]
                   [(--backbone_rmsd FILE --backbone_rg FILE)]
                   [(--c_alpha_rmsd FILE --c_alpha_rg FILE)]
                   [(--rmsf "FILES" <begin> <end>)]
                   [(--trajectory FILE --structure FILE -c <int> -a <list>)]
                   [--dipoles FILE]
                   [--ss FILE]
                   [--info JSON]
                   HDF_FILE

  md_davis collect -h | --help


Arguments:
  <begin>               Starting time for RMSF calculation in nanoseconds
  <end>                 End time for RMSF calculation in nanoseconds

Options:
  --all_atom_rmsd FILE      Add All atom RMSD into HDF_FILE
  --backbone_rmsd FILE      Add Backbone RMSD into HDF_FILE
  --c_alpha_rmsd FILE       Add c_alpha RMSD into HDF_FILE
  --all_atom_rg FILE        Add All atom Rg into HDF_FILE
  --backbone_rg FILE        Add Backbone Rg into HDF_FILE
  --c_alpha_rg FILE         Add c_alpha Rg into HDF_FILE

  -r, --rmsf FILES          Add data from RMSF file into the HDF_FILE
  -d, --dipoles FILE        Add dipole_moment into the HDF_FILE
  -t, --trajectory FILE     MD trajectory file
  -s, --structure FILE      Structure file, preferably in PDB format
  -c, --chunk <int>         Number of frames to read at once
  -a, --atoms <list>        Atoms to use for dihedral calculation
  --ss FILE                 Add secondary sctructure information from
                            'gmx do_dssp' into HDF_FILE

  -h, --help                 Show this screen

  -i, --info JSON           Add labels and other information as attributes
                            in the HDF_FILE

The attributes resuired as JSON with '--info' are given in the example below:

{
    "label": "MD Simulation",
    "short_label": "MD",
    "html": "<i>MD Simulation</i>",
    "short_html": "<i>MD Simulation</i>",
    "protein": "My Protein",
    "scientific_name": "Homo sapiens",
    "common_name": "Human",
    "sequence": "ACDEFGHIKLMNPQRSTVWY"
}

This information is primarily parsed to create labels for plots with
this data file. 'sequence' is required to split the data into
chains. A JSON file can also be supplied instead of a string
"""

import os
import collections
import json
import statistics
import h5py
import docopt
import mdtraj
import numpy
import scipy.stats


SECSTR_CODES = {'H':'α-helix',
                'G':'3_10-helix',
                'I':'π-helix',
                'B':'β-bridge',
                'E':'β strand',
                'T':'Turn',
                'S':'Bend',
                '~':'Loop',
}


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


def split_chains(array, lengths):
    """ """
    assert sum(lengths) == len(array)
    start = 0
    for end in lengths:
        end = start + end
        sub_array = array[start:end]
        start = end
        yield sub_array


def add_info(hdf_file, info):
    """ Add labels and other information as attributes in the HDF file """
    # TODO: Catch error if filename wrong
    if os.path.exists(info):
        attributes = json.load(open(info, 'r'))
    else:
        json_acceptable_string = info.replace("'", "\"")
        attributes = json.loads(json_acceptable_string)
    for key, value in attributes.items():
        hdf_file.attrs[key] = value


def add_rmsd_rg(hdf_file, rmsd, rg, dataset='all_atom',
    comment='RMSD is calculated with respect to last frame',
    time_unit='picosecond', unit='nanometer'):
    """ Make RMSD & Rg dataset in HDF5 file """
    group = hdf_file.require_group('rmsd_rg')
    group.attrs['comment'] = comment
    print(f'Collecting data from {rmsd} and {rg} into HDF5 file.')
    rmsd_data = numpy.loadtxt(rmsd, comments=('#', '@'), dtype=numpy.single)
    rg_data = numpy.loadtxt(rg, comments=('#', '@'), dtype=numpy.single)
    assert numpy.array_equal(rmsd_data[:, 0], rg_data[:, 0]), \
        "The times in RMSD file do not match those in radius of gyration file"
    data = numpy.core.records.fromarrays(numpy.hstack([rmsd_data, rg_data[:, 1:]]).T,
        names='time, rmsd, rg, rg_x, rg_y, rg_z')
    dset = group.create_dataset(dataset, data=data)
    dset.attrs['time_unit'] = time_unit
    dset.attrs['unit'] = unit


def add_rmsf(hdf_file, rmsf, begin, end, unit='nanometer'):
    group = hdf_file.require_group('rmsf')
    if isinstance(rmsf, list):
        rmsf_split_by_chains = []
        for xvg_file in rmsf:
            rmsf_data = numpy.loadtxt(xvg_file, comments=('#', '@'), dtype=numpy.single)
            rmsf_split_by_chains.append(rmsf_data)
    else:
        rmsf_data = numpy.loadtxt(rmsf, comments=('#', '@'), dtype=numpy.single)
        rmsf_split_by_chains = split_increasing(rmsf_data.T)
    for ch, chain_rmsf in enumerate(rmsf_split_by_chains):
        group.create_dataset(f'chain {ch}', data=chain_rmsf)
    group.attrs['unit'] = unit
    group.attrs['begin'] = f'{begin} ns'
    group.attrs['end'] = f'{end} ns'


def parse_dssp(filename):
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
    """ Evaluate secondary structure counts per residue """
    ss_dtype = [(_, numpy.uint32) for _ in SECSTR_CODES.keys()]
    output = numpy.zeros(len(chain_dssp_data.T), dtype=ss_dtype)
    for resSeq, residue in enumerate(chain_dssp_data.T):
        counts = numpy.unique(residue, return_counts=True)
        for structure, count in zip(*counts):
            output[structure.decode('ascii')][resSeq] = count
    return output


def add_dssp(hdf_file, dssp):
    """ Add the secondary structure data exhibited by each residue in
        the output of gmx do_dssp
    """
    dssp_group = hdf_file.require_group('secondary_structure/dssp_data')
    ss_counts_group = hdf_file.require_group('secondary_structure/counts_per_residue')
    dssp_data = parse_dssp(dssp)
    for ch, chain_dssp_data in enumerate(dssp_data):
        dssp_group.create_dataset(f'chain {ch}', data=chain_dssp_data)
        ss_counts_group.create_dataset(f'chain {ch}',
            data=ss_count_per_residue(chain_dssp_data)
        )


def get_dihedrals(hdf_file, trajectory, structure, chunk=1000, atoms=None):
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

    group = hdf_file.require_group('dihedrals')
    phi_dset = group.create_dataset( 'phi', data=numpy.vstack( phi_array ) )
    psi_dset = group.create_dataset( 'psi', data=numpy.vstack( psi_array ) )
    omega_dset = group.create_dataset( 'omega', data=numpy.vstack( omega_array ) )

    phi_dset.attrs['indices'] = phi_indices.astype(numpy.int32)
    psi_dset.attrs['indices'] = psi_indices.astype(numpy.int32)
    omega_dset.attrs['indices'] = omega_indices.astype(numpy.int32)
    return numpy.hstack( time )


def get_dihedral_sd(hdf_file, begin=0, end=1000, step=200):
    """
        hdf_file must contain dihedrals and sequence

        begin = Start frame in ns
        end = End frame in ns
        step = Number of frames per ns
    """
    start = int(begin * step)
    end = int(end * step)

    dih_sd_group = hdf_file.require_group('dihedral_standard_deviation')
    dih_sd_group.attrs['unit'] = 'radians'
    dih_sd_group.attrs['begin'] = f'{begin} ns'
    dih_sd_group.attrs['end'] = f'{end} ns'

    lengths = [len(_) - 1  for _ in hdf_file.attrs['sequence'].split('/')]
    data_type =  numpy.dtype([("phi", numpy.float),
                              ("psi", numpy.float),
                              ("omega", numpy.float)
                            ])

    for ch, length in enumerate(lengths):
        dih_sd_group.create_dataset(f'chain {ch}', (length + 1,), dtype=data_type)

    for angle in ['phi', 'psi', 'omega']:
        _, num_angles = numpy.shape(hdf_file[f'/dihedrals/{angle}'])
        assert sum(lengths) == num_angles
        circ_sd = scipy.stats.circstd(hdf_file[f'/dihedrals/{angle}'][start:end], axis=0, high=numpy.pi)
        for ch, chain_dih_sd in enumerate(split_chains(circ_sd, lengths)):
            dset = hdf_file[f'/dihedral_standard_deviation/chain {ch}']
            if angle == 'phi':
                dset[f'{angle}', 1:] = chain_dih_sd
            else:
                dset[f'{angle}', :-1] = chain_dih_sd


def add_dipoles(hdf_file, dipoles):
    dipole_moment = numpy.loadtxt(dipoles, comments=('#', '@'), dtype=numpy.single)
    data = numpy.core.records.fromarrays(dipole_moment[:, 1:].T, names='mu_x, mu_y, mu_z, mu')
    dset = hdf_file.create_dataset('dipole_moment', data=data)
    dset.attrs['unit'] = 'Debye'


def main(argv):
    if argv:
        args = docopt.docopt(__doc__, argv=argv)
    else:
        args = docopt.docopt(__doc__)

    if not os.path.exists(args['HDF_FILE']):
        print(f"The file {args['HDF_FILE']} does not exist and will be created.")

    with h5py.File(args['HDF_FILE'], 'a') as hdf_file:
        if args['--info']:
            add_info(hdf_file=hdf_file, info=args['--info'])
        if args['--all_atom_rmsd'] and args['--all_atom_rg']:
            add_rmsd_rg(hdf_file=hdf_file,
                        rmsd=args['--all_atom_rmsd'],
                        rg=args['--all_atom_rg'],
                        dataset='all_atom',
            )
        if args['--backbone_rmsd'] and args['--backbone_rg']:
            add_rmsd_rg(hdf_file=hdf_file,
                        rmsd=args['--backbone_rmsd'],
                        rg=args['--backbone_rg'],
                        dataset='backbone',
            )
        if args['--c_alpha_rmsd'] and args['--c_alpha_rg']:
            add_rmsd_rg(hdf_file=hdf_file,
                        rmsd=args['--c_alpha_rmsd'],
                        rg=args['--c_alpha_rg'],
                        dataset='c-alpha',
            )
        if args['--rmsf']:
            add_rmsf(hdf_file=hdf_file, rmsf=args['--rmsf'].strip().split(),
                        begin=args['<begin>'], end=args['<end>'])

        if args['--ss']:
            add_dssp(hdf_file=hdf_file, dssp=args['--ss'])

        if args['--dipoles']:
            add_dipoles(hdf_file=hdf_file, dipoles=args['--dipoles'])

        if args['--atoms']:
            atom_list = eval(args['--atoms'])
        else:
            atom_list = None

        if args['--trajectory'] and args['--structure'] and args['--chunk']:
            time = get_dihedrals(hdf_file=hdf_file,
                                 trajectory=args['--trajectory'],
                                 structure=args['--structure'],
                                 chunk=int(args['--chunk']),
                                 atoms=atom_list,
            )
            hdf_file.create_dataset('time', data=time)
            hdf_file.attrs['time_unit'] = 'picosecond'
            get_dihedral_sd(hdf_file=hdf_file)


if __name__ == '__main__':
    main()
