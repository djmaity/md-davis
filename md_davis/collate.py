#! /usr/bin/env python
"""

"""

import os
import collections
import pandas
import toml
import h5py
import click
import mdtraj
import numpy
import scipy.stats
import warnings

import md_davis


SECSTR_CODES = {
    'H': 'α-helix',
    'G': '3_10-helix',
    'I': 'π-helix',
    'B': 'β-bridge',
    'E': 'β strand',
    'T': 'Turn',
    'S': 'Bend',
    '~': 'Loop',
}


def split_increasing(array, use_col=0):
    """ Split an array when it is not increasing """
    array = numpy.array(array)
    if array.ndim == 1:
        column = array
    elif array.ndim == 2:
        column = array[use_col]
    else:
        raise RuntimeError(
            "Can not handle data with dimension other than 1 or 2")
    # Check if following value is greater
    where_greater = numpy.greater(column[:-1], column[1:])
    # Use the indices where the above is True
    split_indices = numpy.flatnonzero(where_greater)
    # to split the array
    split_arrays = numpy.split(array, split_indices + 1, axis=array.ndim - 1)
    return list(map(numpy.transpose, split_arrays))


def split_chains(array, lengths):
    """ """
    if sum(lengths) != len(array):
        warnings.warn('Number of dihedral angles is inconsistent '
                      'with the number of residues.')
    start = 0
    for end in lengths:
        end = start + end
        sub_array = array[start:end]
        start = end
        yield sub_array


def add_rmsd_rg(hdf_file, rmsd, rg, group_name, dataset,
                unit=None, time_unit=None, comment=None):
    """ Make RMSD & Rg dataset in HDF5 file """
    group = hdf_file.require_group(group_name)
    if comment:
        group.attrs['comment'] = comment
    print('Collecting data into HDF5 file:\n'
          f'RMSD from {rmsd}\n'
          f'Radius of gyration from {rg}\n')
    rmsd_data = numpy.loadtxt(rmsd, comments=('#', '@'))
    rg_data = numpy.loadtxt(rg, comments=('#', '@'))
    if not numpy.array_equal(rmsd_data[:, 0], rg_data[:, 0]):
        print('The times in RMSD file do not match '
              'those in radius of gyration file')
    data = numpy.core.records.fromarrays(
        numpy.hstack([rmsd_data, rg_data[:, 1:]]).T,
        names='time, rmsd, rg, rg_x, rg_y, rg_z'
    )
    dset = group.require_dataset(name=dataset,
                                 shape=data.shape,
                                 dtype=data.dtype,
                                 data=data)
    dset.attrs['unit'] = 'nm' if unit is None else unit
    dset.attrs['time_unit'] = 'ps' if time_unit is None else time_unit


def add_rmsf(hdf_file, rmsf, **attributes):
    group = hdf_file.require_group('residue_property/rmsf')
    if isinstance(rmsf, list) and len(rmsf) > 1:
        rmsf_split_by_chains = []
        for xvg_file in rmsf:
            rmsf_data = numpy.loadtxt(
                xvg_file, comments=('#', '@'), dtype=numpy.single)
            rmsf_split_by_chains.append(rmsf_data)
    elif isinstance(rmsf, list) and len(rmsf) == 1:
        rmsf_data = numpy.loadtxt(
            rmsf[0], comments=('#', '@'), dtype=numpy.single)
        rmsf_split_by_chains = split_increasing(rmsf_data.T)
    else:
        rmsf_data = numpy.loadtxt(rmsf, comments=('#', '@'),
                                  dtype=numpy.single)
        rmsf_split_by_chains = split_increasing(rmsf_data.T)
    for ch, chain_rmsf in enumerate(rmsf_split_by_chains):
        group.require_dataset(f'chain {ch}',
                              shape=chain_rmsf.shape,
                              dtype=chain_rmsf.dtype,
                              data=chain_rmsf)
    for key, value in attributes.items():
        group.attrs['key'] = value


def parse_secondary_structure(filename):
    """ Parse the DSSP data file """
    secondary_structure = []
    with open(filename, 'r') as dat_file:
        next(dat_file)
        for line in dat_file:
            secondary_structure.append(line.strip().split('='))
    output = []
    for chain in zip(*secondary_structure):
        assert len(set([len(_) for _ in chain])) <= 1, \
            "The length of the secondary structure strings are not equal. " \
            "Please check your input .dat file."
        output.append(numpy.array([list(_) for _ in chain], dtype=f'S1'))
    return output


def ss_count_per_residue(chain_ss_data):
    """ Evaluate secondary structure counts per residue """
    ss_dtype = [(_, numpy.uint32) for _ in SECSTR_CODES.keys()]
    output = numpy.zeros(len(chain_ss_data.T), dtype=ss_dtype)
    for resSeq, residue in enumerate(chain_ss_data.T):
        counts = numpy.unique(residue, return_counts=True)
        for structure, count in zip(*counts):
            output[structure.decode('ascii')][resSeq] = count
    return output


def add_secondary_structure(hdf_file, secondary_structure):
    """ Add the secondary structure data exhibited by each residue in
        the output of gmx do_secondary_structure
    """
    ss_group = hdf_file.require_group(
        'secondary_structure/secondary_structure_data')
    ss_counts_group = hdf_file.require_group(
        'secondary_structure/counts_per_residue')
    ss_data = parse_secondary_structure(secondary_structure)
    for ch, chain_ss_data in enumerate(ss_data):
        ss_group.require_dataset(
            f'chain {ch}',
            data=chain_ss_data,
            shape=chain_ss_data.shape,
            dtype=chain_ss_data.dtype,
        )
        counts = ss_count_per_residue(chain_ss_data)
        ss_counts_group.require_dataset(
            f'chain {ch}',
            data=counts,
            shape=counts.shape,
            dtype=counts.dtype,
        )


def get_dihedrals(hdf_file, trajectory, structure,
                  chunk=1000, atoms=None, start=0, stride=1):
    """ Evaluate the dihedral angles and store it as HDF5 dataset """
    trj_iterator = mdtraj.iterload(
        trajectory, top=structure, chunk=chunk, atom_indices=atoms, skip=start,
        stride=stride)
    first_trj_chunk = next(trj_iterator)
    print(f'Calculating dihedrals for trajectory between '
          f'{first_trj_chunk.time[0]} and {first_trj_chunk.time[-1]} ps')
    phi_indices, phi = mdtraj.compute_phi(first_trj_chunk)
    psi_indices, psi = mdtraj.compute_psi(first_trj_chunk)
    omega_indices, omega = mdtraj.compute_omega(first_trj_chunk)

    time = [first_trj_chunk.time]
    phi_array = [phi]
    psi_array = [psi]
    omega_array = [omega]
    for trj_chunk in trj_iterator:
        print(f'Calculating dihedrals for trajectory between '
              f'{trj_chunk.time[0]} and {trj_chunk.time[-1]} ps')
        time.append(trj_chunk.time)
        phi_array.append(mdtraj.compute_phi(trj_chunk)[1])
        psi_array.append(mdtraj.compute_psi(trj_chunk)[1])
        omega_array.append(mdtraj.compute_omega(trj_chunk)[1])

    group = hdf_file.require_group('dihedrals')
    phi_data = numpy.vstack(phi_array)
    phi_dset = group.require_dataset('phi',
                                     data=phi_data,
                                     shape=phi_data.shape,
                                     dtype=phi_data.dtype)
    psi_data = numpy.vstack(psi_array)
    psi_dset = group.require_dataset('psi',
                                     data=psi_data,
                                     shape=psi_data.shape,
                                     dtype=psi_data.dtype)
    omega_data = numpy.vstack(omega_array)
    omega_dset = group.require_dataset('omega',
                                       data=omega_data,
                                       shape=omega_data.shape,
                                       dtype=omega_data.dtype)
    phi_dset.attrs['indices'] = phi_indices.astype(numpy.int32)
    psi_dset.attrs['indices'] = psi_indices.astype(numpy.int32)
    omega_dset.attrs['indices'] = omega_indices.astype(numpy.int32)
    return numpy.hstack(time)


def get_dihedral_sd(hdf_file, chain_lengths, start=0, end=1000, step=200):
    """
        hdf_file must contain dihedrals and sequence

        start = Start frame in ns
        end = End frame in ns
        step = Number of frames per ns
    """
    chain_lengths = [_ - 1 for _ in chain_lengths]

    start = int(start * step)
    end = int(end * step)

    dih_sd_group = hdf_file.require_group('dihedral_standard_deviation')
    dih_sd_group.attrs['unit'] = 'radians'
    dih_sd_group.attrs['start'] = f'{start} ns'
    dih_sd_group.attrs['end'] = f'{end} ns'

    data_type = numpy.dtype([("phi", numpy.float),
                             ("psi", numpy.float),
                             ("omega", numpy.float)
                             ])
    for ch, length in enumerate(chain_lengths):
        dih_sd_group.require_dataset(
            f'chain {ch}', shape=(length + 1,), dtype=data_type)

    for angle in ['phi', 'psi', 'omega']:
        _, num_angles = numpy.shape(hdf_file[f'/dihedrals/{angle}'])

        if sum(chain_lengths) != num_angles:
            warnings.warn(f'Number of {angle} angles is inconsistent '
                          'with the number of residues.')
        circ_sd = scipy.stats.circstd(
            hdf_file[f'/dihedrals/{angle}'][start:end], axis=0, high=numpy.pi)
        for ch, dih_sd in enumerate(split_chains(circ_sd, chain_lengths)):
            dset = hdf_file[f'/dihedral_standard_deviation/chain {ch}']
            if angle == 'phi':
                dset[f'{angle}', 1:] = dih_sd
            else:
                dset[f'{angle}', :-1] = dih_sd


def add_dipoles(hdf_file, dipoles):
    dipole_moment = numpy.loadtxt(
        dipoles, comments=('#', '@'), dtype=numpy.single)
    data = numpy.core.records.fromarrays(
        dipole_moment[:, 1:].T, names='mu_x, mu_y, mu_z, mu')
    dset = hdf_file.require_dataset('dipole_moment',
                                    data=data,
                                    shape=data.shape,
                                    dtype=data.dtype)
    dset.attrs['unit'] = 'Debye'


def add_sasa(hdf_file, sasa, unit='nm^2'):
    group = hdf_file.require_group('residue_property/sasa')
    if isinstance(sasa, list) and len(sasa) > 1:
        sasa_split_by_chains = []
        for xvg_file in sasa:
            sasa_data = numpy.loadtxt(
                xvg_file, comments=('#', '@'), dtype=numpy.single)
            sasa_split_by_chains.append(sasa_data)
    elif isinstance(sasa, list) and len(sasa) == 1:
        sasa_data = numpy.loadtxt(
            sasa[0], comments=('#', '@'), dtype=numpy.single)
        sasa_split_by_chains = split_increasing(sasa_data.T)
    else:
        sasa_data = numpy.loadtxt(
            sasa, comments=('#', '@'), dtype=numpy.single)
        sasa_split_by_chains = split_increasing(sasa_data.T)
    for ch, chain_sasa in enumerate(sasa_split_by_chains):
        data = numpy.core.records.fromarrays(
            chain_sasa.T,
            names='time, average, standard_deviation')
        group.require_dataset(f'chain {ch}',
                              data=data,
                              shape=data.shape,
                              dtype=data.dtype)
    group.attrs['unit'] = unit


def add_surface_potential(hdf_file, surface_potential):
    """ Add potential """

    # TODO if single file directly place into dataframe else combine multiple potential files
    group = hdf_file.require_group('residue_property/surface_potential')

    total_df = {}
    mean_df = {}
    index = 1
    for file in os.listdir(surface_potential):
        if file.endswith('.pot'):
            print('Parsing: ', file)
            potential_file = os.path.join(surface_potential, file)
            parsed_potentials = md_davis.electrostatics.electrostatics.parse_electrostatic_potential(potential_file)
            for chain, potential_df in parsed_potentials.items():
                if chain not in total_df:
                    total_df[chain] = potential_df[['resSeq', 'total']]
                else:
                    total_df[chain] = total_df[chain].merge(potential_df[['resSeq', 'total']], how='outer', left_on='resSeq', right_on='resSeq', suffixes=(None, str(index)))
                if chain not in mean_df:
                    mean_df[chain] = potential_df[['resSeq', 'mean']]
                else:
                    mean_df[chain] = mean_df[chain].merge(potential_df[['resSeq', 'mean']], how='outer', left_on='resSeq', right_on='resSeq', suffixes=(None, str(index)))
            index += 1

    for chain in total_df.keys():
        tot_df = total_df[chain].set_index('resSeq')
        avg_df = mean_df[chain].set_index('resSeq')

        surface_potential_df = pandas.DataFrame()
        surface_potential_df['mean_of_total'] = tot_df.mean(skipna=True, axis=1)
        surface_potential_df['std_of_total'] = tot_df.std(skipna=True, ddof=0, axis=1)
        surface_potential_df['mean_of_mean'] = avg_df.mean(skipna=True, axis=1)
        surface_potential_df['std_of_mean'] = avg_df.std(skipna=True, ddof=0, axis=1)
        surface_potential_df = surface_potential_df.to_records()

        group.require_dataset(chain,
                              shape=surface_potential_df.shape,
                              dtype=surface_potential_df.dtype,
                              data=surface_potential_df)



CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(name='collate', context_settings=CONTEXT_SETTINGS)
@click.argument('input_files', nargs=-1, required=True,
                type=click.Path(exists=True))
def main(input_files):
    """Input TOML files."""
    for toml_file in input_files:
        data = toml.load(toml_file)

        stop = False
        if 'output' not in data:
            print(f"Please specify 'output' in {toml_file}")
            stop = True
        if 'name' not in data:
            print(f"Please specify 'name' in {toml_file}")
            stop = True

        if stop:
            return

        output = data['output']
        if os.path.exists(output):
            print(f"The existing file {output} will be amended.")
        else:
            print(f"The file {output} does not exist and will be created.")

        with h5py.File(output, 'a') as hdf_file:
            hdf_file.attrs['name'] = data['name']

            if 'sequence' in data:
                for ch, seq in data['sequence'].items():
                    resi, resn = zip(*seq)
                    group = hdf_file.require_group(f'sequence/chain {ch[6:]}')  # omit 'Chain '
                    resi = numpy.array(resi, dtype=int)
                    resn = numpy.array(resn, dtype="S3")
                    group.require_dataset('resi', shape=resi.shape, dtype=resi.dtype, data=resi)
                    group.require_dataset('resn', shape=resn.shape, dtype=resn.dtype, data=resn)
            elif 'structure' in data:
                sequences = md_davis.sequence.get_sequence(data['structure'])
                for ch, seq in sequences.items():
                    group = hdf_file.require_group(f'sequence/chain {ch}')
                    resi = numpy.array(seq[0], dtype=int)
                    resn = numpy.array(seq[1], dtype="S3")
                    group.require_dataset('resi', shape=resi.shape, dtype=resi.dtype, data=resi)
                    group.require_dataset('resn', shape=resn.shape, dtype=resn.dtype, data=resn)
            else:
                print("Please provide 'sequence' or 'structure' in the input toml file to parse the amino acid sequence")
                return


            # Add RMSD and radius of gyration into the HDF5 file
            if 'timeseries' in data:
                if 'rmsd' in data['timeseries'] and 'rg' in data['timeseries']:
                    time_unit, unit, comment = None, None, None
                    if 'time_unit' in data['timeseries']:
                        time_unit = data['timeseries']['time_unit']
                    if 'unit' in data['timeseries']:
                        unit = data['timeseries']['unit']
                    if 'comment' in data['timeseries']:
                        comment = data['timeseries']['comment']

                    add_rmsd_rg(hdf_file=hdf_file,
                                rmsd=data['timeseries']['rmsd'],
                                rg=data['timeseries']['rg'],
                                group_name='timeseries',
                                dataset='rmsd_rg',
                                time_unit=time_unit,
                                unit=unit,
                                comment=comment,
                                )
                elif 'rmsd' not in data['timeseries']:
                    print(f"Please provide 'rmsd' under "
                          f"[timeseries] in {toml_file}")
                elif 'rg' not in data['timeseries']:
                    print(f"Please provide 'rg' under "
                          "[timeseries] in {toml_file}")
                else:
                    print(f"Please provide 'rmsd' and 'rg' under "
                          f"[timeseries] in {toml_file}")

            if 'residue_property' in data:
                residue = data['residue_property']

                # Add RMSF into the residue_property group of the HDF5 file
                if 'rmsf' in residue:
                    rmsf = residue['rmsf']

                    # TODO: allow arbitrary key words pairs through toml files
                    start, end = None, None
                    if 'start' in rmsf:
                        start = rmsf['start']
                    if 'end' in rmsf:
                        end = rmsf['end']

                    add_rmsf(
                        hdf_file=hdf_file,
                        rmsf=rmsf['rmsf_files'],
                        start=start,
                        end=end,
                    )

                if 'secondary_structure' in residue:
                    add_secondary_structure(
                        hdf_file=hdf_file,
                        secondary_structure=residue['secondary_structure'])

                if 'sasa' in data['residue_property']:
                    add_sasa(hdf_file=hdf_file, sasa=residue['sasa'])
                #
                # if 'dipoles' in data['residue_property']:
                #     add_dipoles(
                #     hdf_file=hdf_file,
                #     dipoles=data['residue_property']['--dipoles'])

                if 'surface_potential' in data['residue_property']:
                    add_surface_potential(hdf_file=hdf_file, surface_potential=residue['surface_potential'])

            # Evaluate and add dihedral angles into the HDF5 file
            if 'dihedral' in data:
                if 'atoms' in data['dihedral']:
                    atoms = data['dihedral']['atoms']
                else:
                    atoms = None
                if 'chunk' in data['dihedral']:
                    chunk = data['dihedral']['chunk']
                else:
                    chunk = None

                time = get_dihedrals(hdf_file=hdf_file,
                                     trajectory=data['trajectory'],
                                     structure=data['structure'],
                                     chunk=chunk,
                                     atoms=atoms)

                # hdf_file.require_dataset('time',
                #                          data=time,
                #                          shape=time.shape,
                #                          dtype=time.dtype)
                # hdf_file.attrs['time_unit'] = 'picosecond'
                if 'chain_lengths' in data['dihedral']:
                    ch_len = data['dihedral']['chain_lengths']
                else:
                    ch_len = [_.n_residues for _ in
                            mdtraj.load(data['structure']).topology.chains]

                get_dihedral_sd(hdf_file=hdf_file,
                                chain_lengths=ch_len)


if __name__ == '__main__':
    main()
