import mdtraj
import pandas

from md_davis.electrostatics.electrostatics import parse_electrostatic_potential


def save_pdb_with_potentials(df, trajectory, output_filename, save_mean_potentials=True, fillna=False, normalize=False):
    coordinates = pandas.DataFrame(10 * trajectory.xyz[0],
                                   columns=['x', 'y', 'z'])  # Multiply by 10 to convert nm to Angstrom
    df = pandas.concat([df, coordinates], axis=1)

    df['name'] = df['name'].str.replace(' ', '')
    df.loc[df['name'].str.len() < 4, 'name'] = ' ' + df.loc[df['name'].str.len() < 4, 'name']

    df['element'] = df['element'].str.strip()
    df['resName'] = df['resName'].str.replace(' ', '')
    df['segmentID'] = df['segmentID'].str.replace(' ', '')

    df['chainID'] = df['chainID'].mod(26).apply(lambda x: chr(x + 65))

    df['serial'] = df['serial'].mod(100000)
    df['resSeq'] = df['resSeq'].mod(10000)

    if save_mean_potentials:
        potential_col = 'mean'
    else:
        potential_col = 'total'

    if fillna:
        df[potential_col] = df[potential_col].fillna(0.0)

    if normalize:
        max_potential = df[potential_col].abs().max()
        df[potential_col] = 9.99 * df[potential_col].fillna(0.0) / max_potential

    width = df[potential_col].round(decimals=2).astype(str).str.len().max()
    df[potential_col] = df[potential_col].round(decimals=2).map('{:.2f}'.format).astype(str).str.rjust(width)

    with open(output_filename, 'w') as output_pdb_file:
        for index, row in df.iterrows():
            print(f"ATOM  {row['serial']:>5} {row['name']:4} {row['resName']:3} {row['chainID']:1}{row['resSeq']:>4}"
                  "    {row['x']:8.3f}{row['y']:8.3f}{row['z']:8.3f}  1.00 {row[potential_col]} {row['element']}",
                  file=output_pdb_file)


def potential_into_pdb(potential_file, pdb_file, output_filename, atomic_potentials=False, save_mean_potentials=True,
                       normalize=False, fillna=False):
    potential_data = parse_electrostatic_potential(potential_file, atomic_potentials=atomic_potentials)

    potential_dataframes = []
    for chain, df in potential_data.items():
        df['chainID'] = int(chain[6:])
        potential_dataframes.append(df)
    df_potentials = pandas.concat(potential_dataframes, ignore_index=True)

    trajectory = mdtraj.load_pdb(pdb_file, standard_names=False)
    structure = trajectory.topology

    df_pdb = structure.to_dataframe()[0]

    if atomic_potentials:
        df_pdb = df_pdb.merge(df_potentials, how='left', on=['chainID', 'resSeq', 'name'])
    else:
        df_pdb = df_pdb.merge(df_potentials, how='left', on=['chainID', 'resSeq'])

    save_pdb_with_potentials(df=df_pdb,
                             trajectory=trajectory,
                             output_filename=output_filename,
                             save_mean_potentials=save_mean_potentials,
                             normalize=normalize,
                             fillna=fillna)
