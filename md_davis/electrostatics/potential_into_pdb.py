import pandas
import Bio.PDB
import warnings
import collections
import string
import pickle
import argparse
from md_davis.collect_data import collect_residue_wise_data


def parse_atomic_potential(potential_file):
    df = pandas.read_csv(potential_file, skiprows=12, delim_whitespace=True, skipfooter=2,
        dtype={'resSeq': int}, engine='python',
        names=['name', 'resName', 'chainID', 'resSeq', 'potential',
            'reaction', 'coulomb', 'Ex', 'Ey', 'Ez'],
    )
    output = {}
    chain = 0
    for _, data in df.groupby(['chainID'], as_index=False):
        grouped_df = data.groupby(['resSeq', 'resName', 'name'], as_index=False)['potential']
        potential =  grouped_df.sum()
        potential.rename(columns={'potential':'total'}, inplace=True)
        potential['mean'] =  grouped_df.mean()['potential']
        output[f'chain {chain}'] = potential
        chain += 1
    return output


def parse_surface(surface_file):
    vertex = collections.namedtuple('vertex', 'chain name resName resSeq')
    surface = []
    with open(surface_file, 'r') as surf_file:
        lines = surf_file.readlines()
        for line in lines:
            if line[:4] == 'ATOM':
                name = line[12:16]
                resName = line[17:20]
                chain = line[21]
                resSeq = int(line[22:26])
                vertex(name=name,
                       resName=resName,
                       chain=chain,
                       resSeq=resSeq,
                )
                surface.append(vertex)
                # print(f'ATOM        {name:^4} {resName:3} {chain:1}{resSeq:4}')
                # print(line)
    return surface


def add_atomic_potentials(potential, pdb_file, surface, output=None, directory=None, model=0):
    """ Add electrostatics surface potential as B-factor and occupancy in a PDB
    """
    potential_df = parse_atomic_potential(potential)
    vertices = parse_surface(surface)
    pdb_parser = Bio.PDB.PDBParser()
    structure = pdb_parser.get_structure("pdb_file", pdb_file)[model]

    chain_dict = {letter: f'chain {_}' for _, letter in enumerate(string.ascii_uppercase)}

    for residue in structure.get_residues():
        resid = residue.id[1]
        chain = chain_dict[residue.parent._id]
        residue_potentials = potential_df[chain].loc[potential_df[chain]['resSeq'] == resid]
        atoms = list(residue_potentials['name'])
        for atom in residue.get_atoms():
            row = residue_potentials.loc[residue_potentials['name'] == atom.get_fullname().strip()]
            if atom.name in atoms:
                atoms.remove(atom.name)
                atom.set_bfactor(row['total'])
                atom.set_occupancy(row['mean'])
        if len(atoms) != 0:
            warnings.warn(f'All atoms in surface not exhausted')
            print(atoms)
            # for _ in residue.get_atoms():
            #     print(vars(_))
            # print(residue)

    # Save the aligned version of sample_structure
    if output:
        io = Bio.PDB.PDBIO()
        io.set_structure(structure)
        io.save(output)


def main(potential, pdb_file, output=None, pickle_file=None, model=0):
    """ Add electrostatics surface potential as B-factor and occupancy in a PDB
    """
    potential_df = collect_residue_wise_data.parse_potential(potential)
    pdb_parser = Bio.PDB.PDBParser()
    structure = pdb_parser.get_structure("pdb_file", pdb_file)[model]

    sub_df = pickle.load(open(pickle_file, 'rb'))['data']['chain 0']
    std = sub_df['surface_potential']['std_mean']
    resID = sub_df['sequence']['resi']
    std_dict = {i: s for i, s in zip(resID, std)}
    chain_dict = {letter: f'chain {_}' for _, letter in enumerate(string.ascii_uppercase)}

    for residue in structure.get_residues():
        resid = residue.id[1]
        # chain = chain_dict[residue.parent._id]
        chain = 'chain 0'
        residue_potentials = potential_df[chain].loc[potential_df[chain]['resSeq'] == resid, 'mean']
        if not residue_potentials.empty:
            for atom in residue.get_atoms():
                atom.set_bfactor(std_dict[resid])
                atom.set_occupancy(residue_potentials.astype(float))
        else:
            for atom in residue.get_atoms():
                atom.set_bfactor(0)
                atom.set_occupancy(0)
    # Save the aligned version of sample_structure
    if output:
        io = Bio.PDB.PDBIO()
        io.set_structure(structure)
        io.save(output)


if __name__ == "__main__":
    for i in range(0, 1010, 10):
        main(potential=f'/run/media/djmaity/MASS1/dibyajyoti/processed_trajectories/acylphosphatase_output/2GV1_output/uniform_sample_surface_electrostatics/2GV1_sample_{i}.pot',
            pdb_file=f'/run/media/djmaity/MASS1/dibyajyoti/processed_trajectories/acylphosphatase_output/2GV1_output/uniform_sample/2GV1_sample_{i}.pdb',
            output=f'2GV1_sample_{i}.pdb',
            pickle_file='/home/djmaity/Lab_Home/current/collect_residue_wise_data/acylphosphatase/2GV1_residue_wise_data.p',
        )
