""" Script for structural alignment of two PDB files
    Formerly named: biopython_align.py """
import Bio.PDB
import argparse
import numpy
import decimal


def get_CA_atoms(model):
    """ Return a list of C-alpha atoms (in the structures) you wish to align. """
    CA_atoms = []
    for chain in model:
        for residue in chain:
            try:
                CA_atoms.append(residue['CA'])
            except KeyError:
                print('KeyError: for ', residue)
    return CA_atoms


def structure_alignment(sample, reference, output=None,
    sample_model=0, reference_model=0):
    """ Superimpose the sample structure onto the reference structure
        using only the C-alpha atoms.
        
        Previously: superimpose(reference, sample)
    """
    pdb_parser = Bio.PDB.PDBParser()
    ref_structure = pdb_parser.get_structure("reference", reference)
    sample_structure = pdb_parser.get_structure("sample", sample)

    ref_model = ref_structure[reference_model]
    sample_model = sample_structure[sample_model]

    ref_atoms = get_CA_atoms(ref_model)
    sample_atoms = get_CA_atoms(sample_model)

    # Now we initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())

    # Save the aligned version of sample_structure
    if output:
        io = Bio.PDB.PDBIO()
        io.set_structure(sample_structure)
        io.save(output)

    # Return Roation + Translation Matrin & RMSD:
    return super_imposer.rotran


def main():
    parser = argparse.ArgumentParser(description='Structure Alignment')
    parser.add_argument('-r', '--reference', required=True, help='input file', metavar='<PDB File>')
    parser.add_argument('-s', '--sample', required=True, help='input file', metavar='<PDB File>')
    parser.add_argument('-o', '--output', help='output filename', metavar='<PDB File>')
    args = parser.parse_args()
    rotation, translation = structure_alignment(sample=args.sample,
                                  reference=args.reference,
                                  output=args.output,
    )
    print('Rotation Matrix')
    print(rotation)
    print('Translation')
    print(translation)


if __name__ == '__main__':
    main()
