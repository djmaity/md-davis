import argparse
import h5py
import sys
sys.path.append("C:/Users/djmaity/Documents/Repositories")
sys.path.append("/home/djmaity/Documents/Repositories")
from analyze_md.structure import get_pdb_sequence

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-f', '--file', required=True,
                    help='HDF5 file')
parser.add_argument('-p', '--pdb', required=True,
                    help='PDB file')
args = parser.parse_args()
with h5py.File(args.file, 'a') as hdf5_file:
    hdf5_file.attrs['sequence'] = get_pdb_sequence.get_fasta(args.pdb)
