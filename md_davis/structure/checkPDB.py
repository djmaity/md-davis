import argparse

amino_acids = {"ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
               "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
               "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
               "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"}

resColumns = ((19, 22), (23, 26), (27, 30), (31, 34), (35, 38),
              (39, 42), (43, 46), (47, 50), (51, 54), (55, 58),
              (59, 62), (63, 66), (67, 70))


def check_dict(dict):
    out = True
    for key, value in dict.items():
        out = out and value
    return out


def checkPDB(inputFile):
    with open(inputFile) as pdbFile:
        old_res_num = 0
        out_string = ""
        seq = ""
        sequence = ""
        atoms_to_check = {' N  ': False, ' CA ': False,
                          ' C  ': False, ' O  ': False}
        FirstRes = False
        for line in pdbFile:
            if line[0:6] == "SEQRES":
                for i, j in resColumns:
                    if line[i:j] in amino_acids:
                        sequence = sequence + amino_acids[line[i:j]]
            if line[0:6] == "ATOM  ":
                res_num = int(line[22:26])
                res_diff = res_num - old_res_num
                atom = line[12:16]
                if res_diff > 0:
                    if FirstRes and not check_dict(atoms_to_check):
                        out_string = out_string + "Residue %d missing "
                        "backbone atoms\n" % old_res_num
                    FirstRes = True
                    for k in atoms_to_check:
                        atoms_to_check[k] = False
                    if res_diff == 1:
                        seq = seq + amino_acids[line[17:20]]
                    elif res_diff == 2:
                        out_string = out_string + \
                            "Missing Residue numer %d\n" % (old_res_num + 1)
                        seq = seq + "X"
                    else:
                        out_string = out_string + \
                            "Missing Residues from %d to %d\n" % (
                                old_res_num + 1, res_num - 1)
                        seq = seq + "X" * (res_diff - 1)
                if res_diff < 0:
                    print(old_res_num, res_num)
                    print("Wrong residue index")
                    break
                if atom in atoms_to_check:
                    atoms_to_check[atom] = True
                old_res_num = res_num
            if line[0:6] == "TER   ":
                if not check_dict(atoms_to_check):
                    out_string = out_string + "Residue %d missing backbone "
                    "atoms\n" % old_res_num
                break
    print(out_string)
    # print(sequence)
    return seq

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check PDB file and get fasta'
                                     ' sequence')
    parser.add_argument('pdbfile', help='PDB structure file',
                        metavar='FILE.pdb')
    args = parser.parse_args()

    seq_from_str = checkPDB(args.pdbfile)
    pdb = args.pdbfile[0:4]
    print(">", pdb, "\n", seq_from_str, "\n")

