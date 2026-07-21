"""Validate PDB residue numbering and backbone completeness."""

from __future__ import annotations

import argparse

from md_davis.core.amino_acids import THREE_TO_ONE

BACKBONE_ATOMS = (" N  ", " CA ", " C  ", " O  ")


def check_pdb(input_file: str) -> str:
    """Return a one-letter sequence while reporting gaps and missing atoms."""
    sequence: list[str] = []
    issues: list[str] = []
    previous_residue: int | None = None
    atoms_present = dict.fromkeys(BACKBONE_ATOMS, False)

    with open(input_file, encoding="utf-8") as pdb_file:
        for line in pdb_file:
            record = line[:6]
            if record == "TER   ":
                if previous_residue is not None and not all(atoms_present.values()):
                    issues.append(f"Residue {previous_residue} missing backbone atoms")
                break
            if record != "ATOM  ":
                continue

            residue_number = int(line[22:26])
            if residue_number != previous_residue:
                if previous_residue is not None:
                    if not all(atoms_present.values()):
                        issues.append(f"Residue {previous_residue} missing backbone atoms")
                    if residue_number < previous_residue:
                        issues.append(f"Residue numbering decreases from {previous_residue} to {residue_number}")
                        break
                    if residue_number > previous_residue + 1:
                        issues.append(f"Missing residues from {previous_residue + 1} to {residue_number - 1}")
                        sequence.extend("X" * (residue_number - previous_residue - 1))

                residue_name = line[17:20]
                sequence.append(THREE_TO_ONE.get(residue_name, "X"))
                atoms_present = dict.fromkeys(BACKBONE_ATOMS, False)
                previous_residue = residue_number

            atom_name = line[12:16]
            if atom_name in atoms_present:
                atoms_present[atom_name] = True

    if issues:
        print("\n".join(issues))
    return "".join(sequence)


def main() -> None:
    parser = argparse.ArgumentParser(description="Check a PDB file and print its FASTA sequence")
    parser.add_argument("pdbfile", help="PDB structure file", metavar="FILE.pdb")
    args = parser.parse_args()

    sequence = check_pdb(args.pdbfile)
    print(f">{args.pdbfile[:4]}\n{sequence}")


if __name__ == "__main__":
    main()
