"""
Create a pandas dataframe from the information in HDF5 file for plotting

Usage:  md_davis residue dataframe [options] HDF_FILE OUTPUT

Options:
  -h, --help                        Show this screen.
  -v, --version                     Show version.
  --prefix <string>                 Prefix used in the alignment file
  -a, --annotations <JSON_FILE>     JSON file containing annotations to mark on the plot
  --potentials <directory>          Dicrectory containing potential files
  --pdb_potentials <pot_file>       .pot file obtained from electrostatic calculation


"""

from __future__ import annotations

import collections
import pickle
from collections.abc import Mapping, MutableMapping, Sequence
from os import PathLike
from typing import Any

import h5py
import numpy
import pandas

# Import from my packages
from md_davis.analysis import secondary_structure as secStr_counts
from md_davis.core.amino_acids import ONE_TO_THREE
from md_davis.io.alignment import parse_clustal
from md_davis.io.config import load_toml

pandas.option_context("display.max_rows", None)


ResidueFrames = dict[str, pandas.DataFrame]
AlignedResidueFrames = MutableMapping[str, MutableMapping[str, pandas.DataFrame]]


def residue_dataframe(hdf_file: h5py.File | h5py.Group) -> ResidueFrames:
    """Create a pandas dataframe for Residue data to plot"""
    output_dict = {}
    for chain in hdf_file["sequence"]:
        array_to_concat = []
        keys_to_concat = []
        # Add sequence
        sequence = pandas.DataFrame(
            data=zip(hdf_file[f"sequence/{chain}/resi"], hdf_file[f"sequence/{chain}/resn"].asstr(), strict=False),
            columns=["resi", "resn"],
        )
        array_to_concat.append(sequence)
        keys_to_concat.append("sequence")
        # Add scondary structure to dtaframe from HDF5 file
        if f"secondary_structure/counts_per_residue/{chain}" in hdf_file:
            ss_counts = hdf_file[f"secondary_structure/counts_per_residue/{chain}"]
            secondary_structure_dict = {}
            for code in secStr_counts.SECSTR_CODES.keys():
                secondary_structure_dict[code] = ss_counts[code]
            secondary_structure_df = pandas.DataFrame(secondary_structure_dict)
            array_to_concat.append(secondary_structure_df)
            keys_to_concat.append("secondary_structure")
        # Add RMSF to dtaframe from HDF5 file
        if f"residue_property/rmsf/{chain}" in hdf_file:
            rmsf_df = pandas.DataFrame(hdf_file[f"residue_property/rmsf/{chain}"][:, 1])
            array_to_concat.append(rmsf_df)
            keys_to_concat.append("rmsf")
        # Add SASA to dtaframe from HDF5 file
        if f"residue_property/sasa/{chain}" in hdf_file:
            sasa_dict = {}
            for measure in ["average", "standard_deviation"]:
                sasa_dict[measure] = hdf_file[f"residue_property/sasa/{chain}"][measure]
            sasa_df = pandas.DataFrame(sasa_dict)
            array_to_concat.append(sasa_df)
            keys_to_concat.append("sasa")
        # Add dihedral standard deviation to dtaframe from HDF5 file
        if f"dihedral_standard_deviation/{chain}" in hdf_file:
            dih_sd_dict = {}
            for angle in ["phi", "psi", "omega"]:
                dih_sd_dict[angle] = hdf_file[f"dihedral_standard_deviation/{chain}"][angle]
            dih_sd_df = pandas.DataFrame(dih_sd_dict)
            array_to_concat.append(dih_sd_df)
            keys_to_concat.append("dihedral_standard_deviation")

        # Combine all of the above
        output_dict[chain] = pandas.concat(array_to_concat, keys=keys_to_concat, axis=1)
        if f"residue_property/surface_potential/{chain}" in hdf_file:
            potential = hdf_file[f"residue_property/surface_potential/{chain}"]
            pot_df = pandas.DataFrame.from_records(potential)
            columns = [
                ("surface_potential", "resSeq"),
                ("surface_potential", "mean_of_total"),
                ("surface_potential", "std_of_total"),
                ("surface_potential", "mean_of_mean"),
                ("surface_potential", "std_of_mean"),
            ]
            pot_df.columns = pandas.MultiIndex.from_tuples(columns)
            output_dict[chain] = output_dict[chain].merge(
                pot_df, how="left", left_on=[("sequence", "resi")], right_on=[("surface_potential", "resSeq")]
            )
            output_dict[chain].drop([("surface_potential", "resSeq")], axis="columns", inplace=True)
            # for _, pdb_row in pot_df[chain].iterrows():
            #     pot_df.loc[data_frame['sequence']['resi'] == pdb_row.resSeq, 'pdb_total'] = pdb_row['total']
            #     pot_df.loc[data_frame['sequence']['resi'] == pdb_row.resSeq, 'pdb_mean'] = pdb_row['mean']
            # columns += [
            #     ('surface_potential', 'pdb_total'),
            #     ('surface_potential', 'pdb_mean'),
            # ]

        output_dict[chain].set_index(("sequence", "resi"))
    return output_dict


def parse_alignment(alignment_file: str | PathLike[str]) -> dict[str, str]:
    """Parse an alignment file in CLUSTAL format."""
    return parse_clustal(alignment_file)


# def align_residues(residue_dataframes, alignment, mapping=None):
#     """ Align residue dataframes based on sequence alignment """
#
#     for chain, alignment_file in alignment['alignment'].items():
#         aln_dict = parse_alignment(alignment)
#         if aln_dict is None:
#             return
#
#         # Length of sequences in alignment file
#         aln_seq_length = len(next(iter(aln_dict.values())))
#
#         aligned_dataframes = collections.defaultdict(dict)
#         for name, data in residue_dataframes.items():
#             aligned_dataframes[name] = data
#
#
#
#             out_df = pandas.DataFrame(
#                 data=numpy.zeros((aln_seq_length, df.shape[1])),
#                 columns=df.columns)
#             out_df['sequence', 'resn'] = '-'
#             i, j = 0, 0
#             for residue in aln_dict[name]:
#                 if residue != '-':
#                     if residue == df.iloc[j].sequence.resn:
#                         out_df.iloc[i] = df.iloc[j]
#                         j += 1
#                     else:
#                         print('Mismatch between sequence in HDF file'
#                               ' and alignment file.')
#                         return
#                 i += 1
#             aligned_dataframes[name][chain] = out_df
#     return aligned_dataframes


def align_residues(
    residue_dataframes: AlignedResidueFrames,
    alignment: Mapping[str, Any],
    mapping: Mapping[str, str] | None = None,
) -> AlignedResidueFrames | None:
    """Align residue dataframes based on sequence alignment"""
    aligned_dataframes = residue_dataframes
    for chain, alignment_file in alignment["alignment"].items():
        aln_dict = parse_alignment(alignment_file)
        if aln_dict is None:
            return
        # Length of sequences in alignment file
        aln_seq_length = len(next(iter(aln_dict.values())))

        for name, seq_label in alignment["names"].items():
            df = residue_dataframes[name][chain]
            nan_array = numpy.empty((aln_seq_length, df.shape[1]))
            nan_array[:] = numpy.nan
            out_df = pandas.DataFrame(data=nan_array, columns=df.columns)
            out_df["sequence", "resn"] = "-"
            i, j = 0, 0
            for residue in aln_dict[seq_label]:
                if residue != "-":
                    if ONE_TO_THREE[residue] == df.iloc[j].sequence.resn:
                        out_df.iloc[i] = df.iloc[j]
                        j += 1
                    else:
                        print("Mismatch between sequence in HDF file and alignment file.")
                        return
                i += 1
            aligned_dataframes[name][chain] = out_df.reset_index()
    return aligned_dataframes


def main(
    input_files: Sequence[str], output: str, alignment: str | None = None, annotations: str | None = None
) -> None:
    res_dfs = collections.defaultdict(dict)
    for input_file in input_files:
        with h5py.File(input_file, "r") as hdf_file:
            res_df = residue_dataframe(hdf_file=hdf_file)
            res_dfs["data"][hdf_file.attrs["name"]] = res_df
    if annotations:
        with open(annotations) as annotations_file:
            res_dfs["annotations"] = load_toml(annotations_file)
    if alignment:
        alignment = load_toml(alignment)
        res_dfs["data"] = align_residues(res_dfs["data"], alignment=alignment)
    pickle.dump(res_dfs, open(output, "wb"))


if __name__ == "__main__":
    import sys

    from md_davis.cli import main as cli_main

    raise SystemExit(cli_main(["residue", *sys.argv[1:]]))
