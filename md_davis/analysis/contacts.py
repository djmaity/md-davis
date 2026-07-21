"""
Parse contacts evaluated by gmx hbonds

Usage:
  md_davis contacts [options] (--file <.xpm>)
                              (--index <.ndx>)
                              (--structure <.pdb/.gro>)
                              (--group <string>)
  md_davis contacts -h | --help

Options:

  -f, --file <.xpm>             Contact file obtained from GROMACS -i,
  --index <.ndx>                Index file
  -s, --structure <.pdb/.gro>   Structure file
  -g, --group <string>          Group to match from index file to get the list of contacts

  -b, --begin <int>             Frame to start calculation from

  --pickle FILENAME             Save the output to a pickle file
  --hdf FILENAME                Save the output to a HDF file
  --csv FILENAME                Save the output to a CSV file

  --matrix FILENAME             Save a file containing the contact matrix

  --pdb FILENAME                Save a pdb file with percentage of frames
                                in the B-factor column

  -h, --help                    Show this screen.
"""

from __future__ import annotations

import collections
import operator
import pickle
import re
from collections.abc import Iterator, Sequence
from os import PathLike
from typing import Any

import numpy
import pandas
from biopandas.pdb import PandasPdb
from plotly import graph_objects as go
from plotly.offline import plot as plot_html

from md_davis.visualization.contacts import contact_labels

# TODO: Directly use biopandas data frame and remove the atom class


class Atom:
    """Wrapper for MDTraj structure.atom object"""

    def __init__(self, index: int, structure: PandasPdb) -> None:
        atom = structure.df["ATOM"]
        self.index = index
        self.segment = atom.iloc[index]["segment_id"]
        self.chain = atom.iloc[index]["chain_id"]
        self.residue = atom.iloc[index]["residue_name"]
        self.res_id = atom.iloc[index]["residue_number"]
        self.name = atom.iloc[index]["atom_name"]

    def __repr__(self) -> str:
        return f"{self.segment:2} {self.chain} {self.residue} {self.res_id:3} {self.name:4}"

    def __str__(self) -> str:
        return f"{self.segment:2} {self.chain} {self.residue} {self.res_id:3} {self.name:4}"

    def __gt__(self, other: Atom) -> bool:
        for prop1, prop2 in zip(self, other, strict=True):
            if prop1 == prop2:
                continue
            elif prop1 > prop2:
                return True
            else:
                return False
        return False

    def __ge__(self, other: Atom) -> bool:
        if self == other or self > other:
            return True
        return False

    def __iter__(self) -> Iterator[Any]:
        """Segment Chain Residue ResSeq Atom"""
        return iter([self.segment, self.chain, self.residue, self.res_id, self.name])


class Contact:
    """Defines a contact between two atoms"""

    def __init__(self, atom_indices: Sequence[int], structure: PandasPdb) -> None:
        self.atom1 = Atom(atom_indices[0], structure)
        self.atom2 = Atom(atom_indices[-1], structure)
        self.count = 0
        self.time_series = []

    def __repr__(self) -> str:
        return f"\n{repr(self.atom1)} -- {repr(self.atom2)}: {self.count}"

    def __str__(self) -> str:
        return f"\n{str(self.atom1)} -- {str(self.atom2)}: {self.count}"

    def __iter__(self) -> Iterator[Any]:
        return iter(list(self.atom1) + list(self.atom2) + [self.count])


class Contacts:
    """All hydrogen bonds"""

    def __init__(self, structure: str | PathLike[str]) -> None:
        self.structure = PandasPdb().read_pdb(structure)
        self.bonds = []
        self.nframes = 1

    def __repr__(self) -> Any:
        return self.bonds

    def __str__(self) -> str:
        out_str = "Segment Chain Residue ResSeq Atom -- Segment Chain Residue ResSeq Atom : Count\n"
        for bond in self.bonds:
            out_str += bond.__repr__()  # + '\n'
        return out_str

    def __iter__(self) -> Iterator[list[Any]]:
        return iter([list(_) for _ in self.bonds])

    def __len__(self) -> int:
        return len(self.bonds)

    def parse_indices(self, index_file: str | PathLike[str], group: str) -> None:
        save = False
        with open(index_file) as ndx_file:
            for line in ndx_file:
                if save:
                    indices = [int(_) - 1 for _ in line.strip().split()]
                    self.bonds.append(Contact(indices, self.structure))
                if line == f"[ {group} ]\n":
                    save = True
        if not save:
            print("Specified group not found in index file")

    def add_counts(self, xpm_file: str | PathLike[str]) -> None:
        bond_idx = 0
        toggle = False
        with open(xpm_file) as xpm_data:
            for line in xpm_data:
                if line.startswith("/* x-axis:"):
                    toggle = True
                if toggle:
                    if line.startswith('"'):  # Read data
                        contacts = re.search('"(.*)"', line).group(1)
                        bond_seires_array = [1 if char == "o" else 0 for char in contacts]
                        self.bonds[bond_idx].time_series = bond_seires_array
                        self.bonds[bond_idx].count = sum(bond_seires_array)
                        bond_idx += 1

    @property
    def to_df(self) -> pandas.DataFrame:
        columns = [
            "Segment1",
            "Chain1",
            "Residue1",
            "ResSeq1",
            "Atom1",
            "Segment2",
            "Chain2",
            "Residue2",
            "ResSeq2",
            "Atom2",
            "Count",
        ]
        df = pandas.DataFrame(data=list(self), columns=columns)
        return df

    def to_pdb(self, filename: str | PathLike[str]) -> None:
        contact_dict = collections.defaultdict(list)
        for bond in self.bonds:
            contact_dict[bond.atom1.index].append(bond.time_series)
            contact_dict[bond.atom2.index].append(bond.time_series)
        for index, contacts in contact_dict.items():
            contacts = numpy.array(contacts).any(axis=0)
            b_factor = round(numpy.sum(contacts) / len(contacts) * 99.99, 2)
            self.structure.df["ATOM"].loc[index, "b_factor"] = b_factor
        self.structure.to_pdb(path=filename, records=None, gz=False, append_newline=True)

    def plot_matrix(
        self,
        df: pandas.DataFrame | str | bytes | None = None,
        n_frames: int = 100,
        cutoff: float = 0,
        title: str = "Contact matrix",
        output: str | PathLike[str] | None = None,
    ) -> go.Figure:
        """Create a contact heatmap and optionally save it as standalone HTML."""
        if isinstance(df, (str, bytes)):
            df = pandas.read_pickle(df)
        elif df is None:
            df = self.to_df()

        # Convert to percentage of frames; default is to return counts
        df["Count"] = df["Count"] * 100 / n_frames

        # counts or percentage above which to show contacts
        df = df[df["Count"] > cutoff]

        atoms1 = list(df.groupby(["Chain1", "Residue1", "ResSeq1", "Atom1"]).groups.keys())
        atoms1 = sorted(atoms1, key=operator.itemgetter(0, 2, 1))
        atom_group2 = df.groupby(["Chain2", "Residue2", "ResSeq2", "Atom2"])
        atoms2 = list(atom_group2.groups.keys())
        atoms2 = sorted(atoms2, key=operator.itemgetter(0, 2, 1))

        matrix = numpy.empty(shape=(len(atoms1), len(atoms2)))
        matrix[:] = numpy.nan
        for group, atoms in atom_group2:
            j = atoms2.index(group)
            for _, row in atoms.iterrows():
                i = atoms1.index((row["Chain1"], row["Residue1"], row["ResSeq1"], row["Atom1"]))
                matrix[i, j] = row["Count"]

        fig = go.Figure(
            data=go.Heatmap(
                z=matrix,
                x=contact_labels(atoms1),
                y=contact_labels(atoms2),
                colorscale="Greys",
                zmin=0,
                zmax=100,
                transpose=True,
            )
        )
        fig.update_layout(title=title)
        if output:
            plot_html(fig, filename=output, auto_open=False)
        return fig


class Hbond:
    """Defines a hydrogen bond"""

    def __init__(self, atom_indices: Sequence[int], structure: PandasPdb) -> None:
        self.donor = Atom(atom_indices[0], structure)
        self.hydrogen = Atom(atom_indices[1], structure)
        self.recipient = Atom(atom_indices[2], structure)
        self.count = 0
        self.time_series = []

    def __repr__(self) -> str:
        return f"\n{self.donor} -- {self.hydrogen} -> {self.recipient}: {self.count}"

    def __iter__(self) -> Iterator[Any]:
        return iter(list(self.donor) + list(self.hydrogen) + list(self.recipient) + [self.count])

    def __eq__(self, other: Hbond) -> bool:
        if self.donor == other.donor and self.hydrogen == other.hydrogen and self.recipient == other.recipient:
            return True
        else:
            return False

    def __sub__(self, other: Hbond) -> Hbond:
        """Calculates the frames in which the first bond is present and
        not the second bond"""
        output = Hbond()
        output.donor = self.donor
        output.hydrogen = self.hydrogen
        output.recipient = self.recipient
        tmp = numpy.array(self.time_series) - numpy.array(other.time_series)
        output.time_series = tmp
        output.count = sum(tmp)
        return output

    def __neg__(self, other: Hbond) -> Hbond:
        """Calculates the frames in which the first bond is present and
        not the second bond"""
        output = Hbond()
        output.donor = self.donor
        output.hydrogen = self.hydrogen
        output.recipient = self.recipient
        output.time_series = -numpy.array(self.time_series)
        output.count = -self.count
        return output


# Subclass the Contacts to get Hbonds Class
class Hbonds(Contacts):
    """All hydrogen bonds"""

    def __init__(self, structure: str | PathLike[str]) -> None:
        self.structure = PandasPdb().read_pdb(structure)
        self.bonds = []

    def __repr__(self) -> Any:
        return self.bonds

    def __str__(self) -> str:
        out_str = "Seg Ch Res ResSeq Atom -- Seg Ch Res ResSeq Atom --> Seg Ch Res ResSeq Atom : Count\n"
        for bond in self.bonds:
            out_str += bond.__repr__()  # + '\n'
        return out_str

    def __iter__(self) -> Iterator[list[Any]]:
        return iter([list(_) for _ in self.bonds])

    def __sub__(self, other: Hbonds) -> None:
        output = Hbonds()
        output.structure = self.structure
        for bond1 in self.bonds:
            for bond2 in other.bonds:
                if bond1 == bond2:
                    output.bonds.append(bond1 - bond2)
                else:
                    output.bonds.append(bond1)
                    output.bonds.append(-bond2)

    def parse_indices(self, index_file: str | PathLike[str], group: str) -> None:
        save = False
        with open(index_file) as ndx_file:
            for line in ndx_file:
                if save:
                    indices = [int(_) - 1 for _ in line.strip().split()]
                    self.bonds.append(Hbond(indices, self.structure))
                if line == f"[ {group} ]\n":
                    save = True

    @property
    def to_df(self) -> pandas.DataFrame:
        columns = [
            "Segment1",
            "Chain1",
            "Residue1",
            "ResSeq1",
            "Atom1",
            "Segment1-H",
            "Chain1-H",
            "Residue1-H",
            "ResSeq1-H",
            "Atom1-H",
            "Segment2",
            "Chain2",
            "Residue2",
            "ResSeq2",
            "Atom2",
            "Count",
        ]
        df = pandas.DataFrame(data=list(self), columns=columns)
        return df


class Contact_Matrix:
    """All hydrogen bonds"""

    def __init__(self, structure: str | PathLike[str]) -> None:
        self.atoms1 = PandasPdb().read_pdb(structure)
        self.atoms2 = PandasPdb().read_pdb(structure)
        self.matrix = []


def main(
    xpm_file: str,
    index_file: str,
    structure: str,
    group: str,
    contacts: bool,
    silent: bool = False,
    output: str | None = None,
    save_pickle: str | None = None,
    save_hdf: str | None = None,
    save_csv: str | None = None,
    save_pdb: str | None = None,
    save_matrix: str | None = None,
) -> None:
    """Parse contacts evaluated by gmx hbonds"""
    if contacts:
        molecular_contacts = Contacts(structure)
    else:
        molecular_contacts = Hbonds(structure)

    molecular_contacts.parse_indices(index_file=index_file, group=group)
    molecular_contacts.add_counts(xpm_file=xpm_file)

    if not silent:
        print(molecular_contacts)

    if output:
        pickle.dump(molecular_contacts, open(output, "wb"))
    df = molecular_contacts.to_df
    if save_hdf:
        df.to_hdf(save_hdf, key="contacts")
    if save_csv:
        df.to_csv(save_csv)
    if save_pickle:
        df.to_pickle(save_pickle)
    if save_pdb:
        molecular_contacts.to_pdb(save_pdb)
    if save_matrix:
        molecular_contacts.plot_matrix(output=save_matrix)


if __name__ == "__main__":
    import sys

    from md_davis.cli import main as cli_main

    raise SystemExit(cli_main(["hbond", *sys.argv[1:]]))
