"""Sequence-alignment parsers."""

from __future__ import annotations

import re
from os import PathLike

_ALIGNMENT_LINE = re.compile(r"^(?P<label>\S+)\s+(?P<sequence>[A-Za-z-]+)(?:\s+\d+)?\s*$")


def parse_clustal(filename: str | PathLike[str]) -> dict[str, str]:
    """Return complete sequences from a CLUSTAL-style alignment file."""
    sequences: dict[str, str] = {}
    with open(filename, encoding="utf-8") as alignment_file:
        for line in alignment_file:
            match = _ALIGNMENT_LINE.match(line)
            if match is None or match.group("label").upper() == "CLUSTAL":
                continue
            label = match.group("label")
            sequences[label] = sequences.get(label, "") + match.group("sequence")

    if not sequences:
        raise ValueError(f"No aligned sequences found in {filename}")

    lengths = {len(sequence) for sequence in sequences.values()}
    if len(lengths) != 1:
        raise ValueError(f"Alignment sequences have inconsistent lengths in {filename}")
    return sequences
