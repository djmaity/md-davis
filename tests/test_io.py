"""Tests for shared file-format parsers."""

from io import StringIO

import pytest

from md_davis.io.alignment import parse_clustal
from md_davis.io.config import load_toml, parse_literal
from md_davis.io.delphi import (
    read_delphi_potential,
    residue_potential_plot_data,
    residue_potential_series,
)


def test_load_toml_from_path_and_stream(tmp_path):
    filename = tmp_path / "input.toml"
    filename.write_text('name = "simulation"\n[options]\nframes = 10\n', encoding="utf-8")

    expected = {"name": "simulation", "options": {"frames": 10}}
    assert load_toml(filename) == expected
    assert load_toml(StringIO(filename.read_text(encoding="utf-8"))) == expected


def test_parse_literal_validates_type_without_executing_code():
    assert parse_literal("{'x': [0, 1]}", dict) == {"x": [0, 1]}
    with pytest.raises(ValueError, match="Expected dict"):
        parse_literal("[1, 2]", dict)
    with pytest.raises(ValueError, match="Invalid literal"):
        parse_literal("__import__('os').getcwd()")


def test_parse_multiblock_clustal_alignment(tmp_path):
    filename = tmp_path / "alignment.aln"
    filename.write_text(
        "CLUSTAL W multiple sequence alignment\n\n"
        "protein_a  AC-D  4\n"
        "protein_b  ACED  4\n\n"
        "protein_a  EF\n"
        "protein_b  EF\n",
        encoding="utf-8",
    )

    assert parse_clustal(filename) == {
        "protein_a": "AC-DEF",
        "protein_b": "ACEDEF",
    }


def test_parse_clustal_rejects_inconsistent_lengths(tmp_path):
    filename = tmp_path / "bad.aln"
    filename.write_text("a  AAA\nb  AA\n", encoding="utf-8")

    with pytest.raises(ValueError, match="inconsistent lengths"):
        parse_clustal(filename)


def _potential_row(name, residue_name, chain, residue, potential):
    numeric = [potential, 0.0, 0.0, 0.0, 0.0, 0.0]
    return f"{name:<5}{residue_name:<3}{chain:<3}{residue:>9}" + "".join(f"{value:>10.3f}" for value in numeric) + "\n"


def test_delphi_potential_aggregation(tmp_path):
    filename = tmp_path / "potential.pot"
    filename.write_text(
        "".join(f"header {index}\n" for index in range(12))
        + _potential_row("CA", "ALA", "A", 1, 1.25)
        + _potential_row("CB", "ALA", "A", 1, 0.75)
        + _potential_row("CA", "GLY", "A", 2, -0.5)
        + "footer\nfooter\n",
        encoding="utf-8",
    )

    data = read_delphi_potential(filename)
    assert data["resSeq"].tolist() == [1, 1, 2]
    assert residue_potential_series(filename).to_dict() == {1: 2.0, 2: -0.5}
    assert residue_potential_plot_data(filename) == (
        [1, 2],
        ["(1, 'ALA')", "(2, 'GLY')"],
        [2.0, -0.5],
    )
