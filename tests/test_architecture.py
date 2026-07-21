"""Tests that enforce the canonical package hierarchy."""

from importlib import import_module
from importlib.util import find_spec

import pytest

CANONICAL_PACKAGES = (
    "md_davis.analysis",
    "md_davis.analysis.electrostatics",
    "md_davis.analysis.landscape",
    "md_davis.analysis.structure",
    "md_davis.cli",
    "md_davis.core",
    "md_davis.integrations",
    "md_davis.io",
    "md_davis.visualization",
    "md_davis.visualization.plotting",
)

LEGACY_MODULES = (
    "md_davis.amino_acids",
    "md_davis.collate",
    "md_davis.common",
    "md_davis.electrostatics",
    "md_davis.hbond",
    "md_davis.landscape",
    "md_davis.plotting",
    "md_davis.process_gromacs_trajectory",
    "md_davis.residue",
    "md_davis.sequence",
    "md_davis.structure",
    "md_davis.xvg",
)


@pytest.mark.parametrize("module_name", CANONICAL_PACKAGES)
def test_canonical_package_imports(module_name):
    assert import_module(module_name).__name__ == module_name


@pytest.mark.parametrize("module_name", LEGACY_MODULES)
def test_legacy_modules_are_removed(module_name):
    assert find_spec(module_name) is None
