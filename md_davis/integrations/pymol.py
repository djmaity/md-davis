"""Lazy access to the conda-forge ``pymol-open-source`` package."""

from __future__ import annotations

from collections.abc import Sequence
from types import ModuleType
from typing import Any


class PyMOLUnavailableError(ImportError):
    """Raised when a PyMOL feature is used without open-source PyMOL."""


def require_pymol() -> ModuleType:
    """Import PyMOL only when a PyMOL-backed feature is requested."""
    try:
        import pymol
    except ImportError as exc:  # pragma: no cover - depends on optional package
        raise PyMOLUnavailableError(
            "This feature requires pymol-open-source. Install the project "
            "Conda environment or run `conda install -c conda-forge "
            "pymol-open-source`."
        ) from exc
    return pymol


def initialize_pymol(
    window: bool = False,
    extra_args: Sequence[str] = (),
) -> Any:
    """Start open-source PyMOL and return its command API."""
    pymol = require_pymol()
    pymol.pymol_argv = ["pymol"]
    if not window:
        pymol.pymol_argv.append("-cQ")
    pymol.pymol_argv.extend(extra_args)
    pymol.finish_launching()
    return pymol.cmd
