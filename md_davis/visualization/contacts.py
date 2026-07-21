"""Formatting helpers for contact and hydrogen-bond plots."""

from __future__ import annotations

from collections.abc import Iterable


def contact_labels(atoms: Iterable[object]) -> list[str]:
    """Format atom tuples as labels suitable for plot axes."""
    return ["-".join(str(part) for part in atom) if atom else "" for atom in atoms]
