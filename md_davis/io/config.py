"""Configuration readers."""

from __future__ import annotations

import ast
import tomllib
from os import PathLike
from typing import Any, BinaryIO, TextIO, overload


def load_toml(source: str | PathLike[str] | BinaryIO | TextIO) -> dict[str, Any]:
    """Load TOML from a path or an already-open file object."""
    if hasattr(source, "read"):
        content = source.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8")
        return tomllib.loads(content)

    with open(source, "rb") as stream:
        return tomllib.load(stream)


@overload
def parse_literal(value: Any, expected_type: None = None) -> Any: ...


@overload
def parse_literal[T](value: Any, expected_type: type[T] | tuple[type[T], ...]) -> T | None: ...


def parse_literal(value: Any, expected_type: type[Any] | tuple[type[Any], ...] | None = None) -> Any:
    """Parse a command-line Python literal without executing arbitrary code."""
    if value is None or not isinstance(value, str):
        parsed = value
    else:
        try:
            parsed = ast.literal_eval(value)
        except (SyntaxError, ValueError) as error:
            raise ValueError(f"Invalid literal: {value!r}") from error

    if expected_type is not None and parsed is not None and not isinstance(parsed, expected_type):
        expected_types = expected_type if isinstance(expected_type, tuple) else (expected_type,)
        expected_name = " or ".join(candidate.__name__ for candidate in expected_types)
        raise ValueError(f"Expected {expected_name}, received {type(parsed).__name__}")
    return parsed
