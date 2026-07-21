"""Small rolling-statistics helpers."""

from __future__ import annotations

import numpy
from numpy.typing import ArrayLike, NDArray


def _window_values(data: ArrayLike, window: int = 1, shift: int = 1) -> list[NDArray]:
    if window < 1 or shift < 1:
        raise ValueError("window and shift must be positive integers")
    values = numpy.asarray(data)
    if values.size == 0:
        return []
    return [values[start : start + window] for start in range(0, len(values), shift)]


def moving_avg(data: ArrayLike, window: int = 1, shift: int = 1) -> NDArray:
    """Return the mean of each rolling window, including the final partial window."""
    return numpy.asarray([numpy.mean(values) for values in _window_values(data, window, shift)])


def rolling_std(data: ArrayLike, window: int = 1, shift: int = 1) -> NDArray:
    """Return the standard deviation of each rolling window."""
    return numpy.asarray([numpy.std(values) for values in _window_values(data, window, shift)])


def moving_avg_std(data: ArrayLike, window: int = 1, shift: int = 1) -> tuple[NDArray, NDArray]:
    """Return rolling means and standard deviations calculated from shared windows."""
    windows = _window_values(data, window, shift)
    means = numpy.asarray([numpy.mean(values) for values in windows])
    deviations = numpy.asarray([numpy.std(values) for values in windows])
    return means, deviations
