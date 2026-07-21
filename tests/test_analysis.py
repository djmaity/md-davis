"""Tests for pure numerical and parsing behavior."""

import numpy as np
import pytest

from md_davis.analysis.landscape.landscape import Landscape
from md_davis.analysis.statistics import moving_avg, moving_avg_std, rolling_std
from md_davis.core.collation import split_chains, split_increasing
from md_davis.io.xvg import Xvg


def test_split_increasing_and_split_chains():
    values = np.array([1, 2, 3, 1, 2])
    segments = split_increasing(values)
    assert [segment.tolist() for segment in segments] == [[1, 2, 3], [1, 2]]

    chains = split_chains(np.arange(6), [2, 4])
    assert [chain.tolist() for chain in chains] == [[0, 1], [2, 3, 4, 5]]


def test_rolling_statistics_share_window_semantics():
    values = np.array([1.0, 2.0, 3.0, 4.0])
    means, deviations = moving_avg_std(values, window=2, shift=1)

    np.testing.assert_allclose(means, moving_avg(values, window=2, shift=1))
    np.testing.assert_allclose(deviations, rolling_std(values, window=2, shift=1))
    np.testing.assert_allclose(means, [1.5, 2.5, 3.5, 4.0])


def test_landscape_round_trip(tmp_path):
    landscape = Landscape.landscape(
        "sample",
        time=np.arange(4),
        x_data=np.array([0.0, 1.0, 2.0, 3.0]),
        y_data=np.array([3.0, 2.0, 1.0, 0.0]),
        shape=(4, 4),
    )
    assert landscape.zValues.sum() == 4
    assert Landscape.minmax([np.nan, 3.0, -1.0]) == (-1.0, 3.0)

    filename = tmp_path / "landscape.h5"
    landscape.save(filename, "sample")
    loaded = Landscape.open(filename)["sample"]
    np.testing.assert_array_equal(loaded.zValues, landscape.zValues)
    assert sum(len(times) for _, _, times in loaded) == 4


def test_landscape_rejects_unequal_inputs():
    with pytest.raises(ValueError, match="Unequal length"):
        Landscape.landscape("bad", [0, 1], [0, 1], [0])


def test_xvg_parser_handles_headers_legends_and_blank_lines(tmp_path):
    filename = tmp_path / "rmsd.xvg"
    filename.write_text(
        '@ title "RMSD"\n'
        '@ subtitle "simulation"\n'
        '@ xaxis  label "Time"\n'
        '@ yaxis  label "RMSD"\n'
        '@ s0 legend "Backbone"\n\n'
        "0 0.1\n1 0.2\n",
        encoding="utf-8",
    )

    data = Xvg(filename)
    assert data.title == "RMSD"
    assert data.legend == {0: "Time", 1: "Backbone"}
    np.testing.assert_allclose(data.data, [[0.0, 0.1], [1.0, 0.2]])
    assert "title = RMSD" in str(data)
