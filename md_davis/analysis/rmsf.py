from __future__ import annotations

from os import PathLike

import numpy
import pandas


def rmsf_to_df(xvg_file: str | PathLike[str]) -> pandas.DataFrame:
    """Convert . xvg file containing RMSF to Pandas Data Frame"""
    data = numpy.loadtxt(xvg_file, comments=("#", "@"))
    df = pandas.DataFrame(data, columns=("Residue", "RMSF"))
    df["Residue"] = df["Residue"].astype(int)
    return df
