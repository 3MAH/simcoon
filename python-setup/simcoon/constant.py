"""
Constant class to manage simcoon computations constants
"""

import os
from typing import List, NamedTuple, Union
import numpy as np


class Constant(NamedTuple):
    """
    Constant class to manage of set of constant values to be applied during an
    identification and/or DOE test matrix
    :param number: int, number of the constant
    :param ninput_values: number of inputs for the considered constant
    :param key: alphanumeric key to identify the constant in a file
    """

    number: int
    key: str
    input_values: np.ndarray
    value: float
    sim_input_files: List[str]


def read_constants(
    n_consts: int,
    fname: Union[str, os.PathLike] = "data/constants.inp",
    skiprows: int = 1,
) -> List[Constant]:
    """
    Read constants from a simcoon input file.

    :param n_consts: number of input values per constant
    :param fname: path to the ``constants.inp`` file
    :param skiprows: number of leading lines to skip (pandas-style).
                     Defaults to 1 (skip the header). Use ``0`` for files
                     without a header.
    :return: List of Constant
    """

    if not isinstance(fname, (str, os.PathLike)):
        raise TypeError(
            f"Invalid type: {type(fname).__name__}. Expected str or os.PathLike."
        )
    if skiprows < 0:
        raise ValueError(f"skiprows must be >= 0, got {skiprows}")

    if isinstance(fname, os.PathLike):
        fname = os.fspath(fname)

    consts = []
    with open(fname, "r", encoding="utf-8") as constinit:
        lines = constinit.readlines()

        for line in lines[skiprows:]:
            values = line.split()
            array_input_values = np.zeros(n_consts)
            nfiles = int(values[2 + n_consts])
            for j in range(n_consts):
                array_input_values[j] = values[2 + j]

            const = Constant(
                number=int(values[0]),
                key=values[1],
                input_values=array_input_values,
                value=array_input_values[0],
                sim_input_files=[values[3 + n_consts + j] for j in range(nfiles)],
            )
            consts.append(const)
    return consts

