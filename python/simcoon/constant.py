"""
Constant class to manage simcoon computations constants
"""

import os
import shutil
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
    n_consts: int, fname: Union[str, os.PathLike] = "data/constants.inp"
) -> List[Constant]:
    """
    read constants from a simcoon input file
    :param n_const: number of constants
    :param path: path where constants.inp simcoon input file is located
    :return: List of Constant
    """

    if not isinstance(fname, (str, os.PathLike)):
        raise TypeError(
            f"Invalid type: {type(fname).__name__}. Expected str or os.PathLike."
        )

    if isinstance(fname, os.PathLike):
        fname = os.fspath(fname)

    consts = []
    with open(fname, "r", encoding="utf-8") as constinit:
        lines = constinit.readlines()

        for line in lines[1:]:
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


def copy_constants(
    consts: List[Constant],
    src_path: Union[str, os.PathLike],
    dst_path: Union[str, os.PathLike],
) -> None:
    """
    Copy constants file from a source path to a destination path, so that key values can be applied
    :param consts: List of Constant objects
    :param src_path: Source path
    :param dst_path: Destination path
    :return: None
    """
    if not isinstance(src_path, (str, os.PathLike)):
        raise TypeError(
            f"Invalid type: {type(src_path).__name__}. Expected str or os.PathLike."
        )

    if not isinstance(dst_path, (str, os.PathLike)):
        raise TypeError(
            f"Invalid type: {type(dst_path).__name__}. Expected str or os.PathLike."
        )

    for co in consts:
        if not all(isinstance(item, str) for item in co.sim_input_files):
            raise TypeError("All elements in sim_input_files must be strings.")

        for ifiles in co.sim_input_files:
            src_files = os.path.join(src_path, ifiles)
            dst_files = os.path.join(dst_path, ifiles)
            shutil.copy(src_files, dst_files)


def apply_constants(
    consts: List[Constant],
    dst_path: Union[str, os.PathLike],
) -> None:
    """
    Apply constants, i.e. replace in file the value of the key with the value of the Constant
    :param params: List of Constant objects
    :param src_path: Source path
    :param dst_path: Destination path
    :return: None
    """
    if not isinstance(dst_path, (str, os.PathLike)):
        raise TypeError(
            f"Invalid type: {type(dst_path).__name__}. Expected str or os.PathLike."
        )

    for co in consts:
        for ifiles in co.sim_input_files:
            mod_files = os.path.join(dst_path, ifiles)

            with open(mod_files, "r", encoding="utf-8") as in_files:
                content = in_files.read()
            modified_content = content.replace(const.key, str(const.value))
            with open(mod_files, "w", encoding="utf-8") as ou_files:
                ou_files.write(modified_content)
