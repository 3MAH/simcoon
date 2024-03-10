"""
Constant class to manage list of solids belonging to the same phase
"""

import os
import shutil
from typing import List, NamedTuple
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
    input_files: List[str]


def read_constants(n_consts: int) -> List[Constant]:
    """
    read_constants from a simcoon input file
    @param: n_const
    @return : List of Constant
    """
    consts = []
    with open("data/constants.inp", "r", encoding="utf-8") as paraminit:
        lines = paraminit.readlines()

        for line in lines[1:]:
            values = line.split()
            array_input_values = np.zeros(n_consts)
            nfiles = int(values[2 + n_consts])
            for j in range(n_consts):
                array_input_values[j] = values[2 + j]

            co = Constant(
                number=int(values[0]),
                key=values[1],
                input_values=array_input_values,
                value=array_input_values[0],
                input_files=[values[3 + n_consts + j] for j in range(nfiles)],
            )
            consts.append(co)
    return consts


def copy_constants(
    consts: List[Constant],
    src_path: str,
    dst_path: str,
) -> None:
    """
    Copy constants file from a source path to a destination path, so that key values can be applied
    :param consts: List of Constant objects
    :param src_path: Source path
    :param dst_path: Destination path
    :return: None
    """
    for co in consts:
        for ifiles in co.input_files:
            src_files = os.path.join(src_path, ifiles)
            dst_files = os.path.join(dst_path, ifiles)
            shutil.copy(src_files, dst_files)


def apply_constants(
    consts: List[Constant],
    dst_path: str,
) -> None:
    """
    Apply constants, i.e. replace in file the value of the key with the value of the Constant
    :param params: List of Constant objects
    :param src_path: Source path
    :param dst_path: Destination path
    :return: None
    """
    for co in consts:
        for ifiles in co.input_files:
            mod_files = os.path.join(dst_path, ifiles)

            with open(mod_files, "r", encoding="utf-8") as in_files:
                content = in_files.read()
            modified_content = content.replace(co.key, str(co.value))
            with open(mod_files, "w", encoding="utf-8") as ou_files:
                ou_files.write(modified_content)
