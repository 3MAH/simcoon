"""
Constant class to manage list of solids belonging to the same phase
"""

import os
import shutil
from typing import List, NamedTuple
import numpy as np


class Constant:
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
    input_files: List[str]

    def __init__(
        self,
        number: Optional[int] = None,
        ninput_values: int = 1,
        key: Optional[str] = None,
    ) -> None:

        self.number = number
        self.input_values = np.zeros(ninput_values)
        self.key = key
        self.input_files: List[str] = []


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
            const = Constant(ninput_values=n_consts)
            const.number = int(values[0])
            const.key = values[1]
            for j in range(n_consts):
                const.input_values[j] = values[2 + j]
            nfiles = int(values[2 + n_consts])
            for j in range(nfiles):
                const.input_files.append(values[3 + n_consts + j])
            consts.append(const)
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
