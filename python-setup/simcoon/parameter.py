"""
Phase class to manage list of solids belonging to the same phase
"""

import os
import shutil
from typing import List, NamedTuple


class Parameter(NamedTuple):
    """
    Parameter class to manage of set of parameters values to be applied during an identification
    and/or DOE test matrix
    :param number: int, number of the constant
    :param ninput_values: number of inputs for the considered constant
    :param key: alphanumeric key to identify the constant in a file
    """

    number: int
    key: str
    value: float
    min_value: float
    max_value: float
    input_files: List[str]


def read_parameters() -> List[Parameter]:
    """
    read_parameters from a simcoon input file
    @return : List of Parameter
    """
    params = []
    with open("data/parameters.inp", "r", encoding="utf-8") as paraminit:
        lines = paraminit.readlines()

        for line in lines[1:]:
            values = line.split()
            nfiles = int(values[4])
            pa = Parameter(
                number=int(values[0]),
                value=0.5 * float(values[1]) + 0.5 * float(values[2]),
                min_value=float(values[1]),
                max_value=float(values[2]),
                key=values[3],
                input_files=[values[4 + j] for j in range(nfiles)],
            )
            params.append(pa)
    return params


def copy_parameters(
    params: List[Parameter],
    src_path: str,
    dst_path: str,
) -> None:
    """
    Copy parameters file from a source path to a destination path, so that key values can be applied
    :param params: List of Parameter objects
    :param src_path: Source path
    :param dst_path: Destination path
    :return: None
    """
    for pa in params:
        for ifiles in pa.input_files:
            src_files = os.path.join(src_path, ifiles)
            dst_files = os.path.join(dst_path, ifiles)
            shutil.copy(src_files, dst_files)


def apply_parameter(
    params: List[Parameter],
    dst_path: str,
) -> None:
    """
    Apply parameters, i.e. replace in file the value of the key with the value of the parameter
    :param params: List of Parameter objects
    :param src_path: Source path
    :param dst_path: Destination path
    :return: None
    """
    for pa in params:
        for ifiles in pa.input_files:
            mod_files = os.path.join(dst_path, ifiles)

            with open(mod_files, "r", encoding="utf-8") as in_files:
                content = in_files.read()
            modified_content = content.replace(pa.key, str(pa.value))
            with open(mod_files, "w", encoding="utf-8") as ou_files:
                ou_files.write(modified_content)