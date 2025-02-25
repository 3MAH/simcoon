import os
import shutil
from typing import List, Optional, Tuple, Union


class Parameter:
    """
    Parameter class to manage of set of parameters values to be applied during an identification
    and/or DOE test matrix
    :param number: int, number of the parameter
    :param bounds: min and max boundaries in which parameter may vary
    :param key: alphanumeric key to identify the parameter in a file
    :param sim_input_files: simulation input files in which the parameter is used
    """

    def __init__(
        self,
        number: int = 0,
        bounds: Tuple[float, float] = (0.0, 1.0),
        key: str = "",
        sim_input_files: Optional[List[str]] = None,
    ) -> None:
        self.number = number
        self.key = key
        self.bounds = bounds
        self._min_value = bounds[0]
        self._max_value = bounds[1]
        self.sim_input_files = sim_input_files

        self._value = None

    @property
    def value(self) -> float:
        if self._value is not None:
            return self._value
        else:
            return 0.5 * (self._min_value + self._max_value)

    @value.setter
    def value(self, value):
        self._value = value


def read_parameters(
    fname: Union[str, os.PathLike] = "data/parameters.inp",
) -> List[Parameter]:
    """
    read_parameters from a simcoon input file located in path
    :param path: path where parameters.inp simcoon input file is located
    :return: List of Parameter
    """
    if not isinstance(fname, (str, os.PathLike)):
        raise TypeError(
            f"Invalid type: {type(fname).__name__}. Expected str or os.PathLike."
        )

    if isinstance(fname, os.PathLike):
        fname = os.fspath(fname)

    params = []
    with open(fname, "r", encoding="utf-8") as paraminit:
        lines = paraminit.readlines()

        for line in lines[1:]:
            values = line.split()
            nfiles = int(values[4])
            param = Parameter(
                number=int(values[0]),
                bounds=(float(values[1]), float(values[2])),
                key=values[3],
                sim_input_files=[values[5 + j] for j in range(nfiles)],
            )
            params.append(param)
    return params


def copy_parameters(
    params: List[Parameter],
    src_path: Union[str, os.PathLike],
    dst_path: Union[str, os.PathLike],
) -> None:
    """
    Copy parameters file from a source path to a destination path, so that key values can be applied
    :param params: List of Parameter objects
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

    for param in params:
        if not all(isinstance(item, str) for item in param.sim_input_files):
            raise TypeError("All elements in sim_input_files must be strings.")

        for ifiles in param.sim_input_files:
            src_files = os.path.join(src_path, ifiles)
            dst_files = os.path.join(dst_path, ifiles)
            shutil.copy(src_files, dst_files)


def apply_parameters(
    params: List[Parameter],
    dst_path: Union[str, os.PathLike],
) -> None:
    """
    Apply parameters, i.e. replace in file the value of the key with the value of the Parameter
    :param params: List of Parameter objects
    :param src_path: Source path
    :param dst_path: Destination path
    :return: None
    """
    if not isinstance(dst_path, (str, os.PathLike)):
        raise TypeError(
            f"Invalid type: {type(dst_path).__name__}. Expected str or os.PathLike."
        )

    for param in params:
        for ifiles in param.sim_input_files:
            mod_files = os.path.join(dst_path, ifiles)

            with open(mod_files, "r", encoding="utf-8") as in_files:
                content = in_files.read()
            modified_content = content.replace(param.key, str(param.value))
            with open(mod_files, "w", encoding="utf-8") as ou_files:
                ou_files.write(modified_content)
