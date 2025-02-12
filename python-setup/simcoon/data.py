"""
Data Class to manage computation data
"""

from typing import List, Union, Optional
import numpy as np
import numpy.typing as npt
import glob


class Data:
    """
    Data class to manage input or output Simcoon computation data
    """

    def __init__(
        self,
        control: npt.NDArray[np.float_],
        observation: npt.NDArray[np.float_],
        time: Optional[npt.NDArray[np.float_]] = None,
        timestep: Optional[float] = 0.01,
    ):
        self.control = control
        self.observation = observation
        self.timestep = timestep
        self.increments: npt.NDArray[np.float_] = np.asarray(
            [i + 1 for i in range(np.shape(self.control)[0])]
        )
        if time is None:
            self.time = np.asarray(
                [i * timestep for i in range(np.shape(self.control)[0])]
            )
        else:
            self.time = time


def write_input_and_tab_files(
    list_data: List[Data],
    exp_data_path: str = "exp_data/",
    tab_files_path: str = "data/",
) -> None:
    i = 1
    for element in list_data:
        exp_data_array = np.column_stack(
            (element.increments, element.time, element.control, element.observation)
        )
        tab_file_array = np.column_stack(
            (element.increments, element.time, element.control)
        )
        np.savetxt(exp_data_path + f"input_data_{i:02}.txt", exp_data_array)
        np.savetxt(tab_files_path + "tab_file_" + f"{i:02}" + ".txt", tab_file_array)
        i += 1


def write_files_exp(
    list_data: List[Data], path: str = "data/", exp_data_path: str = "exp_data/"
) -> None:
    list_exp_input_files_names = glob.glob("input_data_*.txt", root_dir=exp_data_path)
    list_nb_columns_in_files = []
    list_nb_observation_columns = []
    list_observation_columns_indices = []
    for element in list_data:
        nb_columns_in_files = (
            element.increments.ndim
            + element.time.ndim
            + element.control.ndim
            + element.observation.ndim
        )
        nb_observation_columns = element.observation.ndim
        observation_columns_indices = [
            i
            for i in range(
                nb_columns_in_files - element.observation.ndim, nb_columns_in_files
            )
        ]
        list_nb_columns_in_files.append(nb_columns_in_files)
        list_nb_observation_columns.append(nb_observation_columns)
        list_observation_columns_indices.append(observation_columns_indices)
    with open(path + "files_exp.inp", "w+") as file:
        file.write("#Name_of_the_exp_files\n")
        for file_name in list_exp_input_files_names:
            file.write(file_name + "\n")
        file.write("\n#EXP_Nb_columns_in_files\n")
        for nb_col_file in list_nb_columns_in_files:
            file.write(str(nb_col_file) + "\n")
        file.write("\n#EXP_Nb_columns_to_identify\n")
        for nb_obs_col in list_nb_observation_columns:
            file.write(str(nb_obs_col) + "\n")
        file.write("\n#EXP_colums_to_identify\n")
        for indices_list in list_observation_columns_indices:
            file.write(" ".join(str(val) for val in indices_list))
            file.write("\n")


def write_files_num(
    list_data: List[Data],
    list_columns_to_compare: List[List[Union[int, str]]],
    path: str = "data/",
) -> None:
    n_columns = 24
    column_header_to_index = {
        "phase": 0,
        "block": 1,
        "step": 2,
        "increment": 3,
        "time": 4,
        "temperature": 5,
        "Q": 6, #specific heat
        "r": 7, #-q
        "e11": 8,
        "e22": 9,
        "e33": 10,
        "e12": 11,
        "e13": 12,
        "e23": 13,
        "s11": 14,
        "s22": 15,
        "s33": 16,
        "s12": 17,
        "s13": 18,
        "s23": 19,
        "Wm": 20, #total deformation energy
        "Wm_r": 21, #reversible deformation energy
        "Wm_ir": 22, #irreversible deformation energy
        "Wm_d": 23, #dissipated deformation energy
    }
    if len(list_data) != len(list_columns_to_compare):
        raise IndexError(
            "list_data and list_columns_to_compare must have the same length"
        )
    with open(path + "files_num.inp", "w+") as file:
        file.write("NUMNb_columnsinfiles\n")
        for _ in list_data:
            file.write(str(n_columns) + "\n")  # total number of columns
        file.write("\nNUMNb_colums_to_identify\n")
        for columns_to_compare in list_columns_to_compare:
            try:
                converted_columns = [
                    column_header_to_index[col] if isinstance(col, str) else col
                    for col in columns_to_compare
                ]
            except KeyError as e:
                raise ValueError(f"Invalid column name: {e.args[0]}")
            file.write(" ".join(str(val) for val in converted_columns))
            file.write("\n")
