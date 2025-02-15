"""
Data Class to manage computation data
"""
from typing import List, Optional
import numpy as np
import numpy.typing as npt
import glob
np.float_ = np.float64


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


def write_files_exp(list_data: List[Data],
                    exp_data_path: str = "exp_data/",
                    path: str = "data/",
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
            file.write(f"{file_name}\n")
        file.write("\n#EXP_Nb_columns_in_files\n")
        for nb_col_file in list_nb_columns_in_files:
            file.write(f"{nb_col_file}\n")
        file.write("\n#EXP_Nb_columns_to_identify\n")
        for nb_obs_col in list_nb_observation_columns:
            file.write(f"{nb_obs_col}\n")
        file.write("\n#EXP_colums_to_identify\n")
        for indices_list in list_observation_columns_indices:
            file.write(" ".join(str(val) for val in indices_list))
            file.write("\n")
        file.close()
