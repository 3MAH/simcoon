"""
Data Class to manage computation data
"""


from typing import List, Union, Optional
import numpy as np
import numpy.typing as npt


class Data:
    """
    Data class to manage input or output Simcoon computation data
    """

    def __init__(self, control: npt.NDArray[np.float_],
                 observation: npt.NDArray[np.float_],
                 time: Optional[npt.NDArray[np.float_]],
                 timestep: Optional[np.float] = 0.01):

        self.control = control
        self.observation = observation
        self.timestep = timestep
        self.increment = np.array([i+1 for i in range(np.shape(self.control)[0])])
        if time is None:
            self.time = np.array([i*timestep for i in range(np.shape(self.control)[0])])
        else:
            self.time = time


def write_input_and_tab_files(self, exp_data_path: str = "exp_data/", tab_files_path: str = "data/") -> None:
    for element in self.list_data:  # order: strain, stress
        i = 1
        increments = np.array([j + 1 for j in range(len(element))])
        time = np.array([j * 0.01 for j in range(len(element))])
        strain = element[:, 0]
        exp_data_array = np.column_stack((increments, time, element))
        tab_file_array = np.column_stack((increments, time, strain))
        np.savetxt(exp_data_path + "input_data_" + str(i) + ".txt", exp_data_array)
        np.savetxt(tab_files_path + "tab_file_" + str(i) + ".txt", tab_file_array)
        i += 1


def write_files_exp(self, list_nb_columns_in_files: List[int], list_nb_columns_to_identify: List[int], list_columns_to_identify: List[int], path: str = "data/", ) -> None:
    pass


def write_files_num(self, list_nb_columns_in_files: List[int], list_columns_to_identify: List[int], path: str = "data/", ) -> None:
    pass