"""
Data Class to manage Simcoon computation data
"""


from typing import List, Union, Optional
import numpy as np
import numpy.typing as npt

class Data:
    """
    Data class to manage input or output Simcoon computation data
    """

    def __init__(self, list_data: Union[List[npt.NDArray[np.float_]], List[str]]):

        if type(list_data) == List[str]:
            tmp_list = []
            for i in range(len(list_data)):
                tmp_list[i] = np.loadtxt(list_data[i])
            self.list_data = tmp_list
        elif type(list_data) == List[npt.NDArray[np.float_]]:
            self.list_data = list_data
        else:
            raise TypeError("list_data must be a list of numpy arrays or of strings")

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