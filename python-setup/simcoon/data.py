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

    def write_input_files(self, path: str = "exp_data/") -> None:
        for element in self.list_data:  # order: strain, stress
            i = 1
            increments = np.array([j + 1 for j in range(len(element))])
            time = np.array([j * 0.01 for j in range(len(element))])
            data_array = np.column_stack((increments, time, element))
            np.savetxt(path + "input_data_" + str(i) + ".txt", data_array)
            i += 1


    def write_tab_files(self, path: str = "data/") -> None:
        for element in self.list_data:
            i = 1
            increments = np.array([j + 1 for j in range(len(element))])
            time = np.array([j * 0.01 for j in range(len(element))])
            strain = element[:, 0]
            data_array = np.column_stack((increments, time, strain))
            np.savetxt(path + "tab_file_" + str(i) + ".txt", data_array)
            i += 1
