"""
Data Class to manage Simcoon computation data
"""


from typing import List, Union
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

    def write_input_file(self):
        pass

    def write_tab_file(self):
        pass