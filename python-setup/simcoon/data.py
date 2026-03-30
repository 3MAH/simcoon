from typing import List, Optional
import numpy as np
import numpy.typing as npt
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


