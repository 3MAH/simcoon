from typing import List, Optional


class Step:
    def __init__(self, control_type: int):
        self.control_type = control_type

    def __str__(self):
        return f"Step with control type {self.control_type}"


class StepMeca(Step):
    def __init__(self, control_type: int):
        super().__init__(control_type)

    def __str__(self):
        return f"Mechanical Step with control type {self.control_type}"


class StepThermoMeca(Step):
    def __init__(self, control_type: int):
        super().__init__(control_type)

    def __str__(self):
        return f"Thermo-Mechanical Step with control type {self.control_type}"


class Block:
    def __init__(
        self,
        number: int = 0,
        nstep: int = 0,
        ncycle: int = 0,
        type_: str = "mechanical",
        control_type: str = "linear",
        steps: Optional[List[Step]] = None,
    ):
        assert nstep >= 0 and ncycle >= 0, "Steps and cycles must be non-negative."
        self.number = number
        self.nstep = nstep
        self.ncycle = ncycle
        self.type = type_
        self.control_type = control_type
        self.steps = steps if steps is not None else []

    def generate(self):
        assert (
            self.number >= 0
            and self.nstep > 0
            and self.ncycle > 0
            and self.type in ["mechanical", "thermo-mechanical", "thermal"]
            and self.control_type in ["linear", "sinusoidal"]
        ), "Invalid block parameters."

        if self.type == "mechanical":
            self.steps = [StepMeca(self.control_type) for _ in range(self.nstep)]
        elif self.type == "thermo-mechanical":
            self.steps = [StepThermoMeca(self.control_type) for _ in range(self.nstep)]
        else:
            print(
                f"Please provide a consistent loading type for the block {self.number}"
            )

    def __str__(self):
        output = [
            f"Display info on the block {self.number}",
            f"Number of steps: {self.nstep}",
            f"Number of cycles: {self.ncycle}",
            f"Type of block: {self.type}",
            f"Control type of block: {self.control_type}",
        ]

        for step in self.steps:
            output.append(str(step))

        return "\n".join(output)

    def __eq__(self, other):
        if not isinstance(other, Block):
            return False
        return (
            self.number == other.number
            and self.nstep == other.nstep
            and self.ncycle == other.ncycle
            and self.type == other.type
            and self.control_type == other.control_type
            and self.steps == other.steps
        )
