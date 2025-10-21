from typing import List, Optional

class Problem(ProblemBase):
    """Base class to define a problem that generate a linear system and to solve
    the linear system with some defined boundary conditions.

    The linear problem is written under the form:
    A*X = B+D
    where:
     * A is a square matrix build with the associated assembly object calling
         assembly.get_global_matrix()
     * X is the column vector containing the degrees of freedom (solution after solving)
     * B is a column vector used to set Neumann boundary conditions
     * D is a column vector build with the associated assembly object calling
         assembly.get_global_vector()

    Parameters
    ----------
    A: scipy sparse matrix
        Matrix that define the discretized linear system to solve.
    B: np.ndarray or 0
        if 0, B is initialized to a zeros array with de adequat shape.
    D: np.ndarray or 0
        if 0, D is ignored.
    mesh: fedoo Mesh
        mesh associated to the problem.
    name: str (default = "MainProblem")
        name of the problem.
    space: ModelingSpace(Optional)
        ModelingSpace on which the problem is defined.
    name: str
        name of the problem.
    """