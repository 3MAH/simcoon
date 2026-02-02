#pragma once

namespace simcoon_docs {

constexpr auto solver_optimized = R"pbdoc(
    Run the optimized C++ solver for material point simulations.

    This is a high-performance alternative to the Python ``Solver`` class.
    Both solvers accept the same ``Block``/``Step`` objects and return the same
    ``List[HistoryPoint]`` output format, making them interchangeable.

    The C++ solver provides significant performance improvements through:

    - Static UMAT dispatch (function pointer map built once at startup)
    - Pre-allocated Newton-Raphson buffers (reused across all increments)
    - No Python interpreter overhead in the inner loop

    Parameters
    ----------
    blocks : List[Block]
        List of Block objects defining the simulation. Each Block contains:

        - ``steps``: List of Step objects (StepMeca or StepThermomeca)
        - ``umat_name``: Name of the UMAT (e.g., 'ELISO', 'EPICP')
        - ``props``: Material properties array
        - ``nstatev``: Number of state variables
        - ``control_type``: Strain measure ('small_strain', 'green_lagrange', 'logarithmic')
        - ``corate_type``: Objective rate ('jaumann', 'green_naghdi', etc.)
        - ``ncycle``: Number of cycles to repeat the steps

    max_iter : int, optional
        Maximum Newton-Raphson iterations per increment (default: 10).
        Increase for highly nonlinear problems or tight tolerances.

    tol : float, optional
        Convergence tolerance for Newton-Raphson (default: 1e-9).
        The solver converges when the residual norm is below this value.

    lambda_solver : float, optional
        Penalty stiffness for strain-controlled components (default: 10000.0).
        Higher values enforce strain constraints more strictly but may
        cause numerical issues if too large.

    Returns
    -------
    List[HistoryPoint]
        History of state at each converged increment. Each HistoryPoint contains:

        - ``Etot``: Total strain tensor (6,) in Voigt notation
        - ``sigma``: Cauchy stress tensor (6,) in Voigt notation
        - ``Wm``: Mechanical work measures (4,): [Wm, Wm_r, Wm_ir, Wm_d]
        - ``statev``: Internal state variables (nstatev,)
        - ``R``: Rotation matrix (3, 3)
        - ``T``: Temperature

    Raises
    ------
    RuntimeError
        If the UMAT is not found or if Newton-Raphson fails to converge
        after reaching the maximum number of sub-increments.

    See Also
    --------
    simcoon.solver.Solver : Python solver class with identical interface
    simcoon.solver.Block : Block configuration dataclass
    simcoon.solver.StepMeca : Mechanical step configuration
    simcoon.solver.HistoryPoint : Output history point dataclass

    Notes
    -----
    The solver uses automatic sub-incrementation: if Newton-Raphson fails
    to converge, the time step is halved and retried. If it converges easily,
    the time step is doubled (up to the initial increment size).

    The C++ and Python solvers follow the same ``to_start()``/``set_start()``
    architecture from the C++ ``state_variables`` class:

    - ``to_start()``: Reset trial values to start-of-increment (NR rollback)
    - ``set_start(corate_type)``: Save converged values and advance strain

    Examples
    --------
    Basic elastic simulation with uniaxial tension:

    .. code-block:: python

        import numpy as np
        from simcoon.solver import Block, StepMeca
        import simcoon._core as scc

        # Material properties for isotropic elasticity
        E = 70000.0  # Young's modulus (MPa)
        nu = 0.3     # Poisson's ratio
        props = np.array([E, nu])

        # Define uniaxial tension step (strain-controlled in direction 1)
        step = StepMeca(
            DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
            control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
            Dn_init=100
        )

        # Create block with ELISO (isotropic elastic) UMAT
        block = Block(
            steps=[step],
            umat_name='ELISO',
            props=props,
            nstatev=1
        )

        # Run C++ solver
        history = scc.solver_optimized(blocks=[block])

        # Access final stress
        final_stress = history[-1].sigma
        print(f"Final stress: {final_stress[0]:.2f} MPa")

    Plasticity with mixed control (stress-controlled loading):

    .. code-block:: python

        import numpy as np
        from simcoon.solver import Block, StepMeca
        import simcoon._core as scc

        # EPICP (isotropic hardening plasticity) properties
        # [E, nu, sigma_y, H, delta, C, s0, m]
        props = np.array([70000, 0.3, 300, 1000, 0, 0, 0, 0])

        # Stress-controlled loading to 500 MPa
        step = StepMeca(
            Dsigma_end=np.array([500, 0, 0, 0, 0, 0]),
            control=['stress', 'stress', 'stress', 'stress', 'stress', 'stress'],
            Dn_init=200
        )

        block = Block(
            steps=[step],
            umat_name='EPICP',
            props=props,
            nstatev=7
        )

        history = scc.solver_optimized(blocks=[block])

    Compare Python and C++ solvers:

    .. code-block:: python

        from simcoon.solver import Solver, Block, StepMeca
        import simcoon._core as scc
        import numpy as np

        block = Block(
            steps=[StepMeca(DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]), Dn_init=100)],
            umat_name='ELISO',
            props=np.array([70000, 0.3]),
            nstatev=1
        )

        # Python solver
        history_py = Solver(blocks=[block]).solve()

        # C++ solver
        history_cpp = scc.solver_optimized(blocks=[block])

        # Results should match
        assert np.allclose(history_py[-1].sigma, history_cpp[-1].sigma)
)pbdoc";

} // namespace simcoon_docs
