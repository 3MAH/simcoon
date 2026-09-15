"""
Chaboche Cyclic Plasticity Identification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Identify 7 elasto-plastic Chaboche parameters from 3 cyclic uniaxial tests.

Material model: ``EPCHA`` UMAT — linear elasticity + Voce isotropic hardening
+ two non-linear kinematic backstresses.

==================  =====================  ============================
Symbol              Parameter              Bounds
==================  =====================  ============================
``sigmaY``          initial yield          50 – 300 MPa
``Q``, ``b``        Voce isotropic         100 – 10000 MPa, 0.01 – 10
``C_1``, ``D_1``    1st backstress         1e3 – 1e5 MPa, 10 – 1000
``C_2``, ``D_2``    2nd backstress         1e4 – 1e6 MPa, 10 – 10000
==================  =====================  ============================

Fixed: :math:`E = 140000` MPa, :math:`\\nu = 0.3`, :math:`\\alpha = 10^{-6}`.

The three tests are cyclic strain-controlled tensile experiments at increasing
amplitudes (~1%, ~1.5%, ~2%). Each one needs a **pre-cycling** stage so the
numerical model arrives at the comparison window with realistic accumulated
backstress, then an **initial-state alignment** so it starts at the same
residual strain as the experiment, then a **replay** of the experimental
loading path. This is encoded in three blocks of the ``path_id_N.json`` file:

1. Block 1 (mode 1, linear) — virtual pre-cycle (±1%, ±1.5%, ±2%)
2. Block 2 (mode 1, linear) — set initial residual strain (first row of exp)
3. Block 3 (mode 3, tabular) — replay of the experimental table, read from the
   ``path_id_N_tab1.csv`` file the JSON references

The ``path_id_N.json`` / ``path_id_N_tab1.csv`` pairs are provided in ``data/``
because they are tricky to construct manually (they were converted from the legacy
``path_id_N.txt`` + ``tab_file_N.txt`` pairs with ``scripts/legacy_to_json.py``).
A later step will
replace this scaffolding by Python helpers that build the steps from the experimental
data.

Forward model: :func:`simcoon.solver` (UMAT material-point integrator).
Optimization: :func:`simcoon.identification` (wraps ``differential_evolution``).
Cost: ``nmse_per_response`` — normalises each test's stress column by its own
sum of squares, balancing the three tests despite different stress magnitudes.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import simcoon as sim
from simcoon.parameter import Parameter
from simcoon.identify import identification, calc_cost


# ---------------------------------------------------------------------------
# Test catalogue — file naming mirrors the legacy ``03 - Identification``
# layout (numbering is intentional: 1, 1.5, 2 strain amplitudes).
# ---------------------------------------------------------------------------
TESTS = [
    # name      path file (its mode-3 table embedded)   exp file
    ("test1", "path_id_1.json", "exp_file_1.txt"),
    ("test2", "path_id_2.json", "exp_file_15.txt"),
    ("test3", "path_id_3.json", "exp_file_2.txt"),
]

UMAT_NAME = "EPCHA"
NSTATEV = 33
SOLVER_TYPE = 0
CORATE_TYPE = 2

# Fixed (not identified)
E_FIXED = 140000.0
NU_FIXED = 0.3
ALPHA_FIXED = 1.0e-6

# Identified — order matches the EPCHA props vector after E, nu, alpha.
PARAMS = [
    Parameter(1, bounds=(50,    300),    key="@1p"),  # sigmaY
    Parameter(2, bounds=(100,   10000),  key="@2p"),  # Q
    Parameter(3, bounds=(0.01,  10.0),   key="@3p"),  # b
    Parameter(4, bounds=(1000,  100000), key="@4p"),  # C_1
    Parameter(5, bounds=(10,    1000),   key="@5p"),  # D_1
    Parameter(6, bounds=(10000, 1.0e6),  key="@6p"),  # C_2
    Parameter(7, bounds=(10,    10000),  key="@7p"),  # D_2
]
PARAM_NAMES = ["sigmaY", "Q", "b", "C_1", "D_1", "C_2", "D_2"]

def build_props(x):
    """Assemble the EPCHA props vector from the optimizer's parameter array."""
    return np.array([E_FIXED, NU_FIXED, ALPHA_FIXED, *x])


def run_one_test(props, programme):
    """Run one case and return the predicted σ11 trajectory.

    ``programme`` is the ``(blocks, T_init)`` pair ``load_path_json`` read once in main():
    an identification evaluates this thousands of times, and none of them touches the
    disk or re-parses the path.

    Only the **last** block is returned. The first two blocks are the virtual
    pre-cycle and the initial-state alignment; the experiment corresponds to the
    third one, the mode-3 replay of the table its third block embeds. The legacy result file
    carried exactly that window, so returning the whole history (501 increments
    against 201 experimental points) would break the cost function.
    """
    blocks, T_init = programme
    res = sim.solver.solve(
        blocks, UMAT_NAME, props, NSTATEV, T_init=T_init,
        solver_type=SOLVER_TYPE, corate=CORATE_TYPE,
    )
    block = np.asarray(res["Block"])
    return np.asarray(res["Stress"][0])[block == block.max()]


def cost(x, exp_stresses, programmes):
    """NMSE-per-response cost across the three tests."""
    props = build_props(x)
    y_num = []
    for programme in programmes:
        try:
            sigma11 = run_one_test(props, programme)
        except Exception:
            return 1e12
        y_num.append(sigma11.reshape(-1, 1))
    return calc_cost(exp_stresses, y_num, metric="nmse_per_response")


def main():
    # sim.solver reads/writes relative to cwd
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        script_dir = os.getcwd()
    os.chdir(script_dir)

    path_data = "data"
    path_exp = "exp_data"

    # Experimental σ11 — exp file columns: incr, time, strain, stress
    exp_stresses = []
    for _, _, expfile in TESTS:
        exp = np.loadtxt(os.path.join(path_exp, expfile))
        exp_stresses.append(exp[:, 3].reshape(-1, 1))

    print("=" * 60)
    print(" CHABOCHE CYCLIC PLASTICITY IDENTIFICATION")
    print(" 7 params from 3 cyclic tests, NMSE-per-response cost")
    print("=" * 60)
    for i, (name, pathfile, expfile) in enumerate(TESTS):
        print(f"  {name}: {pathfile} vs {expfile}  ({len(exp_stresses[i])} pts)")

    # The loading programmes, parsed once for the whole identification
    programmes = [sim.solver.load_path_json(os.path.join(path_data, pathfile))[:2]
                  for _, pathfile, _ in TESTS]

    # Gallery budget (~1-2 min). Bump popsize/maxiter for tighter fits.
    result = identification(
        cost, PARAMS,
        args=(exp_stresses, programmes),
        seed=42,
        popsize=15, maxiter=80, tol=1e-6,
        disp=False,
    )

    print()
    print("=" * 60)
    print(" IDENTIFIED PARAMETERS")
    print("=" * 60)
    for n, p in zip(PARAM_NAMES, PARAMS):
        print(f"  {n:8s} = {p.value:>12.3f}    (bounds {p.bounds})")
    print(f"\n  Final cost (NMSE/response) = {result.fun:.4e}")

    # All three tests on one plot — dashed = experiment, solid = identified
    fig, ax = plt.subplots(figsize=(9, 7))
    final_props = build_props(np.array([p.value for p in PARAMS]))
    colors = ["tab:blue", "tab:orange", "tab:green"]
    for (name, _pf, expfile), programme, color in zip(TESTS, programmes, colors):
        exp = np.loadtxt(os.path.join(path_exp, expfile))
        sigma_num = run_one_test(final_props, programme)
        ax.plot(exp[:, 2], exp[:, 3], color=color, linestyle="--",
                linewidth=1.5, label=f"{name} — experiment")
        ax.plot(exp[:, 2], sigma_num, color=color, linestyle="-",
                linewidth=1.5, label=f"{name} — identified")
    ax.set_xlabel(r"strain $\varepsilon_{11}$")
    ax.set_ylabel(r"stress $\sigma_{11}$ [MPa]")
    ax.set_title("Chaboche Cyclic Plasticity — Identified vs Experimental",
                 fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", framealpha=0.9)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
