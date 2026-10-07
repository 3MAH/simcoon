"""
Identification through keys in a JSON material file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The key system of :mod:`simcoon.parameter` and :mod:`simcoon.constant` is a plain
text substitution: a template file holds placeholders (``@sigmaY``, ``@k``, ``@m``),
and at every evaluation the optimizer's values replace the parameter keys in a
working copy, the constant keys taking their fixed value. Since 2.0 the material of the
solver is a JSON file, so the template is one too. Note that ``keys/material.json``
is **not** valid JSON before the substitution (the placeholders are bare, not quoted);
it becomes valid once the numbers are in, which is the whole point of a template.

The chain at each evaluation:

1. ``copy_parameters`` copies ``keys/material.json`` to the working directory,
2. ``apply_parameters`` and ``apply_constants`` write the values in place of the keys,
3. ``load_simulation_json`` reads the material back with the loading path,
4. ``sim.solver.solve`` runs the tension test in memory,
5. ``calc_cost`` compares the stress history with the experiment.

``sim.identification`` wraps ``scipy.optimize.differential_evolution`` with the
bounds of the :class:`~simcoon.parameter.Parameter` objects read from
``data/parameters.inp``. The "experiment" is synthetic: the same model at the
reference values, with noise added, so the identified values can be checked.

Forward model: ``EPICP`` (isotropic power-law hardening
:math:`\\sigma_Y + k\\,p^m`), identified parameters ``sigmaY`` and ``k``; the
exponent ``m`` is a :class:`~simcoon.constant.Constant` read from
``data/constants.inp`` (on a monotonic tension test ``k`` and ``m`` compensate each
other, so the exponent is fixed from another source); ``E``, ``nu`` and ``alpha``
are written in the template.
"""

import os
import shutil
import tempfile

import numpy as np
import matplotlib.pyplot as plt

import simcoon as sim
from simcoon.constant import apply_constants, read_constants
from simcoon.identify import calc_cost, identification
from simcoon.parameter import apply_parameters, copy_parameters, read_parameters
from simcoon.solver import Block, StepMeca, load_simulation_json, save_path_json

UMAT_NAME = "EPICP"
REFERENCE = {"@sigmaY": 300.0, "@k": 1000.0}   # the values to recover (m = 0.3 is a constant)


def build_working_dir(script_dir):
    """A scratch directory holding the loading path; the material comes from the keys."""
    work = tempfile.mkdtemp(prefix="simcoon_keys_")
    step = StepMeca(control=["strain"] + ["stress"] * 5, value=[0.03, 0, 0, 0, 0, 0],
                    time=1.0, ninc=60)
    save_path_json(os.path.join(work, "path.json"), [Block(steps=[step])], T_init=293.15)
    return work


def run_from_files(params, consts, keys_dir, work):
    """Substitute the keys, read the JSON pair back, run the test: sigma_11(t)."""
    copy_parameters(params, src_path=keys_dir, dst_path=work)
    apply_parameters(params, dst_path=work)
    apply_constants(consts, dst_path=work)
    kwargs = load_simulation_json(os.path.join(work, "material.json"),
                                  os.path.join(work, "path.json"))
    res = sim.solver.solve(**kwargs)
    return np.asarray(res["Strain"][0]), np.asarray(res["Stress"][0])


def main():
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:   # executed by the docs gallery, from the example's directory
        script_dir = os.getcwd()
    keys_dir = os.path.join(script_dir, "keys")
    params = read_parameters(os.path.join(script_dir, "data", "parameters.inp"))
    consts = read_constants(1, os.path.join(script_dir, "data", "constants.inp"))
    work = build_working_dir(script_dir)

    # Synthetic experiment: the reference values, plus 0.5 % multiplicative noise
    for p in params:
        p.value = REFERENCE[p.key]
    strain, sigma_ref = run_from_files(params, consts, keys_dir, work)
    rng = np.random.default_rng(0)
    sigma_exp = sigma_ref * (1.0 + 0.005 * rng.standard_normal(sigma_ref.size))
    y_exp = [sigma_exp.reshape(-1, 1)]

    def cost(x):
        for p, value in zip(params, x):
            p.value = value
        try:
            _, sigma = run_from_files(params, consts, keys_dir, work)
        except Exception:
            return 1e12
        return calc_cost(y_exp, [sigma.reshape(-1, 1)], metric="nmse_per_response")

    result = identification(cost, params, seed=1, maxiter=30, popsize=10, tol=1e-8,
                            polish=True, disp=False)
    print(f"differential_evolution: {result.nfev} evaluations, cost {result.fun:.3e}")
    for p in params:
        print(f"  {p.key:8s} = {p.value:10.3f}   (reference {REFERENCE[p.key]:g})")

    # The material file the last evaluation left behind is the identified one
    _, sigma_id = run_from_files(params, consts, keys_dir, work)
    with open(os.path.join(work, "material.json")) as f:
        print("identified material.json:\n" + f.read())
    shutil.rmtree(work, ignore_errors=True)

    plt.figure(figsize=(6, 4))
    plt.plot(100 * strain, sigma_exp, "k.", ms=4, label="synthetic experiment")
    plt.plot(100 * strain, sigma_ref, "k--", lw=1, label="reference")
    plt.plot(100 * strain, sigma_id, "r-", lw=1.5, label="identified")
    plt.xlabel("strain (%)")
    plt.ylabel("stress (MPa)")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
