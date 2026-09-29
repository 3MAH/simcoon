"""Parallel batch UMATs must not wait on the GIL.

The worker threads of the batch loops allocate armadillo memory; on Windows libsimcoon routes
it through numpy's allocator, whose tracemalloc hook takes the GIL. If the calling thread kept
the GIL while joining the workers, a traced process would deadlock. Run a batch past the
parallel cutoff under ``-X tracemalloc`` in a subprocess, with a timeout, so a regression fails
instead of hanging the suite. Windows only: elsewhere libsimcoon uses the system allocator.
"""

import os
import subprocess
import sys
import textwrap

import pytest

_SCRIPT = textwrap.dedent("""
    import numpy as np
    import simcoon as sim

    n = 400                                     # past the parallel cutoff (100)
    col = lambda a: np.asfortranarray(np.tile(np.asarray(a, dtype=float).reshape(-1, 1), (1, n)))
    eye = np.asfortranarray(np.tile(np.eye(3)[:, :, None], (1, 1, n)))
    statev = col([290.0, 0, 0, 0, 0, 0, 0, 0])
    stress, _, _, _ = sim.umat("EPICP", col(np.zeros(6)), col([2e-3, -6e-4, -6e-4, 0, 0, 0]),
                               eye, eye, col(np.zeros(6)), eye,
                               col([200000.0, 0.3, 0.0, 300.0, 1000.0, 0.5]), statev,
                               0.5, 1.0, col(np.zeros(4)), n_threads=4)
    assert np.isfinite(stress).all()
    stress_T = sim.umat_T("ELISO", col(np.zeros(6)), col([1e-3, 0, 0, 0, 0, 0]), col(np.zeros(6)),
                          eye, col([7800.0, 460.0, 200000.0, 0.3, 1e-5]), col([290.0]),
                          0.5, 1.0, col(np.zeros(4)), col(np.zeros(3)),
                          np.full(n, 290.0), np.zeros(n), n_threads=4)[0]
    assert np.isfinite(stress_T).all()
    quats = np.tile([0.0, 0.0, np.sin(0.3), np.cos(0.3)], (n, 1))
    for rotate in (sim._core._batch_voigt_stress_rotation, sim._core._batch_voigt_strain_rotation):
        assert np.isfinite(rotate(quats)).all()
    print("ok")
""")


@pytest.mark.skipif(sys.platform != "win32",
                    reason="only Windows routes libsimcoon's kernel allocations through numpy")
def test_batch_umat_does_not_wait_on_the_gil():
    cmd = [sys.executable, "-X", "tracemalloc", "-c", _SCRIPT]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=120, env=os.environ.copy())
    except subprocess.TimeoutExpired:
        pytest.fail("batch sim.umat hung: the parallel loop waited on the GIL")
    assert out.returncode == 0, out.stderr[-2000:]
    assert out.stdout.strip().endswith("ok")
