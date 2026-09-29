"""Input buffers stay owned by NumPy, even with an incompatible custom allocator.

Requires a C compiler and setuptools (the dev extra). Each scenario runs in a
subprocess because an allocator mismatch can crash rather than raise Python errors.
"""

import gc
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


@pytest.fixture(scope="module")
def allocator_probe(tmp_path_factory):
    from setuptools import Distribution, Extension
    from setuptools.command.build_ext import build_ext

    build = tmp_path_factory.mktemp("numpy_allocator")
    source = Path(__file__).with_name("numpy_allocator_probe.c")
    distribution = Distribution({"ext_modules": [Extension(
        "_numpy_allocator_probe", [str(source)], include_dirs=[np.get_include()],
    )]})
    command = build_ext(distribution)
    command.build_lib = str(build)
    command.build_temp = str(build / "temp")
    command.ensure_finalized()
    command.run()
    return build


def _layout(a, kind):
    a = np.asarray(a, dtype=float)
    if kind in ("C", "F"):
        return np.array(a, order=kind, copy=True)
    if kind == "strided":
        backing = np.empty(a.shape[:-1] + (2 * a.shape[-1],))
        view = backing[..., ::-2]
        view[...] = a
        return view
    if kind == "unaligned":
        backing = np.empty(a.nbytes + 1, dtype=np.uint8)
        view = np.ndarray(a.shape, dtype=float, buffer=backing, offset=1, order="F")
        view[...] = a
        return view
    result = np.array(a, order="F", copy=True)
    result.flags.writeable = False
    return result


def _exercise(kind):
    import simcoon as sim

    n = 128  # exercise workers as well as serial conversion
    col = lambda a: _layout(np.tile(np.asarray(a).reshape(-1, 1), (1, n)), kind)
    eye = _layout(np.repeat(np.eye(3)[:, :, None], n, axis=2), kind)
    inputs = [
        col(np.zeros(6)), col([2e-3, -6e-4, -6e-4, 0, 0, 0]),
        eye, eye, col(np.zeros(6)), eye,
        col([200000., 0.3, 0., 300., 1000., 0.5]),
        col([290., 0, 0, 0, 0, 0, 0, 0]),
    ]
    wm = col(np.zeros(4))
    originals = [a.copy() for a in inputs + [wm]]
    metadata = [(a.shape, a.strides, a.flags.owndata, a.flags.writeable)
                for a in inputs + [wm]]
    outputs = sim.umat("EPICP", *inputs, 0.5, 1., wm, n_threads=4)
    for a, original, meta in zip(inputs + [wm], originals, metadata):
        np.testing.assert_array_equal(a, original)
        assert (a.shape, a.strides, a.flags.owndata, a.flags.writeable) == meta
    # Small vector copy; matrix/cube view fallbacks; zero-copy output capsules.
    elastic = sim.L_iso(_layout([70000., 0.3], kind), "Enu")
    tangent = _layout(np.repeat(elastic[:, :, None], n, axis=2), kind)
    converted = sim.Lt_convert(tangent, eye, col(np.zeros(6)), "DsigmaDe_2_DSDE")
    strain = sim.Log_strain(eye)
    return (*outputs, elastic, converted, strain)


@pytest.mark.parametrize("layout", ["F", "C", "strided", "unaligned", "readonly"])
def test_numpy_allocator_ownership(allocator_probe, layout):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(allocator_probe) + os.pathsep + env.get("PYTHONPATH", "")
    completed = subprocess.run(
        [sys.executable, "-X", "tracemalloc", str(Path(__file__).resolve()), layout],
        env=env, capture_output=True, text=True, timeout=120,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "allocator ownership OK" in completed.stdout


if __name__ == "__main__":
    import _numpy_allocator_probe as probe

    # Avoid an OS crash dialog in the subprocess if an old build frees a shifted pointer.
    if sys.platform == "win32":
        import ctypes
        ctypes.windll.kernel32.SetErrorMode(0x0001 | 0x0002)
    expected = _exercise(sys.argv[1])  # warm imports before installing the handler
    before = probe.total()
    previous = probe.install()
    try:
        actual = _exercise(sys.argv[1])
    finally:
        probe.restore(previous)
    assert probe.total() > before, "NumPy did not use the test allocator"
    # Inputs have died and the handler has changed; outputs must still be usable.
    for got, want in zip(actual, expected):
        np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-10)
    del actual, got
    gc.collect()
    assert probe.count() == 0, f"{probe.count()} NumPy allocations lost their owner"
    print("allocator ownership OK")
