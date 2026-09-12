"""Batched Lt_convert must not leak memory.

Each point builds a 6x6 tangent inside libsimcoon that is freed in the _core
extension. The two modules do not share an allocator on Windows (Armadillo's
_aligned_malloc vs carma's numpy allocator), so a silently refused free would
show up here as memory growth proportional to the number of points.
"""

import ctypes
import gc
import os
import sys

import numpy as np
import pytest

import simcoon as sim


def _memory_bytes():
    """Private memory (Windows), resident set (Linux) or peak resident set (macOS)."""
    if sys.platform == "win32":
        class ProcessMemoryCounters(ctypes.Structure):
            _fields_ = [
                ("cb", ctypes.c_ulong),
                ("PageFaultCount", ctypes.c_ulong),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t),
                ("PrivateUsage", ctypes.c_size_t),
            ]

        counters = ProcessMemoryCounters()
        counters.cb = ctypes.sizeof(ProcessMemoryCounters)
        kernel32 = ctypes.WinDLL("kernel32")
        psapi = ctypes.WinDLL("psapi")
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        psapi.GetProcessMemoryInfo.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_ulong]
        if not psapi.GetProcessMemoryInfo(kernel32.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
            pytest.skip("GetProcessMemoryInfo failed")
        return counters.PrivateUsage
    if sys.platform.startswith("linux"):
        with open("/proc/self/statm") as statm:
            return int(statm.read().split()[1]) * os.sysconf("SC_PAGE_SIZE")
    if sys.platform == "darwin":
        import resource

        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss  # bytes on macOS
    pytest.skip("no memory probe on this platform")


def test_batched_lt_convert_does_not_leak():
    n_points, calls = 1000, 100
    rng = np.random.default_rng(0)
    F = (np.tile(np.eye(3)[:, :, None], (1, 1, n_points))
         + 0.02 * rng.standard_normal((3, 3, n_points))).copy(order="F")
    stress = np.asfortranarray(rng.standard_normal((6, n_points)))
    L = np.asarray(sim.L_iso([70000.0, 0.3], "Enu"))
    Lt = np.asfortranarray(np.tile(L[:, :, None], (1, 1, n_points)))

    def fedoo_pair():
        # fedoo's stress_equilibrium: box tangent -> dS/dE, then -> dsigma/dD
        dsde = sim.Lt_convert(Lt, F, stress, "DsigmaDe_2_DSDE")
        sim.Lt_convert(dsde, F, stress, "DSDE_2_Dsigma_logarithmicDD")

    for _ in range(10):
        fedoo_pair()
    gc.collect()
    start = _memory_bytes()
    for _ in range(calls):
        fedoo_pair()
    gc.collect()
    growth = _memory_bytes() - start

    leaked_if_leaking = 2 * calls * n_points * 36 * 8
    assert growth < 0.2 * leaked_if_leaking, (
        f"memory grew by {growth / 2**20:.1f} MB over {calls} batched calls; "
        f"leaking every per-point tangent would add {leaked_if_leaking / 2**20:.0f} MB"
    )
