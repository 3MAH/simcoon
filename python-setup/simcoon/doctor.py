"""Environment doctor: which OpenMP runtimes does this Python process load?

Run it as ``python -m simcoon.doctor``. It imports, one after the other, numpy,
scipy, simcoon and (when installed) torch, records the native libraries each
import brings into the process and reports the OpenMP runtimes found
(``libomp``, ``libgomp``, ``libiomp5``, ``vcomp``). More than one runtime in a
process is the classical source of ``OMP: Error #15`` at import and of random
crashes in threaded kernels; ``KMP_DUPLICATE_LIB_OK=TRUE`` only silences the
guard.

The rule on macOS and Linux with conda: **one environment = one OpenMP
runtime**, conda-forge's ``llvm-openmp``. Build simcoon against the
environment's Armadillo (the default when ``CONDA_PREFIX`` is set) and install
PyTorch from conda-forge rather than from PyPI wheels (which bundle their own
``libomp``).

The report is produced in a subprocess started with ``KMP_DUPLICATE_LIB_OK=TRUE``
so that it survives a duplicate runtime and can describe it.
"""

from __future__ import annotations

import ctypes
import json
import os
import re
import subprocess
import sys
from typing import Dict, List, Optional

# OpenMP runtimes only (not libomptarget / libgomp-plugin-* / libompd helpers)
_OMP_PATTERN = re.compile(
    r"^(libomp\d*(\.[a-z0-9_]+)?|libgomp|libiomp5(md)?|vcomp\d*|libvcomp\d*)(-\d+)?(\.\d+)*\.(dylib|so|dll)(\.\d+)*$",
    re.I)
# linear-algebra stack, matched on the full path (simcoon's own _core, not any package's)
_BLAS_PATTERN = re.compile(
    r"(openblas|libblas|liblapack|mkl_rt|libmkl|Accelerate|armadillo|libsimcoon|[/\\]simcoon[/\\]_core\.)", re.I)

#: Imports probed, in order (missing optional packages are skipped).
DEFAULT_MODULES = ("numpy", "scipy.linalg", "simcoon", "torch")


# ---------------------------------------------------------------------------
# loaded native images of the current process
# ---------------------------------------------------------------------------

def loaded_libraries() -> List[str]:
    """Paths of the native libraries currently loaded in this process."""
    if sys.platform == "darwin":
        libsys = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
        libsys._dyld_image_count.restype = ctypes.c_uint32
        libsys._dyld_get_image_name.restype = ctypes.c_char_p
        libsys._dyld_get_image_name.argtypes = [ctypes.c_uint32]
        n = libsys._dyld_image_count()
        names = (libsys._dyld_get_image_name(i) for i in range(n))
        return [p.decode(errors="replace") for p in names if p]     # None if an image vanished
    if sys.platform.startswith("linux"):
        paths = set()
        try:
            with open("/proc/self/maps") as f:
                for line in f:
                    parts = line.rstrip("\n").split(None, 5)        # the path may contain spaces
                    if len(parts) >= 6 and parts[5].startswith("/"):
                        paths.add(parts[5])
        except OSError:
            pass
        return sorted(paths)
    if sys.platform == "win32":
        from ctypes import wintypes
        psapi = ctypes.WinDLL("psapi")
        kernel32 = ctypes.WinDLL("kernel32")
        kernel32.GetCurrentProcess.restype = wintypes.HANDLE
        psapi.EnumProcessModulesEx.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.HMODULE),
                                               wintypes.DWORD, ctypes.POINTER(wintypes.DWORD), wintypes.DWORD]
        psapi.GetModuleFileNameExW.argtypes = [wintypes.HANDLE, wintypes.HMODULE, wintypes.LPWSTR, wintypes.DWORD]
        h = kernel32.GetCurrentProcess()
        mods = (wintypes.HMODULE * 4096)()
        needed = wintypes.DWORD()
        if not psapi.EnumProcessModulesEx(h, mods, ctypes.sizeof(mods), ctypes.byref(needed), 0x03):
            return []
        count = min(needed.value // ctypes.sizeof(wintypes.HMODULE), 4096)
        buf = ctypes.create_unicode_buffer(1024)
        out = []
        for i in range(count):
            if psapi.GetModuleFileNameExW(h, mods[i], buf, 1024):
                out.append(buf.value)
        return out
    return []


def _probe(modules) -> Dict:
    """Import the modules one by one and attribute newly loaded libraries to each."""
    seen = set(loaded_libraries())
    report = {"python": sys.executable, "platform": sys.platform, "imports": [], "libraries": []}
    for name in modules:
        entry = {"module": name, "status": "ok", "new_libraries": []}
        try:
            __import__(name)
        except ImportError as exc:
            entry["status"] = f"not installed ({exc.__class__.__name__})"
        except Exception as exc:  # noqa: BLE001 - report anything, the doctor must not die
            entry["status"] = f"import failed: {exc!r}"
        now = loaded_libraries()
        new = [p for p in now if p not in seen]
        seen.update(new)
        entry["new_libraries"] = new
        report["imports"].append(entry)
    report["libraries"] = sorted(seen)
    return report


def collect(modules=DEFAULT_MODULES) -> Dict:
    """Run the probe in a fresh subprocess (survives a duplicate OpenMP runtime)."""
    env = dict(os.environ)
    env["KMP_DUPLICATE_LIB_OK"] = "TRUE"        # the child must survive the duplicate it reports
    # Run this file as a script (not `import simcoon.doctor`): importing the simcoon package
    # would load _core and its BLAS/OpenMP libraries before the probe records its baseline.
    proc = subprocess.run([sys.executable, os.path.abspath(__file__), "--probe-json", *modules],
                          capture_output=True, text=True, errors="replace", env=env)
    # the JSON line is enough even if the child then dies at interpreter teardown
    report = None
    for line in reversed(proc.stdout.strip().splitlines()):
        if line.startswith("{"):
            try:
                report = json.loads(line)
            except ValueError:
                report = None
            break
    if report is None:
        return {
            "python": sys.executable, "platform": sys.platform, "imports": [], "libraries": [],
            "error": f"probe process failed (exit {proc.returncode}): {proc.stderr[-2000:]}",
        }
    report["openmp"] = _openmp_runtimes(report)
    omp_paths = {r["path"] for r in report["openmp"]}
    report["blas"] = [p for p in report["libraries"] if _BLAS_PATTERN.search(p) and p not in omp_paths]
    return report


def _openmp_runtimes(report: Dict) -> List[Dict]:
    out = []
    attributed = set()
    for entry in report["imports"]:
        for p in entry["new_libraries"]:
            if _OMP_PATTERN.search(os.path.basename(p)):
                out.append({"path": p, "brought_by": entry["module"]})
                attributed.add(p)
    # safety net: runtimes present before the first probed import
    for p in report["libraries"]:
        if p not in attributed and _OMP_PATTERN.search(os.path.basename(p)):
            out.append({"path": p, "brought_by": "(loaded before the probe)"})
    return out


# ---------------------------------------------------------------------------
# human report
# ---------------------------------------------------------------------------

def _advice(report: Dict) -> List[str]:
    omp = report.get("openmp", [])
    tips = []
    if len(omp) <= 1:
        return tips
    conda = os.environ.get("CONDA_PREFIX")
    for r in omp:
        p, by = r["path"], r["brought_by"]
        if by == "torch" and ("site-packages" in p and "torch" in p):
            tips.append("torch brings its own libomp (PyPI wheel): install PyTorch from conda-forge "
                        "(`conda install -c conda-forge pytorch`) so it shares the environment's runtime.")
        elif by == "simcoon" and conda and not p.startswith(conda):
            tips.append(f"simcoon loads an OpenMP runtime from outside the environment ({p}): rebuild "
                        "simcoon against the environment's Armadillo (`conda install -c conda-forge "
                        "armadillo`, delete build/ or build/cp*, then `pip install -e . --no-build-isolation`).")
        elif conda and not p.startswith(conda):
            tips.append(f"{by} loads {p}, outside the environment: prefer the conda-forge build of {by}.")
    tips.append("Do not rely on KMP_DUPLICATE_LIB_OK=TRUE: it hides the error, the crashes remain.")
    return tips


def format_report(report: Dict) -> str:
    lines = [f"simcoon doctor — {report.get('python')} ({report.get('platform')})"]
    if "error" in report:
        lines.append("  " + report["error"])
        return "\n".join(lines)
    for entry in report["imports"]:
        lines.append(f"  import {entry['module']:<14} {entry['status']}")
    omp = report.get("openmp", [])
    lines.append("")
    lines.append(f"OpenMP runtimes loaded: {len(omp)}")
    for r in omp:
        lines.append(f"  {r['path']}   (brought by {r['brought_by']})")
    if report.get("blas"):
        lines.append("BLAS / linear-algebra libraries:")
        for p in report["blas"]:
            lines.append(f"  {p}")
    lines.append("")
    if len(omp) <= 1:
        lines.append("OK: at most one OpenMP runtime in the process.")
    else:
        lines.append("PROBLEM: several OpenMP runtimes in one process (OMP Error #15 / random crashes).")
        for t in _advice(report):
            lines.append("  - " + t)
    return "\n".join(lines)


def main(argv: Optional[List[str]] = None) -> int:
    import argparse

    argv = sys.argv[1:] if argv is None else list(argv)
    if argv[:1] == ["--probe-json"]:            # child process of collect()
        print(json.dumps(_probe(tuple(argv[1:]))))
        return 0

    parser = argparse.ArgumentParser(prog="python -m simcoon.doctor",
                                     description="Report the OpenMP runtimes loaded by numpy/scipy/simcoon/torch.")
    parser.add_argument("--json", action="store_true", help="machine-readable output")
    parser.add_argument("--modules", nargs="*", default=list(DEFAULT_MODULES),
                        help="modules to import, in order")
    args = parser.parse_args(argv)
    report = collect(args.modules)
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(format_report(report))
    return 0 if len(report.get("openmp", [])) <= 1 and "error" not in report else 1


if __name__ == "__main__":
    sys.exit(main())
