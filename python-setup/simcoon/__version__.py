from __future__ import annotations

import re
from pathlib import Path

_VERSION_RE = re.compile(r'(?ms)^\[project\].*?^version\s*=\s*"([^"]+)"')
_PYPROJECT_PATH = Path(__file__).resolve().parents[2] / "pyproject.toml"


def _read_version() -> str:
    try:
        pyproject_text = _PYPROJECT_PATH.read_text(encoding="utf-8")
    except OSError:
        return "0+unknown"
    match = _VERSION_RE.search(pyproject_text)
    return match.group(1) if match else "0+unknown"


__version__ = _read_version()
