from __future__ import annotations

from pathlib import Path

_PYPROJECT_PATH = Path(__file__).resolve().parents[2] / "pyproject.toml"


def _read_version() -> str:
    try:
        pyproject_text = _PYPROJECT_PATH.read_text(encoding="utf-8")
    except OSError:
        return "0+unknown"

    in_project_section = False
    for raw_line in pyproject_text.splitlines():
        line = raw_line.strip()

        if line.startswith("[") and line.endswith("]"):
            in_project_section = line == "[project]"
            continue

        if in_project_section and line.startswith("version ="):
            _, value = line.split("=", 1)
            return value.strip().strip('"').strip("'")

    return "0+unknown"


__version__ = _read_version()
