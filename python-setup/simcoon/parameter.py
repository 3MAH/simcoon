import functools
import os
import re
import warnings
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np


class Parameter:
    """
    Parameter class to manage of set of parameters values to be applied during an identification
    and/or DOE test matrix
    :param number: int, number of the parameter
    :param bounds: min and max boundaries in which parameter may vary
    :param key: alphanumeric key to identify the parameter in a file
    :param sim_input_files: simulation input files in which the parameter is used
    """

    def __init__(
        self,
        number: int = 0,
        bounds: Tuple[float, float] = (0.0, 1.0),
        key: str = "",
        sim_input_files: Optional[List[str]] = None,
    ) -> None:
        self.number = number
        self.key = key
        self.bounds = bounds
        self._min_value = bounds[0]
        self._max_value = bounds[1]
        self.sim_input_files = sim_input_files

        self._value = None

    def __repr__(self) -> str:
        val = self.value
        return (
            f"Parameter(number={self.number}, value={val}, "
            f"bounds={self.bounds}, key='{self.key}', "
            f"sim_input_files={self.sim_input_files})"
        )

    @property
    def value(self) -> float:
        if self._value is not None:
            return self._value
        else:
            return 0.5 * (self._min_value + self._max_value)

    @value.setter
    def value(self, value):
        self._value = value


def read_parameters(
    fname: Union[str, os.PathLike] = "data/parameters.inp",
    skiprows: int = 1,
) -> List[Parameter]:
    """
    Read parameters from a simcoon input file.

    :param fname: path to the ``parameters.inp`` file
    :param skiprows: number of leading lines to skip (pandas-style).
                     Defaults to 1 (skip the ``#Number ...`` header).
                     Use ``0`` for files without a header.
    :return: List of Parameter
    """
    if not isinstance(fname, (str, os.PathLike)):
        raise TypeError(
            f"Invalid type: {type(fname).__name__}. Expected str or os.PathLike."
        )
    if skiprows < 0:
        raise ValueError(f"skiprows must be >= 0, got {skiprows}")

    if isinstance(fname, os.PathLike):
        fname = os.fspath(fname)

    params = []
    with open(fname, "r", encoding="utf-8") as paraminit:
        lines = paraminit.readlines()

        for line in lines[skiprows:]:
            values = line.split()
            nfiles = int(values[4])
            param = Parameter(
                number=int(values[0]),
                bounds=(float(values[1]), float(values[2])),
                key=values[3],
                sim_input_files=[values[5 + j] for j in range(nfiles)],
            )
            params.append(param)
    return params


def _input_files(items) -> List[str]:
    """Distinct ``sim_input_files`` of *items*, in order of first appearance."""
    files = []
    for item in items:
        if not all(isinstance(f, str) for f in item.sim_input_files):
            raise TypeError("All elements in sim_input_files must be strings.")
        for f in item.sim_input_files:
            if f not in files:
                files.append(f)
    return files


@functools.lru_cache(maxsize=256)
def _key_pattern(key: str) -> "re.Pattern":
    # Trailing ``(?!\w)`` lookahead prevents ``@E`` from matching ``@Et``.
    # ``re.escape`` keeps user keys with regex metachars literal.
    return re.compile(re.escape(key) + r"(?!\w)")


def substitute_keys(content: str, items, fname: Optional[str] = None) -> str:
    """
    Replace the keys of *items* by their values in the text *content*.

    :param content: text of a template file
    :param items: Parameter or Constant objects (anything with ``key``,
                  ``value`` and ``sim_input_files``)
    :param fname: when given, only the items listing this file in their
                  ``sim_input_files`` are substituted
    :return: the substituted text
    """
    # Longest keys first so e.g. ``@Et`` is replaced before ``@E`` and a
    # shorter key cannot clobber a longer one that contains it.
    ordered = sorted(items, key=lambda it: -len(it.key))
    for item in ordered:
        if fname is not None and fname not in item.sim_input_files:
            continue
        replacement = str(item.value)
        # the ``lambda`` repl avoids backref interpretation in the value
        content = _key_pattern(item.key).sub(lambda _m: replacement, content)
    return content


class KeyTemplate:
    """
    Template files held in memory, written with the keys substituted.

    The templates (e.g. the files of ``keys/``) are read once, at construction.
    Each :meth:`write` substitutes the keys of the parameters, then those of the
    constants, by their current values and writes every working file once. A key
    is replaced only in the files its object lists in ``sim_input_files``.

    :param src_path: directory holding the template files (e.g. ``keys``)
    :param params: Parameter objects whose ``sim_input_files`` name the templates
    :param consts: Constant objects, substituted after the parameters. A
                   constant holds one value per test (``input_values``):
                   ``write(dst, test=i)`` substitutes the value of test *i*.

    Example::

        template = KeyTemplate("keys", params, consts)
        def cost(x):
            for p, v in zip(params, x):
                p.value = v
            for i in range(n_tests):          # e.g. one fibre angle per specimen
                template.write(f"data/test_{i}", test=i)
            ...
    """

    def __init__(
        self,
        src_path: Union[str, os.PathLike],
        params: Sequence = (),
        consts: Sequence = (),
    ) -> None:
        if not isinstance(src_path, (str, os.PathLike)):
            raise TypeError(
                f"Invalid type: {type(src_path).__name__}. Expected str or os.PathLike."
            )
        self.params = list(params)
        self.consts = list(consts)
        self._single_value_warned = set()
        self.templates = {}
        for ifiles in _input_files(self.params + self.consts):
            with open(os.path.join(src_path, ifiles), "r", encoding="utf-8") as f:
                self.templates[ifiles] = f.read()

    def _consts_of_test(self, test: Optional[int]) -> list:
        """Constants with the value of test *test* (``input_values[test]``);
        their own ``value`` when *test* is None."""
        if test is None:
            return self.consts
        consts = []
        for c in self.consts:
            values = np.atleast_1d(c.input_values)
            if values.size == 1:
                if c.key not in self._single_value_warned:
                    self._single_value_warned.add(c.key)
                    warnings.warn(
                        f"constant {c.key} has a single value for several tests: "
                        "it is not a constant of the identification but a fixed "
                        "parameter of the simulation, write it in the template "
                        "instead of a key",
                        stacklevel=3,
                    )
                consts.append(c._replace(value=values[0]))
            elif not -values.size <= test < values.size:
                raise IndexError(
                    f"test {test} out of range for constant {c.key}, which has "
                    f"{values.size} values"
                )
            else:
                consts.append(c._replace(value=values[test]))
        return consts

    def render(self, fname: str, test: Optional[int] = None) -> str:
        """Text of the template *fname* with the current values substituted.

        :param test: index of the test: each constant takes its value
                     ``input_values[test]``; by default, its ``value``
        """
        return self._render(fname, self._consts_of_test(test))

    def _render(self, fname: str, consts: list) -> str:
        content = substitute_keys(self.templates[fname], self.params, fname)
        return substitute_keys(content, consts, fname)

    def write(self, dst_path: Union[str, os.PathLike], test: Optional[int] = None) -> None:
        """Write every template, substituted with the current values, to *dst_path*.

        :param test: index of the test: each constant takes its value
                     ``input_values[test]`` (e.g. the fibre angle of that
                     specimen); by default, its ``value``. A constant with a
                     single value raises a warning: it is the same in every
                     test, so it belongs in the template, not in a key.
        """
        if not isinstance(dst_path, (str, os.PathLike)):
            raise TypeError(
                f"Invalid type: {type(dst_path).__name__}. Expected str or os.PathLike."
            )
        consts = self._consts_of_test(test)
        for fname in self.templates:
            with open(os.path.join(dst_path, fname), "w", encoding="utf-8") as f:
                f.write(self._render(fname, consts))
