"""Name maps between the Python solver API and the C++ solver integer codes."""

from simcoon._core import (  # noqa: F401  (re-exported)
    tangent_none,
    tangent_continuum,
    tangent_algorithmic,
    tangent_closest_point,
    tangent_default,
    corate_type_total as _corate_type_total,
)

#: Block types (block.type in the C++ solver)
BLOCK_TYPES = {
    "mechanical": 1,
    "thermomechanical": 2,
}

#: Loading control of a block (control_type in the C++ solver)
CONTROL_TYPES = {
    "small_strain": 1,     # small-strain / Kirchhoff stress control
    "green_lagrange": 2,   # Green-Lagrange strain / PKII stress control
    "logarithmic": 3,      # logarithmic strain / Kirchhoff stress control
    "biot": 4,             # right stretch U (identity = 1, not the Biot strain U - I) / Biot stress control
    "F": 5,                # full deformation gradient control
    "gradU": 6,            # displacement gradient control
}

#: Objective (corotational) rate choice (corate_type in the C++ solver)
CORATE_TYPES = {
    "jaumann": 0,
    "green_naghdi": 1,
    "logarithmic": 2,      # XBM logarithmic rate
    "logarithmic_R": 3,    # log_R rate (A_R strain concentration)
    "truesdell": 4,
    "logarithmic_F": 5,    # log_F rate (A_F strain concentration)
}

#: Finite-strain formulation of the constitutive update. "rate": the law is integrated with an
#: objective rate (CORATE_TYPES). "total": the strain is the state function ln U of F in
#: reference components and the law returns its dual stress (Miehe-Apel-Lambrecht strain space);
#: no rate is involved. The C++ solver carries the total formulation in its corate_type slot.
FORMULATIONS = ("rate", "total")

#: Thermal boundary condition of a thermomechanical step (cBC_T in the C++ solver)
THERMAL_CONTROL = {
    "temperature": 0,      # prescribed temperature ramp
    "heat_flux": 1,        # prescribed heat flux Q
    "convection": 3,       # 0D convection: Q = -q_conv (T - T_init)
}

#: Loading interpolation of a step (mode in the C++ solver)
STEP_MODES = {
    "linear": 1,
    "sinusoidal": 2,
    "tabular": 3,
}

#: Tangent operator modes (tangent_mode carried by the state variables)
TANGENT_MODES = {
    "none": tangent_none,
    "continuum": tangent_continuum,
    "algorithmic": tangent_algorithmic,
    "closest_point": tangent_closest_point,
}


def corate_code(corate, formulation):
    """Code the C++ solver reads in its corate_type slot for a formulation and an objective rate.

    Parameters
    ----------
    corate : str, int or None
        Objective rate (see CORATE_TYPES); None is 'logarithmic_R'. Only under "rate".
    formulation : str
        "rate" or "total" (see FORMULATIONS).
    """
    if formulation not in FORMULATIONS:
        raise ValueError(f"unknown formulation '{formulation}'; valid names: {list(FORMULATIONS)}")
    if formulation == "total":
        if corate is not None:
            raise ValueError("the total formulation involves no objective rate: leave corate unset")
        return _corate_type_total
    return as_code("logarithmic_R" if corate is None else corate, CORATE_TYPES, "corate")


def as_code(value, mapping, what):
    """Convert a name or integer code to the integer code of `mapping`.

    Parameters
    ----------
    value : str or int
        Name (key of `mapping`) or already-valid integer code.
    mapping : dict
        One of the maps of this module.
    what : str
        Human-readable name of the parameter, for error messages.
    """
    if isinstance(value, str):
        try:
            return mapping[value]
        except KeyError:
            raise ValueError(
                f"unknown {what} '{value}'; valid names: {sorted(mapping)}"
            ) from None
    code = int(value)
    if code not in mapping.values():
        raise ValueError(
            f"unknown {what} code {code}; valid codes: {sorted(mapping.values())}"
        )
    return code
