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

#: Kinematic form of the finite-strain constitutive update: ONE value, which names the form and,
#: for the rate form, its objective rate. "total": the strain is the state function ln U of F in
#: reference components and the law returns its dual stress (Miehe-Apel-Lambrecht strain space);
#: no rate is involved. A name of CORATE_TYPES: the rate form with that objective rate. "rate":
#: the rate form with its default rate, 'logarithmic_R'. The C++ solver carries the total form
#: in its corate_type slot.
KINEMATIC_FORMS = ("total", "rate", *CORATE_TYPES)
_KINEMATIC_CODES = {"total": _corate_type_total, "rate": CORATE_TYPES["logarithmic_R"], **CORATE_TYPES}

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


def kinematic_code(kinematic_form=None, corate=None):
    """Code the C++ solver reads in its corate_type slot for a kinematic form.

    ``kinematic_form`` and ``corate`` are two names of the same value (``corate`` is the former
    one): "total", "rate", or an objective rate by name or code (see KINEMATIC_FORMS). Neither
    given: 'logarithmic_R'. Both given: they must agree, ``kinematic_form="rate"`` leaving the
    choice of the rate to ``corate``.
    """
    def code(value):
        if isinstance(value, str):
            return as_code(value, _KINEMATIC_CODES, "kinematic_form")
        return as_code(value, CORATE_TYPES, "kinematic_form")   # the total form has a name only

    if kinematic_form is None:
        return code("rate" if corate is None else corate)
    if corate is None or code(corate) == code(kinematic_form):
        return code(kinematic_form)
    if kinematic_form == "rate" and code(corate) != _corate_type_total:
        return code(corate)
    raise ValueError(f"kinematic_form and corate (its former name) are one value; got "
                     f"{kinematic_form!r} and {corate!r}")


def kinematic_name(kinematic_form=None, corate=None):
    """The documented name of a kinematic form: "total" or the name of its objective rate."""
    code = kinematic_code(kinematic_form, corate)
    return next(name for name, c in {"total": _corate_type_total, **CORATE_TYPES}.items() if c == code)


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
