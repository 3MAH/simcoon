#pragma once

namespace simcoon_docs {

constexpr auto stress_convert = R"pbdoc(
Convert between different stress measures using a deformation gradient.

This function provides a single, flexible Python wrapper that converts
between common stress measures used in continuum mechanics. The
conversion is performed using the provided deformation gradient `F` and
an optional Jacobian `J` when required by the transformation.

Parameters
----------
sigma : pybind11::array_t<double>
        The input stress. Accepted shapes:
        - 1D array with 6 elements: Voigt vector [s11, s22, s33, s23, s13, s12].
        - 2D array with shape (3,3): full 3x3 stress matrix.
        - 2D array with shape (6,N): list of N Voigt stress vectors (columns).
F : pybind11::array_t<double>
        The deformation gradient. Accepted shapes:
        - (3,3) for a single conversion.
        - (3,3,N) for N conversions when `sigma` is provided as (6,N) or (3,3,N).
converter_key : str
        A string key selecting the conversion to perform. Supported keys are:

        - 'Cauchy2PKI'   : Cauchy -> First Piola-Kirchhoff (P)
        - 'Cauchy2PKII'  : Cauchy -> Second Piola-Kirchhoff (S)
        - 'Cauchy2Kirchoff' : Cauchy -> Kirchoff (tau)
        - 'Kirchoff2Cauchy' : Kirchoff -> Cauchy
        - 'Kirchoff2PKI'  : Kirchoff -> First Piola-Kirchhoff
        - 'Kirchoff2PKII' : Kirchoff -> Second Piola-Kirchhoff
        - 'PKI2Kirchoff'  : First Piola-Kirchhoff -> Kirchoff
        - 'PKII2Kirchoff' : Second Piola-Kirchhoff -> Kirchoff
        - 'PKI2Cauchy'   : First Piola-Kirchhoff -> Cauchy
        - 'PKII2Cauchy'  : Second Piola-Kirchhoff -> Cauchy

J : float, optional
        The determinant of `F` (Jacobian). If not provided (or zero) the
        function may compute it from `F` when required by the chosen
        conversion. Default is 0.0.
copy : bool, optional
        If true, a copy of the returned numpy array is made. Default is true.

Returns
-------
pybind11::array_t<double>
        The converted stress. The return shape mirrors the input `sigma`:
        - If `sigma` was a 1D Voigt vector -> returns a 1D Voigt vector (6,).
        - If `sigma` was a 3x3 matrix -> returns a 3x3 matrix.
        - If `sigma` was (6,N) or (3,3,N) -> returns the same shape with N converted columns/slices.

Notes
-----
- Voigt ordering used throughout is [11, 22, 33, 23, 13, 12].
- Shear components follow the engineering convention when converting
        between Voigt and full 3x3 representations.
- The wrapper dispatches to the underlying C++ implementations
        (e.g., `Cauchy2PKI`, `PKII2Cauchy`, etc.).

Shape and orientation notes:

When providing batches, `sigma` must be shaped (6, N) and `F` must be shaped
(3, 3, N). The implementation indexes `F.shape[2]` for N. Do not provide `F`
as (N, 3, 3) — the wrapper expects (3, 3, N). Single conversions expect `sigma`
shape (6,) and `F` shape (3, 3).

J handling:

If `J` is left as the default 0.0, the C++ conversion routines may compute
det(F) internally when required. You may pass a precomputed `J` for
performance or to override the computed determinant.

Converter key lookup and errors:

Valid converter keys are listed above. Passing an unknown key may currently
cause a lookup that inserts a default entry in the map and can lead to
undefined `select` behaviour. It's recommended to use one of the documented
keys. For stricter behaviour, the implementation can be changed to use
`map::at()` to throw on unknown keys.

Exceptions and errors:

The function throws `std::invalid_argument` for invalid input shapes or
dimensions (e.g., wrong sigma length, mismatched batch sizes). Other
exceptions may be thrown by the underlying C++ converters.

Examples
--------
Single 6-element Voigt -> first Piola-Kirchhoff:

.. code-block:: python

        import numpy as np
        import simcoon as sim

        sigma = np.array([100.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        F = np.eye(3)
        P = sim.stress_convert(sigma, F, 'Cauchy2PKI')
        print(P)

Batch conversion of N Voigt vectors to Cauchy using matching F slices:

.. code-block:: python

        import numpy as np
        import simcoon as sim

        N = 5
        sigma_list = np.zeros((6, N))
        F_list = np.repeat(np.eye(3)[:, :, np.newaxis], N, axis=2)
        cauchy = sim.stress_convert(sigma_list, F_list, 'PKII2Cauchy')
        print(cauchy.shape)  # (6, N)
)pbdoc";

} // namespace simcoon_docs
