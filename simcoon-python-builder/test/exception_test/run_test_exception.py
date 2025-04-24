"""Test numpy array to matrix conversion function."""

import numpy as np
import numpy.typing as npt
import test_exception as exception
import pytest


@pytest.fixture(scope="session")
def Mat_fail_eigen_sym() -> npt.NDArray[np.float64]:
    mat_fail_eigen_sym = np.array(
        [
            [1.0e1222, 2, 3],
            [9, 1, 4],  # 9 ≠ 2 → asymmetric in lower triangle
            [0, 0, 1],
        ]
    )
    return mat_fail_eigen_sym


def test_eig_sym(Mat_fail_eigen_sym):
    """Test eigen_sym exception."""

    try:
        eigval = exception.test_eig_sym_val_affect(Mat_fail_eigen_sym, True)
    except exception.CppExceptionEigSym as e:
        print("Caught C++ exception in test_eig_sym_val:", e)

    try:
        eigval = exception.test_eig_sym_val_affect(Mat_fail_eigen_sym, True)
    except RuntimeError as e:
        print("Caught C++ std::runtime_error exception in test_eig_sym_val:", e)

    try:
        eigval = exception.test_eig_sym_val_modify(Mat_fail_eigen_sym, True)
    except exception.CppExceptionEigSym as e:
        print("Caught C++ exception in test_eig_sym_val_vec:", e)

    try:
        (eigval, eigvec) = exception.test_eig_sym_val_vec(Mat_fail_eigen_sym, True)
    except exception.CppExceptionEigSym as e:
        print("Caught C++ exception in test_eig_sym_val_vec:", e)
