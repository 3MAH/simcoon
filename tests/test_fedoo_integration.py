"""Integration tests for fedoo + simcoon.

Tests both non-linear umat (OpenMP) and non-linear geometry
(tangent modulus transfer via deformation gradient).
"""

import numpy as np
import pytest

import fedoo
import simcoon
from simcoon import simmit as sim


@pytest.fixture(autouse=True)
def _reset_fedoo():
    """Reset fedoo global state between tests."""
    fedoo.ModelingSpace("3D")
    yield


def _make_cube(name, nx=2):
    """Create a unit cube mesh with standard Dirichlet BCs on min faces."""
    return fedoo.mesh.box_mesh(
        nx=nx, ny=nx, nz=nx,
        x_min=0, x_max=1, y_min=0, y_max=1, z_min=0, z_max=1,
        elm_type="hex8", name=name,
    )


def _apply_compression_bc(pb, mesh, disp_z):
    """Fix min faces, prescribe DispZ on top face."""
    pb.bc.add("Dirichlet", mesh.find_nodes("X", mesh.bounding_box.xmin), "DispX", 0)
    pb.bc.add("Dirichlet", mesh.find_nodes("Y", mesh.bounding_box.ymin), "DispY", 0)
    pb.bc.add("Dirichlet", mesh.find_nodes("Z", mesh.bounding_box.zmin), "DispZ", 0)
    pb.bc.add("Dirichlet", mesh.find_nodes("Z", mesh.bounding_box.zmax), "DispZ", disp_z)


def test_epicp_small_strain():
    """Elastoplastic EPICP umat, small strain (exercises OpenMP)."""
    mesh = _make_cube("epicp_ss")
    props = np.array([200e3, 0.3, 0.0, 300.0, 1000.0, 0.5])
    mat = fedoo.constitutivelaw.Simcoon("EPICP", props, name="epicp_ss_mat")

    wf = fedoo.weakform.StressEquilibrium(mat, name="wf_epicp_ss")
    asm = fedoo.Assembly.create(wf, mesh, elm_type="hex8", name="asm_epicp_ss")

    pb = fedoo.problem.NonLinear(asm, name="pb_epicp_ss")
    _apply_compression_bc(pb, mesh, disp_z=0.01)

    pb.nlsolve(dt=0.5, tmax=1.0, update_dt=True, print_info=0)

    disp = pb.get_disp("DispZ")
    assert np.max(np.abs(disp)) == pytest.approx(0.01, abs=1e-8)
    p = asm.sv["Statev"][1]
    assert np.max(p) > 0, "Expected plastic deformation"


def test_neohookean_nlgeom():
    """Neo-Hookean hyperelastic, large strain (tangent modulus transfer)."""
    mesh = _make_cube("neohc_nl")
    props = np.array([80e3, 400e3])  # mu, kappa
    mat = fedoo.constitutivelaw.Simcoon("NEOHC", props, name="neohc_nl_mat")

    wf = fedoo.weakform.StressEquilibrium(mat, name="wf_neohc_nl", nlgeom=True)
    asm = fedoo.Assembly.create(wf, mesh, elm_type="hex8", name="asm_neohc_nl")

    pb = fedoo.problem.NonLinear(asm, nlgeom=True, name="pb_neohc_nl")
    _apply_compression_bc(pb, mesh, disp_z=0.1)

    pb.nlsolve(dt=0.5, tmax=1.0, update_dt=True, print_info=0)

    disp = pb.get_disp("DispZ")
    assert np.max(np.abs(disp)) == pytest.approx(0.1, abs=1e-8)
    F = asm.sv["F"]
    assert not np.allclose(F, np.eye(3).reshape(3, 3, 1))


def test_epicp_nlgeom():
    """EPICP + large strain (OpenMP umat + tangent modulus + F transfer)."""
    mesh = _make_cube("epicp_nl")
    props = np.array([200e3, 0.3, 0.0, 300.0, 1000.0, 0.5])
    mat = fedoo.constitutivelaw.Simcoon("EPICP", props, name="epicp_nl_mat")

    wf = fedoo.weakform.StressEquilibrium(mat, name="wf_epicp_nl", nlgeom=True)
    asm = fedoo.Assembly.create(wf, mesh, elm_type="hex8", name="asm_epicp_nl")

    pb = fedoo.problem.NonLinear(asm, nlgeom=True, name="pb_epicp_nl")
    _apply_compression_bc(pb, mesh, disp_z=0.05)

    pb.nlsolve(dt=0.25, tmax=1.0, update_dt=True, print_info=0)

    disp = pb.get_disp("DispZ")
    assert np.max(np.abs(disp)) == pytest.approx(0.05, abs=1e-8)
    p = asm.sv["Statev"][1]
    assert np.max(p) > 0, "Expected plastic deformation"
    F = asm.sv["F"]
    assert not np.allclose(F, np.eye(3).reshape(3, 3, 1))
    Lt = asm.sv["TangentMatrix"]
    sym_err = np.abs(Lt - Lt.transpose((1, 0, 2))).max()
    assert sym_err < 1e-3, f"Tangent matrix not symmetric (err={sym_err})"


@pytest.mark.parametrize("n_threads", [1, 2, 4])
def test_umat_openmp_threads(n_threads):
    """Direct sim.umat with explicit n_threads to verify OpenMP."""
    n_gauss = 64
    props = np.asfortranarray(
        np.tile([200e3, 0.3, 0.0, 300.0, 1000.0, 0.5], (n_gauss, 1)).T
    )
    etot = np.zeros((6, n_gauss), order="F")
    Detot = np.zeros((6, n_gauss), order="F")
    Detot[2, :] = 5e-3

    sigma = np.zeros((6, n_gauss), order="F")
    statev = np.zeros((8, n_gauss), order="F")
    Wm = np.zeros((4, n_gauss), order="F")
    DR = np.empty((3, 3, n_gauss), order="F")
    DR[...] = np.eye(3).reshape(3, 3, 1)

    stress, sv, wm, Lt = sim.umat(
        "EPICP", etot, Detot,
        np.array([]), np.array([]),
        sigma, DR, props, statev,
        0.0, 1.0, Wm, n_threads=n_threads,
    )

    assert np.allclose(stress[:, 0:1], stress, atol=1e-6), \
        f"Stress not uniform across Gauss points with n_threads={n_threads}"
    assert sv[1, 0] > 0, f"Expected plastic strain with n_threads={n_threads}"
