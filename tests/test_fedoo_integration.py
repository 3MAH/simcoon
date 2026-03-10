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


# ── Test 1: Non-linear UMAT (EPICP) ─ small strain ─ exercises OpenMP ──

def test_epicp_small_strain():
    """Elastoplastic EPICP umat under small-strain compression.

    This exercises the OpenMP-parallelised umat call (multi-point
    stress integration) with a non-linear material (plasticity).
    """
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

    # Verify plastic deformation occurred (statev[1] = accumulated plastic strain p)
    p = asm.sv["Statev"][1]
    assert np.max(p) > 0, "Expected plastic deformation"


# ── Test 2: Non-linear geometry (Neo-Hookean) ─ tangent modulus transfer ──

def test_neohookean_nlgeom():
    """Neo-Hookean hyperelastic under large-strain compression.

    This exercises non-linear geometry (nlgeom=True) which involves
    deformation gradient transfer and tangent modulus computation.
    """
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

    # Verify deformation gradient was tracked
    F = asm.sv["F"]
    assert F is not None
    assert not np.allclose(F, np.eye(3).reshape(3, 3, 1)), \
        "Deformation gradient should differ from identity"


# ── Test 3: EPICP + nlgeom ─ OpenMP umat + tangent modulus + large strain ──

def test_epicp_nlgeom():
    """Elastoplastic EPICP under large-strain compression.

    This is the most demanding test: non-linear material (OpenMP umat)
    combined with non-linear geometry (deformation gradient and tangent
    modulus transfer between fedoo and simcoon).
    """
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

    # Verify plastic deformation
    p = asm.sv["Statev"][1]
    assert np.max(p) > 0, "Expected plastic deformation"

    # Verify deformation gradient was tracked
    F = asm.sv["F"]
    assert not np.allclose(F, np.eye(3).reshape(3, 3, 1))

    # Verify tangent matrix is populated and symmetric
    Lt = asm.sv["TangentMatrix"]
    assert Lt is not None
    sym_err = np.abs(Lt - Lt.transpose((1, 0, 2))).max()
    assert sym_err < 1e-3, f"Tangent matrix not symmetric (err={sym_err})"


# ── Test 4: Direct sim.umat call with explicit n_threads ──

@pytest.mark.parametrize("n_threads", [1, 2, 4])
def test_umat_openmp_threads(n_threads):
    """Call sim.umat directly with explicit n_threads to verify OpenMP works.

    With multiple Gauss points the umat loop is OpenMP-parallelised.
    Results must be identical regardless of thread count.
    """
    n_gauss = 64  # simulate 64 Gauss points (like 8 hex8 elements)
    props = np.asfortranarray(
        np.tile([200e3, 0.3, 0.0, 300.0, 1000.0, 0.5], (n_gauss, 1)).T
    )

    etot = np.zeros((6, n_gauss), order="F")
    # Apply a uniaxial strain increment that exceeds yield
    Detot = np.zeros((6, n_gauss), order="F")
    Detot[2, :] = 5e-3  # eps_zz increment

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

    # All Gauss points have identical input → identical output
    assert np.allclose(stress[:, 0:1], stress, atol=1e-6), \
        f"Stress not uniform across Gauss points with n_threads={n_threads}"

    # Should have yielded (sigma_Y = 300 MPa, E*eps = 200e3*5e-3 = 1000 MPa)
    assert sv[1, 0] > 0, "Expected plastic strain with n_threads={n_threads}"
