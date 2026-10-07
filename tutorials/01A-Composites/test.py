import numpy as np
import simcoon as sim
from simcoon.solver.micromechanics import Ellipsoid

nstatev = 0

int1 = 50  # Integration points of the Eshelby integrals, first direction
int2 = 50  # Integration points of the Eshelby integrals, second direction
n_matrix = 0  # Index of the matrix phase in the list of phases

# Mori-Tanaka props: [int1, int2, n_matrix]; the phases themselves are passed as objects
props = np.array([int1, int2, n_matrix], dtype="float")

matrix = Ellipsoid(
    number=0, umat_name="ELISO", save=1, concentration=0.8, nstatev=1,
    props=np.array([2250.0, 0.19, 8.8e-5]),
)
reinforcement = Ellipsoid(
    number=1, umat_name="ELISO", save=1, concentration=0.2, nstatev=1,
    props=np.array([73000.0, 0.19, 0.5e-6]),
)

psi_rve = 0.0
theta_rve = 0.0
phi_rve = 0.0

umat_name = "MIMTN"

L = sim.L_eff(
    umat_name, props, nstatev, orientation=(psi_rve, theta_rve, phi_rve),
    phases=[matrix, reinforcement],
)
p = sim.L_iso_props(L)
print(p)
