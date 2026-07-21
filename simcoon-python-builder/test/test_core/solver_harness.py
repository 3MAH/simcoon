"""Shared file-solver harness for solver-level tests.

One copy of the path-file grammar, the output layout, and the temp-dir runner
used by test_modular_finite.py and test_solver_robustness.py (not collected by
pytest: no test_ prefix).
"""

import os

import numpy as np

import simcoon as sim

# Column layout of res_global-0.txt with the OUTPUT_DAT block below:
# 3 inc | 4 time | 8:14 log strain | 14:20 Kirchhoff stress | 20:29 R |
# 29:38 F | 38:42 Wm, Wm_r, Wm_ir, Wm_d
C_TIME = 4
S_STRAIN = slice(8, 14)
S_STRESS = slice(14, 20)
S_WM = slice(38, 42)

OUTPUT_DAT = """#Output_values
strain_type 3
nb_strain   6
0   1   2   3   4   5
stress_type\t3
nb_stress   6
0   1   2   3   4   5

Rotation_type\t1
Tangent_type\t0
T   1
T   1

Number_of_wanted_internal_variables\t0

#Block #type_1_N_2_T    #every
1      1                1
"""


def path_file(targets, control_type):
    """Uniaxial path, one mode per target: ("S", v) drives the 11 stress to v,
    ("E", v) drives the 11 strain to v; lateral components stay stress-free."""
    txt = f"""#Initial_temperature
290
#Number_of_blocks
1

#Block
1
#Loading_type
1
#Control_type(NLGEOM)
{control_type}
#Repeat
1
#Steps
{len(targets)}
"""
    spin = """#spin
0. 0. 0.
0. 0. 0.
0. 0. 0.
""" if 2 <= control_type <= 4 else ""
    for ctrl, target in targets:
        txt += f"""
#Mode
1
#Dn_init 1.
#Dn_mini 0.001
#Dn_inc 0.02
#time
1.
#mechanical_state
{ctrl} {target}
S 0 S 0
S 0 S 0 S 0
{spin}#temperature_state
T 290
"""
    return txt


def run_path(base_dir, umat_name, props, nstatev, corate, path_text):
    """Run the file solver in base_dir (a pytest tmp_path); return the global
    history array of res_global-0.txt."""
    dd, rd = os.path.join(base_dir, "data"), os.path.join(base_dir, "results")
    os.makedirs(dd, exist_ok=True)
    os.makedirs(rd, exist_ok=True)
    with open(os.path.join(dd, "path.txt"), "w") as f:
        f.write(path_text)
    with open(os.path.join(dd, "output.dat"), "w") as f:
        f.write(OUTPUT_DAT)
    sim._core.solver(
        umat_name, np.asarray(props, dtype=float), nstatev,
        0.0, 0.0, 0.0, 0, corate, dd, rd, "path.txt", "res.txt",
    )
    return np.loadtxt(os.path.join(rd, "res_global-0.txt"), ndmin=2)
