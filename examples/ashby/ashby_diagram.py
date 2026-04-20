#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Example: Ashby material property chart with simcoon bridge.

Demonstrates the full workflow:
1. Load the built-in curated material dataset
2. Create a Young's modulus vs density Ashby diagram with envelopes
3. Add performance-index guide lines (E/rho = const for light-stiff design)
4. Pick a candidate material and convert to a simcoon stiffness tensor
5. (Optional) Fetch additional materials from Materials Project
"""

import numpy as np

# =====================================================================
# 1. Load built-in dataset
# =====================================================================

from simcoon.ashby import load_builtin, Material, MaterialCollection

mats = load_builtin()
print(f"Loaded {len(mats)} materials")
print(f"Categories: {sorted(set(m.category for m in mats))}")

# Quick peek at the first few metals
metals = mats.filter(category="Metal")
print(f"\nMetals ({len(metals)}):")
for m in metals[:5]:
    print(f"  {m.name:30s}  E={m.E:7.1f} GPa  rho={m.density:.0f} kg/m3")

# =====================================================================
# 2. Ashby diagram: E vs density with ellipse envelopes
# =====================================================================

from simcoon.ashby import ashby_plot
import matplotlib.pyplot as plt

ax = ashby_plot(
    mats,
    x_prop="density",
    y_prop="E",
    envelope="ellipse",
    guidelines=[
        {
            "slope": 1.0,
            "label": "E/ρ = const  (tie rod)",
            "intercepts": [-1.0, 0.0, 1.0],
        },
        {
            "slope": 0.5,
            "label": "E^0.5/ρ = const  (beam)",
            "intercepts": [-0.5, 0.5, 1.5],
        },
    ],
)
ax.set_title("Ashby Diagram — Young's Modulus vs Density")
plt.tight_layout()
plt.savefig("ashby_E_vs_density.png", dpi=150)
print("\nSaved ashby_E_vs_density.png")
plt.show()

# =====================================================================
# 3. Pick a candidate and convert to simcoon stiffness tensor
# =====================================================================

from simcoon.ashby import to_stiffness, to_solver_props

# Select Al 6061-T6
al = [m for m in mats if "6061" in m.name][0]
print(f"\nSelected material: {al.name}")
print(f"  E  = {al.E} GPa")
print(f"  nu = {al.nu}")
print(f"  G  = {al.G} GPa")

# Stiffness tensor (6x6 in MPa)
L = to_stiffness(al)
print(f"\nStiffness tensor L (MPa):\n{np.array2string(L, precision=1, suppress_small=True)}")

# Solver properties
sp = to_solver_props(al)
print(f"\nSolver props: umat_name={sp['umat_name']!r}, props={sp['props']}")

# =====================================================================
# 4. (Optional) Fetch from Materials Project
# =====================================================================
# Uncomment the block below if you have mp-api installed and an API key.
#
# from simcoon.ashby import fetch_materials
#
# mp_mats = fetch_materials(elements=["Al", "O"], limit=20)
# print(f"\nFetched {len(mp_mats)} materials from Materials Project")
# for m in mp_mats[:5]:
#     print(f"  {m.source_id}  {m.name:20s}  E={m.E:.1f} GPa")
#
# # Merge with built-in dataset and re-plot
# combined = MaterialCollection(list(mats) + list(mp_mats))
# ashby_plot(combined, "density", "E")
# plt.title("Enriched Ashby Diagram")
# plt.show()
