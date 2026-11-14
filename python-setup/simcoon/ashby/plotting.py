# -*- coding: utf-8 -*-

"""Ashby diagram plotting for the simcoon Ashby module.

Includes the original legacy visualization helpers (verbatim) and a new
high-level ``ashby_plot()`` function.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import math as m

if TYPE_CHECKING:
    from simcoon.ashby.material import MaterialCollection

# ======================================================================
# Legacy helper functions (preserved verbatim from the original ashby.py)
# ======================================================================

def unit_vector(pt_a, pt_b):
    b_a = [pt_b[0]-pt_a[0],pt_b[1]-pt_a[1]]
    distance = m.sqrt(b_a[0]**2+b_a[1]**2)
    return np.array([b_a[0]*(1.0/distance),b_a[1]*(1.0/distance)])

def length(pt_a, pt_b):
    b_a = [pt_b[0]-pt_a[0],pt_b[1]-pt_a[1]]
    distance = m.sqrt(b_a[0]**2+b_a[1]**2)
    return distance

def uv_2(uv1, uv2):
    b_a = [uv1[0]-uv2[0],uv1[1]-uv2[1]]
    distance = m.sqrt(b_a[0]**2+b_a[1]**2)
    return np.array([b_a[0]*(1.0/distance),b_a[1]*(1.0/distance)])

def poly_convexHull(points, color, coef_multi=0.1, rad=0.3, lw=2):
    """
    Plot the convex hull around a set of points as a
    shaded polygon.
    """
    from scipy.spatial import ConvexHull
    from matplotlib.path import Path
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt

    pt_env = points
    for i in range(0,len(points)):
        pt_env = np.append(pt_env, points[i] + (2.0*coef_multi*np.random.rand(10,2)-coef_multi)*points[i], axis=0)

    hull_env = ConvexHull(pt_env)
    hull_indices_env = hull_env.vertices

    u_v = np.zeros((len(hull_indices_env),2))

    verts = np.zeros((len(hull_indices_env),2))
    dist = np.zeros((len(hull_indices_env),2))

    for i in range(0,len(hull_indices_env)):
        verts[i] = pt_env[hull_indices_env[i]]
        u_v[i] = unit_vector(pt_env[hull_indices_env[i-1]], pt_env[hull_indices_env[i]])
        dist[i] = length(pt_env[hull_indices_env[i-1]], pt_env[hull_indices_env[i]])

    verts2 = np.zeros((3*len(verts)+1,2))
    for i in range(0,len(hull_indices_env)-1):
        verts2[i*3] = verts[i] - u_v[i]*dist[i]*rad
        verts2[i*3+1] = verts[i]
        verts2[i*3+2] = verts[i] + u_v[i+1]*dist[i+1]*rad
    verts2[-4] = verts[-1] - u_v[-1]*dist[-1]*rad
    verts2[-3] = verts[-1]
    verts2[-2] = verts[-1] + u_v[0]*dist[0]*rad
    verts2[-1] = verts2[0]

    codes = [Path.MOVETO,]
    for j in range(len(verts)):
        codes.extend([Path.CURVE3, Path.CURVE3, Path.LINETO,])

    path = Path(verts2, codes)
    patch = patches.PathPatch(path, facecolor=color, lw=0, alpha=0.2)
    edge = patches.PathPatch(path, edgecolor=color, facecolor='none', lw=lw)
    plt.gca().add_patch(patch)
    plt.gca().add_patch(edge)


def poly_enclose(points, color, inc=1.2, rad=0.3, lw=2):
    """
        Plot the convex hull around a set of points as a
        shaded polygon.
        """
    from scipy.spatial import ConvexHull
    from matplotlib.path import Path
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt

    hull = ConvexHull(points)


    cent = np.mean(points, 0)
    pts = []
    for pt in points[hull.vertices]:
        pts.append(pt.tolist())

    pts.insert(len(pts), pts[0])


    verts = inc*(np.array(pts)- cent) + cent
    verts2 = np.zeros((3*verts.shape[0]-2,2))
    verts2[0::3] = verts
    verts2[1::3,:] = (1-rad)*verts[0:-1,:] + rad*verts[1:,:]
    verts2[2::3,:] = rad*verts[0:-1,:] + (1-rad)*verts[1:,:]
    verts2[0:-1] = verts2[1:]
    verts2[-1] = verts2[0]

    codes = [Path.MOVETO, Path.LINETO, Path.CURVE3,]
    for j in range(len(pts)-2):
        codes.extend([Path.CURVE3, Path.LINETO, Path.CURVE3,])
    codes.append(Path.CURVE3)


    path = Path(verts2, codes)
    patch = patches.PathPatch(path, facecolor=color, lw=0, alpha=0.2)
    edge = patches.PathPatch(path, edgecolor=color, facecolor='none', lw=lw)
    plt.gca().add_patch(patch)
    plt.gca().add_patch(edge)

def ellip_enclose(points, color, inc=1.2, lw=2, nst=2):
    """
        Plot the minimum ellipse around a set of points.

        Based on:
        https://github.com/joferkington/oost_paper_code/blob/master/error_ellipse.py
        """
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt

    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    x = points[:,0]
    y = points[:,1]
    cov = np.cov(x, y)
    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    w, h = 2 * nst * np.sqrt(vals)
    center = np.mean(points, 0)
    ell = patches.Ellipse(center, width=inc*w, height=inc*h, angle=theta,
                          facecolor=color, alpha=0.2, lw=0)
    edge = patches.Ellipse(center, width=inc*w, height=inc*h, angle=theta,
                                                 facecolor='none', edgecolor=color, lw=lw)
    plt.gca().add_artist(ell)
    plt.gca().add_artist(edge)

# ======================================================================
# New high-level API
# ======================================================================

CATEGORY_COLORS: dict[str, str] = {
    "Metal": "#1f77b4",
    "Ceramic": "#d62728",
    "Polymer": "#2ca02c",
    "Composite": "#ff7f0e",
    "Foam": "#9467bd",
    "Wood": "#8c564b",
    "Elastomer": "#e377c2",
    "Natural": "#7f7f7f",
}

_PROP_LABELS: dict[str, str] = {
    "E": "Young's modulus (GPa)",
    "G": "Shear modulus (GPa)",
    "K": "Bulk modulus (GPa)",
    "nu": "Poisson's ratio",
    "density": "Density (kg/m³)",
    "CTE": "CTE (1/K)",
    "yield_strength": "Yield strength (MPa)",
    "tensile_strength": "Tensile strength (MPa)",
}


def ashby_plot(
    data,
    x_prop: str = "density",
    y_prop: str = "E",
    *,
    group_by: str = "category",
    envelope: str | None = "ellipse",
    log: bool = True,
    guidelines: list[dict] | None = None,
    ax=None,
    figsize: tuple[float, float] = (10, 7),
    label_groups: bool = True,
):
    """Create an Ashby-style material property chart.

    Parameters
    ----------
    data : MaterialCollection or dict[str, MaterialCollection]
        Materials to plot.  If a single collection, it will be grouped
        using *group_by*.
    x_prop, y_prop : str
        Material attribute names for the x and y axes.
    group_by : str
        Field used to group a single collection (default ``"category"``).
    envelope : str or None
        ``"ellipse"``, ``"convex_hull"``, ``"enclose"``, or ``None``.
    log : bool
        Use log-log axes (default ``True``).
    guidelines : list of dict, optional
        Performance-index guide lines.  Each dict should contain
        ``slope`` (in log-log space) and optionally ``label`` and
        ``intercepts`` (list of floats, log10 of the y-intercept).
    ax : matplotlib Axes, optional
        Axes to draw on; a new figure is created if ``None``.
    figsize : tuple
        Figure size when creating a new figure.
    label_groups : bool
        Place group name labels at the centroid of each group.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    from simcoon.ashby.material import MaterialCollection as _MC

    # Organise into groups
    if isinstance(data, dict):
        groups = data
    else:
        groups = data.group_by(group_by)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    for gname, gcoll in groups.items():
        x, y, names = gcoll.get_property_arrays(x_prop, y_prop)
        if len(x) == 0:
            continue

        color = CATEGORY_COLORS.get(gname, None)
        ax.scatter(x, y, label=gname, c=color, s=30, zorder=3)

        # Envelope (need ≥3 points for hull/ellipse)
        if envelope and len(x) >= 3:
            pts = np.column_stack([x, y])
            if log:
                pts = np.column_stack([np.log10(x), np.log10(y)])
            _color = color or ax._get_lines.get_next_color()
            if envelope == "ellipse":
                _draw_ellipse(ax, pts, _color, log=log)
            elif envelope == "convex_hull":
                _draw_convex_hull(ax, pts, _color, log=log)
            elif envelope == "enclose":
                _draw_enclose(ax, pts, _color, log=log)

        # Group label at centroid
        if label_groups and len(x) > 0:
            cx, cy = np.median(x), np.median(y)
            ax.annotate(
                gname,
                (cx, cy),
                fontsize=9,
                fontweight="bold",
                ha="center",
                va="center",
                alpha=0.7,
                zorder=4,
            )

    if log:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_xlabel(_PROP_LABELS.get(x_prop, x_prop))
    ax.set_ylabel(_PROP_LABELS.get(y_prop, y_prop))
    ax.legend(loc="best", framealpha=0.9)
    ax.grid(True, which="both", ls=":", alpha=0.4)

    # Performance index guidelines
    if guidelines:
        _draw_guidelines(ax, guidelines)

    return ax


# ------------------------------------------------------------------
# Internal drawing helpers
# ------------------------------------------------------------------

def _draw_ellipse(ax, pts_log, color, log=True):
    """Draw a covariance-based ellipse around *pts_log*."""
    import matplotlib.patches as mpatches

    cov = np.cov(pts_log[:, 0], pts_log[:, 1])
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]

    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    w, h = 2 * 2.0 * np.sqrt(np.maximum(vals, 0))
    center = np.mean(pts_log, axis=0)

    if log:
        # Draw in data (log) coordinates via a transform
        from matplotlib.transforms import Affine2D

        tr = Affine2D().rotate_deg(theta).translate(*center)
        tr = tr + ax.transData
        # We cannot use log-scale patches directly; instead draw in
        # log-space then transform back.
        cx, cy = 10 ** center[0], 10 ** center[1]
        # Approximate: just draw on log axes via the axis transform
        ell = mpatches.Ellipse(
            (10 ** center[0], 10 ** center[1]),
            width=10 ** (center[0] + w / 2) - 10 ** (center[0] - w / 2),
            height=10 ** (center[1] + h / 2) - 10 ** (center[1] - h / 2),
            angle=theta,
            facecolor=color,
            alpha=0.15,
            lw=0,
            zorder=1,
        )
        edge = mpatches.Ellipse(
            (10 ** center[0], 10 ** center[1]),
            width=10 ** (center[0] + w / 2) - 10 ** (center[0] - w / 2),
            height=10 ** (center[1] + h / 2) - 10 ** (center[1] - h / 2),
            angle=theta,
            facecolor="none",
            edgecolor=color,
            lw=1.5,
            zorder=1,
        )
    else:
        ell = mpatches.Ellipse(
            center, w, h, angle=theta, facecolor=color, alpha=0.15, lw=0, zorder=1
        )
        edge = mpatches.Ellipse(
            center, w, h, angle=theta, facecolor="none", edgecolor=color, lw=1.5, zorder=1
        )
    ax.add_patch(ell)
    ax.add_patch(edge)


def _draw_convex_hull(ax, pts_log, color, log=True):
    """Draw a convex hull around *pts_log*."""
    from scipy.spatial import ConvexHull
    import matplotlib.patches as mpatches
    from matplotlib.path import Path

    hull = ConvexHull(pts_log)
    verts = pts_log[hull.vertices]
    if log:
        verts = 10 ** verts
    verts = np.vstack([verts, verts[0]])
    path = Path(verts)
    patch = mpatches.PathPatch(path, facecolor=color, alpha=0.15, lw=0, zorder=1)
    edge = mpatches.PathPatch(path, edgecolor=color, facecolor="none", lw=1.5, zorder=1)
    ax.add_patch(patch)
    ax.add_patch(edge)


def _draw_enclose(ax, pts_log, color, log=True):
    """Draw a smooth enclosing hull via the legacy ``poly_enclose``."""
    from scipy.spatial import ConvexHull
    from matplotlib.path import Path
    import matplotlib.patches as mpatches

    hull = ConvexHull(pts_log)
    cent = np.mean(pts_log, 0)
    pts = [pts_log[v].tolist() for v in hull.vertices]
    pts.append(pts[0])

    inc = 1.2
    rad = 0.3
    verts = inc * (np.array(pts) - cent) + cent
    verts2 = np.zeros((3 * verts.shape[0] - 2, 2))
    verts2[0::3] = verts
    verts2[1::3, :] = (1 - rad) * verts[:-1, :] + rad * verts[1:, :]
    verts2[2::3, :] = rad * verts[:-1, :] + (1 - rad) * verts[1:, :]
    verts2[:-1] = verts2[1:]
    verts2[-1] = verts2[0]

    codes = [Path.MOVETO, Path.LINETO, Path.CURVE3]
    for _ in range(len(pts) - 2):
        codes.extend([Path.CURVE3, Path.LINETO, Path.CURVE3])
    codes.append(Path.CURVE3)

    if log:
        verts2 = 10 ** verts2

    path = Path(verts2, codes)
    patch = mpatches.PathPatch(path, facecolor=color, alpha=0.15, lw=0, zorder=1)
    edge = mpatches.PathPatch(path, edgecolor=color, facecolor="none", lw=1.5, zorder=1)
    ax.add_patch(patch)
    ax.add_patch(edge)


def _draw_guidelines(ax, guidelines):
    """Draw performance-index lines on a log-log Ashby plot."""
    xlim = ax.get_xlim()
    log_x = np.linspace(np.log10(xlim[0]), np.log10(xlim[1]), 200)

    for gl in guidelines:
        slope = gl["slope"]
        label = gl.get("label", "")
        intercepts = gl.get("intercepts", [0.0])
        for c in intercepts:
            log_y = slope * log_x + c
            ax.plot(
                10 ** log_x,
                10 ** log_y,
                "--",
                color="gray",
                alpha=0.5,
                lw=1,
                zorder=0,
            )
        if label and intercepts:
            mid = len(log_x) // 2
            c = intercepts[0]
            ax.annotate(
                label,
                (10 ** log_x[mid], 10 ** (slope * log_x[mid] + c)),
                fontsize=8,
                color="gray",
                rotation=np.degrees(np.arctan(slope)),
                ha="center",
                va="bottom",
                zorder=0,
            )
