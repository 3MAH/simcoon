"""Gallery thumbnails rendered from the equations of the examples that draw no figure.

Each entry of ``EQUATIONS`` is written to ``_static/thumbnails/eq_<stem>.png`` at
the start of the build; the example points at it with
``# sphinx_gallery_thumbnail_path = '_static/thumbnails/eq_<stem>.png'``.
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

EQUATIONS = {
    "constitutive_relations": [
        r"$\mathbf{L} = 3\kappa\,\mathbb{J} + 2\mu\,\mathbb{K}$",
        r"$\boldsymbol{\sigma} = \mathbf{L} : \boldsymbol{\varepsilon}, \quad \mathbf{M} = \mathbf{L}^{-1}$",
    ],
    "rotation": [
        r"$\boldsymbol{\sigma}' = \mathbf{R}\,\boldsymbol{\sigma}\,\mathbf{R}^T$",
        r"$\mathbf{L}' = \mathbf{Q}_\sigma\,\mathbf{L}\,\mathbf{Q}_\sigma^T$",
    ],
    "tensor": [
        r"$\boldsymbol{\tau} = \mathbf{F}\,\mathbf{S}\,\mathbf{F}^T$",
        r"$\mathbf{e} = \mathbf{F}^{-T}\,\mathbf{E}\,\mathbf{F}^{-1}$",
    ],
    "tensor_basis": [
        r"$\mathbf{g}_i = \mathbf{F}\,\mathbf{G}_i, \quad g_{ij} = \mathbf{g}_i \cdot \mathbf{g}_j$",
        r"$\mathrm{tr}\,\boldsymbol{\tau} = g_{ij}\,\tau^{ij}$",
    ],
    "stress_transfer_helpers": [
        r"$\boldsymbol{\tau} = J\boldsymbol{\sigma}, \quad \mathbf{P} = J\boldsymbol{\sigma}\,\mathbf{F}^{-T}$",
        r"$\mathbf{S} = J\,\mathbf{F}^{-1}\boldsymbol{\sigma}\,\mathbf{F}^{-T}$",
    ],
    "eshelby_tensors": [
        r"$\boldsymbol{\varepsilon}^c = \mathbb{S} : \boldsymbol{\varepsilon}^*$",
        r"$\mathbb{S}^{\mathrm{sph}} = \frac{1+\nu}{3(1-\nu)}\,\mathbb{J} + \frac{2(4-5\nu)}{15(1-\nu)}\,\mathbb{K}$",
    ],
}

_NAVY = "#24445c"
_FONTSIZE = 26
_MAX_WIDTH = 0.94


def build(outdir):
    """Write one PNG per entry of ``EQUATIONS`` into ``outdir``."""
    os.makedirs(outdir, exist_ok=True)
    rc = {
        "text.usetex": True,
        "text.latex.preamble": r"\usepackage{amsmath,amssymb}",
        "font.family": "serif",
    }
    with plt.rc_context(rc):
        for stem, lines in EQUATIONS.items():
            fig = plt.figure(figsize=(4.0, 2.8), dpi=100)
            n = len(lines)
            renderer = fig.canvas.get_renderer()
            for k, line in enumerate(lines):
                text = fig.text(0.5, (n - k) / (n + 1), line, ha="center", va="center",
                                fontsize=_FONTSIZE, color=_NAVY)
                width = text.get_window_extent(renderer).width
                if width > _MAX_WIDTH * fig.bbox.width:
                    text.set_fontsize(_FONTSIZE * _MAX_WIDTH * fig.bbox.width / width)
            fig.savefig(os.path.join(outdir, f"eq_{stem}.png"), facecolor="white")
            plt.close(fig)
