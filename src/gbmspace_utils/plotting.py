"""Consistent spatial plotting helpers for the GBM-Space Level 2 (Visium) notebook."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from anndata import AnnData
from matplotlib.figure import Figure


def plot_gene_on_tissue(
    adata: AnnData,
    gene: str,
    spatial_key: str = "spatial",
    layer: str | None = None,
    cmap: str = "Reds",
    point_size: float = 8.0,
    figsize: tuple[float, float] = (6, 6),
    ax: "plt.Axes | None" = None,
) -> Figure:
    """Single-gene spatial expression scatter. Caps `vmax` at the 95th percentile of
    nonzero expression so a handful of outlier spots don't wash out the color scale.
    Plots into `ax` if given (e.g. one panel of an outer `plt.subplots` grid), otherwise
    creates its own standalone figure.
    """
    if gene not in adata.var_names:
        raise ValueError(f"Gene '{gene}' not found in adata.var_names")

    expr = adata[:, gene].layers[layer] if layer else adata[:, gene].X
    expr = np.asarray(expr.todense()).flatten() if hasattr(expr, "todense") else np.asarray(expr).flatten()
    coords = adata.obsm[spatial_key]

    nonzero = expr[expr > 0]
    vmax = np.percentile(nonzero, 95) if len(nonzero) else None

    fig = ax.figure if ax is not None else plt.subplots(figsize=figsize)[0]
    ax = ax if ax is not None else fig.axes[0]
    sc_plot = ax.scatter(coords[:, 0], coords[:, 1], c=expr, cmap=cmap, s=point_size, vmax=vmax)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_title(gene)
    ax.axis("off")
    fig.colorbar(sc_plot, ax=ax, shrink=0.7, label="expression")
    return fig


def plot_spatial_categories(
    adata: AnnData,
    category_key: str,
    spatial_key: str = "spatial",
    categories: list[str] | None = None,
    point_size: float = 8.0,
    figsize: tuple[float, float] = (7, 6),
    ax: "plt.Axes | None" = None,
) -> Figure:
    """Spatial scatter colored by a categorical `.obs` column (cell type, niche, cluster, ...),
    with a consistent tab20-based color map and an external legend. Plots into `ax` if given
    (e.g. one panel of an outer `plt.subplots` grid), otherwise creates its own standalone figure.
    """
    cats = categories or list(adata.obs[category_key].astype("category").cat.categories)
    colors = dict(zip(cats, plt.cm.tab20.colors[: len(cats)]))
    coords = adata.obsm[spatial_key]

    fig = ax.figure if ax is not None else plt.subplots(figsize=figsize)[0]
    ax = ax if ax is not None else fig.axes[0]
    for cat in cats:
        mask = (adata.obs[category_key] == cat).to_numpy()
        ax.scatter(coords[mask, 0], coords[mask, 1], c=[colors[cat]], s=point_size, label=cat)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.axis("off")
    fig.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8, markerscale=2)
    return fig
