"""Quality-control metrics and filtering for AnnData.

These operate in place on an :class:`~anndata.AnnData`, writing per-cell QC
columns into ``adata.obs`` and stamping provenance. The gene-group percentages
(e.g. mitochondrial, apicoplast) are the Plasmodium-relevant QC signals.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

import numpy as np
import scipy.sparse as sp
from anndata import AnnData

from plasmopack._utils.provenance import record_provenance


def _counts_per_cell(x: Any) -> np.ndarray:
    """Total counts per cell (row sums), as a dense 1-D array."""
    sums = x.sum(axis=1)
    return np.asarray(sums).ravel()


def _genes_per_cell(x: Any) -> np.ndarray:
    """Number of genes with non-zero counts per cell."""
    if sp.issparse(x):
        nz = x.getnnz(axis=1)
        return np.asarray(nz).ravel()
    return np.asarray((x > 0).sum(axis=1)).ravel()


def qc_metrics(
    adata: AnnData,
    *,
    gene_groups: Mapping[str, Iterable[str]] | None = None,
    layer: str | None = None,
) -> AnnData:
    """Compute per-cell QC metrics in place.

    Always writes:
    - ``adata.obs["total_counts"]`` — total counts per cell
    - ``adata.obs["n_genes_by_counts"]`` — number of detected genes per cell

    For each entry in ``gene_groups`` also writes
    ``adata.obs["pct_counts_<name>"]`` — the percentage of each cell's counts
    coming from that gene group (0 for cells with no counts).

    Parameters
    ----------
    adata
        Annotated data; modified in place and also returned for chaining.
    gene_groups
        Mapping of group name -> gene identifiers (matched against
        ``adata.var_names``). Typical for Plasmodium::

            {"mito": ["PBANKA_MIT..."], "apicoplast": [...]}
    layer
        Use this layer instead of ``adata.X`` (e.g. raw counts).

    Returns
    -------
    AnnData
        The same object, for chaining.
    """
    x = adata.layers[layer] if layer is not None else adata.X

    total = _counts_per_cell(x)
    adata.obs["total_counts"] = total
    adata.obs["n_genes_by_counts"] = _genes_per_cell(x)

    groups_recorded: dict[str, int] = {}
    if gene_groups:
        var_index = adata.var_names
        with np.errstate(invalid="ignore", divide="ignore"):
            safe_total = np.where(total > 0, total, 1.0)
            for name, genes in gene_groups.items():
                mask = var_index.isin(list(genes))
                n_matched = int(mask.sum())
                groups_recorded[name] = n_matched
                sub = x[:, mask]
                group_counts = (
                    _counts_per_cell(sub) if n_matched else np.zeros_like(total)
                )
                pct = np.where(total > 0, 100.0 * group_counts / safe_total, 0.0)
                adata.obs[f"pct_counts_{name}"] = pct

    record_provenance(
        adata,
        "pp.qc_metrics",
        params={"gene_groups": groups_recorded, "layer": layer},
    )
    return adata


def filter_cells(
    adata: AnnData,
    *,
    min_genes: int | None = None,
    min_counts: float | None = None,
    max_counts: float | None = None,
) -> AnnData:
    """Filter cells in place by detected genes and/or total counts.

    Thresholds are inclusive lower bounds (``min_*``) and an inclusive upper
    bound (``max_counts``). ``None`` disables that threshold.

    Returns
    -------
    AnnData
        The filtered object (same instance), for chaining.
    """
    x = adata.X
    total = _counts_per_cell(x)
    n_genes = _genes_per_cell(x)

    keep = np.ones(adata.n_obs, dtype=bool)
    if min_genes is not None:
        keep &= n_genes >= min_genes
    if min_counts is not None:
        keep &= total >= min_counts
    if max_counts is not None:
        keep &= total <= max_counts

    n_before = adata.n_obs
    adata._inplace_subset_obs(keep)

    record_provenance(
        adata,
        "pp.filter_cells",
        params={
            "min_genes": min_genes,
            "min_counts": min_counts,
            "max_counts": max_counts,
            "n_removed": int(n_before - adata.n_obs),
        },
    )
    return adata


def filter_genes(
    adata: AnnData,
    *,
    min_cells: int | None = None,
    min_counts: float | None = None,
) -> AnnData:
    """Filter genes in place by the number of cells they are detected in
    and/or their total counts across cells.

    Returns
    -------
    AnnData
        The filtered object (same instance), for chaining.
    """
    x = adata.X
    if sp.issparse(x):
        cells_per_gene = np.asarray(x.getnnz(axis=0)).ravel()
    else:
        cells_per_gene = np.asarray((x > 0).sum(axis=0)).ravel()
    counts_per_gene = np.asarray(x.sum(axis=0)).ravel()

    keep = np.ones(adata.n_vars, dtype=bool)
    if min_cells is not None:
        keep &= cells_per_gene >= min_cells
    if min_counts is not None:
        keep &= counts_per_gene >= min_counts

    n_before = adata.n_vars
    adata._inplace_subset_var(keep)

    record_provenance(
        adata,
        "pp.filter_genes",
        params={
            "min_cells": min_cells,
            "min_counts": min_counts,
            "n_removed": int(n_before - adata.n_vars),
        },
    )
    return adata
