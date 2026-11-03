"""Pathway enrichment — over-representation analysis (ORA).

This is "Mode B" from the enrichment strategy: a local, reproducible
over-representation test. Given a query gene list (typically the significant
DEGs for one group) and a set of gene sets (a :class:`~plasmopack.io.GeneSets`,
e.g. built from a PlasmoDB GAF), it asks for each gene set: *are the query
genes over-represented in this set, beyond chance?*

The test is Fisher's exact test (one-sided, "greater") on the 2x2 table of
(in query / not) x (in set / not), restricted to a background gene universe.
P-values are Benjamini-Hochberg adjusted across the tested sets. This is the
same statistic PlasmoDB's web tool uses, so results are directly comparable.

Results go to ``adata.uns["enrichment"][label]`` as a tidy DataFrame.

Modes A (PlasmoDB API) and C (GSEA) are planned; see
``publish/decisions/pathway_enrichment_strategy.md``.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats

from plasmopack._utils.provenance import record_provenance
from plasmopack.io import GeneSets, read_gmt


def _bh(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0.0, 1.0)
    return out


def _query_from_deg(
    adata: AnnData,
    groupby: str,
    group: str,
    deg_method: str,
    direction: str,
    padj_cutoff: float,
    lfc_cutoff: float,
) -> list[str]:
    store = adata.uns.get("deg", {})
    if groupby not in store or deg_method not in store[groupby]:
        raise KeyError(
            f"No DEG result at adata.uns['deg'][{groupby!r}][{deg_method!r}]. "
            f"Run pp.tl.deg(adata, {groupby!r}, method={deg_method!r}) first."
        )
    df = store[groupby][deg_method]
    df = df[df["group"].astype(str) == str(group)]
    sig = df["padj"] < padj_cutoff
    if direction == "up":
        sig &= df["logFC"] > lfc_cutoff
    elif direction == "down":
        sig &= df["logFC"] < -lfc_cutoff
    elif direction == "both":
        sig &= df["logFC"].abs() > lfc_cutoff
    else:
        raise ValueError("direction must be 'up', 'down', or 'both'.")
    return [str(g) for g in df.loc[sig, "gene"]]


def enrich(
    adata: AnnData,
    gene_sets: GeneSets | str | Path,
    *,
    genes: Sequence[str] | None = None,
    groupby: str | None = None,
    group: str | None = None,
    deg_method: str = "wilcoxon",
    direction: str = "up",
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
    background: Iterable[str] | None = None,
    min_set: int = 5,
    max_set: int = 500,
    key_added: str | None = None,
) -> AnnData:
    """Over-representation analysis of a query gene list against gene sets.

    Provide the query genes one of two ways:

    1. explicitly via ``genes=[...]``, or
    2. from a stored DEG result via ``groupby`` + ``group`` (filtered by
       ``direction``/``padj_cutoff``/``lfc_cutoff``).

    Parameters
    ----------
    adata
        Annotated data; result written to ``adata.uns["enrichment"]``.
    gene_sets
        A :class:`~plasmopack.io.GeneSets`, or a path to a ``.gmt`` file.
    genes
        Explicit query gene list. Overrides the DEG path if given.
    groupby, group
        Pull the query from ``adata.uns["deg"][groupby][deg_method]`` for this
        group when ``genes`` is not provided.
    deg_method, direction, padj_cutoff, lfc_cutoff
        How to select significant genes from the DEG table. ``direction`` is
        ``"up"``, ``"down"``, or ``"both"``.
    background
        The gene universe (default: ``adata.var_names``). The test and set-size
        filtering are restricted to this universe.
    min_set, max_set
        Keep only gene sets whose size (within the background) is in this range
        (standard GSEA practice).
    key_added
        Label under ``adata.uns["enrichment"]``. Defaults to a descriptive key.

    Returns
    -------
    AnnData
        The same object; results in ``adata.uns["enrichment"][key]`` as a
        DataFrame (set_id, name, set_size, n_overlap, odds_ratio, pval, padj,
        overlap_genes), sorted by padj.
    """
    if isinstance(gene_sets, (str, Path)):
        gene_sets = read_gmt(gene_sets)

    if genes is not None:
        query = [str(g) for g in genes]
        source_label = key_added or "custom"
    elif groupby is not None and group is not None:
        query = _query_from_deg(
            adata, groupby, group, deg_method, direction, padj_cutoff, lfc_cutoff
        )
        source_label = key_added or f"{groupby}:{group}:{direction}"
    else:
        raise ValueError("Provide either genes=[...] or both groupby= and group=.")

    universe = (
        {str(g) for g in background}
        if background is not None
        else set(map(str, adata.var_names))
    )
    query_set = set(query) & universe
    n_query = len(query_set)
    big_n = len(universe)

    rows = []
    for gs in gene_sets:
        members = set(gs.genes) & universe
        k_set = len(members)
        if k_set < min_set or k_set > max_set:
            continue
        overlap = query_set & members
        k = len(overlap)
        # 2x2: [[in query & in set, in query & not set],
        #       [not query & in set, not query & not set]]
        table = [
            [k, n_query - k],
            [k_set - k, big_n - n_query - (k_set - k)],
        ]
        odds, pval = stats.fisher_exact(table, alternative="greater")
        rows.append(
            {
                "set_id": gs.id,
                "name": gs.name,
                "description": gs.description,
                "set_size": k_set,
                "n_overlap": k,
                "odds_ratio": float(odds),
                "pval": float(pval),
                "overlap_genes": ",".join(sorted(overlap)),
            }
        )

    result = pd.DataFrame(rows)
    if not result.empty:
        result["padj"] = _bh(result["pval"].to_numpy())
        result = result.sort_values("padj", kind="stable").reset_index(drop=True)
    else:
        result = pd.DataFrame(
            columns=[
                "set_id",
                "name",
                "description",
                "set_size",
                "n_overlap",
                "odds_ratio",
                "pval",
                "padj",
                "overlap_genes",
            ]
        )

    store = adata.uns.setdefault("enrichment", {})
    if not isinstance(store, dict):
        store = {}
        adata.uns["enrichment"] = store
    store[source_label] = result

    record_provenance(
        adata,
        "tl.enrich",
        params={
            "method": "ora",
            "label": source_label,
            "n_query": n_query,
            "n_background": big_n,
            "min_set": min_set,
            "max_set": max_set,
            "direction": direction if genes is None else None,
        },
        sources={"gmt_version": str(gene_sets.metadata.get("version", "unknown"))},
    )
    return adata
