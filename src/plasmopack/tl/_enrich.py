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

_COLLECTIONS = {"GO", "BP", "CC", "MF"}


def _resolve_gene_sets(
    adata: AnnData, gene_sets: GeneSets | str | Path, organism: str | None
) -> GeneSets:
    """Turn the ``gene_sets`` argument into a GeneSets object.

    Accepts a GeneSets (returned as-is), a collection name ("GO"/"BP"/"CC"/"MF"
    -> bundled defaults for the resolved organism), or a path to a .gmt file.
    """
    if isinstance(gene_sets, GeneSets):
        return gene_sets
    if isinstance(gene_sets, str) and gene_sets in _COLLECTIONS:
        from plasmopack import datasets as _ds

        org = organism or adata.uns.get("organism") or _ds.detect_organism(adata)
        if org is None:
            raise ValueError(
                f"gene_sets={gene_sets!r} needs an organism, but it could not be "
                f"detected from the gene IDs. Pass organism='pberghei' (or set "
                f"adata.uns['organism'])."
            )
        return _ds.gene_sets(str(org), collection=gene_sets)
    return read_gmt(gene_sets)


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
    top_n: int | None,
    rank_by: str,
) -> list[str]:
    store = adata.uns.get("deg", {})
    if groupby not in store or deg_method not in store[groupby]:
        raise KeyError(
            f"No DEG result at adata.uns['deg'][{groupby!r}][{deg_method!r}]. "
            f"Run pp.tl.deg(adata, {groupby!r}, method={deg_method!r}) first."
        )
    df = store[groupby][deg_method].copy()
    df = df[df["group"].astype(str) == str(group)]

    # 1. direction filter (always applied)
    if direction == "up":
        df = df[df["logFC"] > 0]
    elif direction == "down":
        df = df[df["logFC"] < 0]
    elif direction != "both":
        raise ValueError("direction must be 'up', 'down', or 'both'.")

    if top_n is None:
        # 2a. CUTOFF selection: combined padj AND |logFC| thresholds (default)
        sig = (df["padj"] < padj_cutoff) & (df["logFC"].abs() >= lfc_cutoff)
        return [str(g) for g in df.loc[sig, "gene"]]

    # 2b. TOP-N selection: rank then take the best N
    nlp = -np.log10(np.clip(df["padj"].to_numpy(), 1e-300, 1.0))
    alfc = df["logFC"].abs().to_numpy()
    if rank_by == "padj":
        score = nlp
    elif rank_by == "logfc":
        score = alfc
    elif rank_by == "combined":
        score = nlp * alfc  # significance weighted by effect size
    else:
        raise ValueError("rank_by must be 'padj', 'logfc', or 'combined'.")
    order = np.argsort(score)[::-1][:top_n]
    return [str(g) for g in df["gene"].to_numpy()[order]]


def enrich(
    adata: AnnData,
    gene_sets: GeneSets | str | Path,
    *,
    organism: str | None = None,
    genes: Sequence[str] | None = None,
    groupby: str | None = None,
    group: str | None = None,
    deg_method: str = "wilcoxon",
    direction: str = "up",
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
    top_n: int | None = None,
    rank_by: str = "combined",
    background: Iterable[str] | None = None,
    min_set: int = 5,
    max_set: int = 500,
    key_added: str | None = None,
) -> AnnData:
    """Over-representation analysis of a query gene list against gene sets.

    Provide the query genes one of two ways:

    1. explicitly via ``genes=[...]``, or
    2. from a stored DEG result via ``groupby`` + ``group``.

    When pulling from DEG, choose **how** genes are selected:

    - **Cutoff (default):** significant *and* strong — ``padj < padj_cutoff``
      **and** ``|logFC| >= lfc_cutoff`` (after the ``direction`` filter).
    - **Top-N ranked:** set ``top_n`` to take the best N genes ranked by
      ``rank_by``: ``"padj"`` (most significant), ``"logfc"`` (largest fold
      change), or ``"combined"`` (``-log10(padj) * |logFC|`` — significance
      weighted by effect size; the recommended single ranking).

    Parameters
    ----------
    adata
        Annotated data; result written to ``adata.uns["enrichment"]``.
    gene_sets
        One of: a :class:`~plasmopack.io.GeneSets`; a path to a ``.gmt`` file;
        or a **collection name** string — ``"GO"`` (all aspects), ``"BP"``,
        ``"CC"`` or ``"MF"`` — which loads the bundled default gene sets for the
        organism. The organism is taken from ``organism=``, else
        ``adata.uns["organism"]``, else auto-detected from the gene-ID pattern
        in ``adata.var_names``.
    organism
        Organism key/alias for resolving a collection-name ``gene_sets``. If not
        given, it is detected automatically.
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
    gene_sets = _resolve_gene_sets(adata, gene_sets, organism)

    if genes is not None:
        query = [str(g) for g in genes]
        source_label = key_added or "custom"
    elif groupby is not None and group is not None:
        query = _query_from_deg(
            adata,
            groupby,
            group,
            deg_method,
            direction,
            padj_cutoff,
            lfc_cutoff,
            top_n,
            rank_by,
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
