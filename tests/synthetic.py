"""Synthetic P. berghei AnnData generator for the test suite.

IMPORTANT — scope of this data:
    This is a *software* test fixture, NOT scientific data. It exists so the
    test suite can check that functions run and read/write the correct AnnData
    slots. It must never be used as scientific evidence, and it never appears
    in published results. Scientific validation uses the real PSA dataset and
    public datasets — see publish/decisions/.

Design (Level 2 in our discussion):
    - Real PBANKA_ gene-ID *format* (so ID-parsing code is exercised realistically)
    - Synthetic counts with *deliberate structure*: a subset of "marker" genes
      is up-regulated at specific timepoints, so DEG / enrichment functions
      actually detect signal instead of noise.
    - obs columns mirror the real PSA dataset layout (time, lane, leiden, phase)
      so tests match production usage.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

TIMEPOINTS = ("24h", "48h", "62h")
_PHASES = ("G1", "S", "G2M")


def make_dummy_pberghei(
    n_cells: int = 300,
    n_genes: int = 800,
    n_markers_per_time: int = 20,
    seed: int = 0,
) -> AnnData:
    """Build a small synthetic P. berghei :class:`AnnData` with real ID format.

    Parameters
    ----------
    n_cells
        Number of cells (split roughly evenly across timepoints).
    n_genes
        Number of genes.
    n_markers_per_time
        How many genes are up-regulated in each timepoint (the planted signal
        that DEG/enrichment tests should recover).
    seed
        RNG seed for reproducibility.

    Returns
    -------
    AnnData
        ``X`` holds integer counts; ``layers["counts"]`` is a copy; ``obs`` has
        ``time``/``lane``/``leiden``/``phase``; ``var`` has ``gene_id`` with
        valid ``PBANKA_`` identifiers and an ``is_marker`` flag.
    """
    rng = np.random.default_rng(seed)

    # --- gene IDs in valid PBANKA_ format (PBANKA_dddddd) ---
    gene_ids = [f"PBANKA_{i:06d}" for i in range(1, n_genes + 1)]

    # --- assign cells to timepoints ---
    time = rng.choice(TIMEPOINTS, size=n_cells)
    # two "lanes" (replicates) per timepoint, mirroring PSA structure
    lane = np.array([f"{t}_lane{rng.integers(1, 3)}" for t in time], dtype=object)

    # --- baseline counts: negative binomial-ish via Poisson(gamma) ---
    base_rate = rng.gamma(shape=2.0, scale=1.0, size=n_genes)
    counts = rng.poisson(lam=np.broadcast_to(base_rate, (n_cells, n_genes)))
    counts = counts.astype(np.float32)

    # --- plant timepoint-specific marker signal ---
    is_marker = np.zeros(n_genes, dtype=bool)
    marker_time = np.array([""] * n_genes, dtype=object)
    marker_pool = rng.permutation(n_genes)
    cursor = 0
    for t in TIMEPOINTS:
        idx = marker_pool[cursor : cursor + n_markers_per_time]
        cursor += n_markers_per_time
        is_marker[idx] = True
        marker_time[idx] = t
        # boost those genes in cells of this timepoint
        cell_mask = time == t
        counts[np.ix_(cell_mask, idx)] += rng.poisson(
            lam=15.0, size=(cell_mask.sum(), idx.size)
        ).astype(np.float32)

    obs = pd.DataFrame(
        {
            "time": pd.Categorical(time, categories=list(TIMEPOINTS)),
            "lane": pd.Categorical(lane),
            "leiden": pd.Categorical(rng.integers(0, 4, size=n_cells).astype(str)),
            "phase": pd.Categorical(rng.choice(_PHASES, size=n_cells)),
        },
        index=[f"cell_{i:04d}" for i in range(n_cells)],
    )

    var = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "is_marker": is_marker,
            "marker_time": marker_time,
        },
        index=gene_ids,
    )

    adata = AnnData(X=counts, obs=obs, var=var)
    adata.layers["counts"] = counts.copy()
    adata.uns["organism"] = "pberghei"
    return adata
