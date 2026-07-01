"""Phase 1 tests — the synthetic data generator itself.

The fixture must be reproducible, have valid gene-ID format, and actually
contain the planted marker signal (otherwise downstream DEG/enrichment tests
would be testing nothing).
"""

from __future__ import annotations

import numpy as np

from plasmopack.datasets import get_organism
from tests.synthetic import TIMEPOINTS, make_dummy_pberghei


def test_shape_and_slots() -> None:
    adata = make_dummy_pberghei(n_cells=120, n_genes=200, seed=1)
    assert adata.shape == (120, 200)
    assert "counts" in adata.layers
    assert set(adata.obs["time"].cat.categories) == set(TIMEPOINTS)
    assert {"gene_id", "is_marker", "marker_time"} <= set(adata.var.columns)


def test_gene_ids_valid_pberghei_format() -> None:
    adata = make_dummy_pberghei(n_cells=50, n_genes=100, seed=2)
    pb = get_organism("pberghei")
    assert all(pb.matches_gene_id(g) for g in adata.var["gene_id"])


def test_reproducible_with_seed() -> None:
    a = make_dummy_pberghei(n_cells=80, n_genes=150, seed=42)
    b = make_dummy_pberghei(n_cells=80, n_genes=150, seed=42)
    assert np.array_equal(a.X, b.X)


def test_planted_signal_is_detectable() -> None:
    """Marker genes should be higher in their own timepoint than elsewhere."""
    adata = make_dummy_pberghei(n_cells=600, n_genes=400, seed=3)
    for t in TIMEPOINTS:
        markers = adata.var["marker_time"] == t
        in_time = adata.obs["time"] == t
        # mean expression of this timepoint's markers, inside vs outside
        inside = adata[in_time, markers].X.mean()
        outside = adata[~in_time, markers].X.mean()
        assert inside > outside, f"planted signal missing for {t}"
