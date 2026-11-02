"""Phase 4 tests — normalization."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from anndata import AnnData

import plasmopack as pp
from plasmopack._utils.provenance import get_history


def _toy() -> AnnData:
    x = np.array([[1, 1, 2, 0], [0, 0, 5, 5]], dtype=np.float32)
    return AnnData(X=x)


def test_normalize_writes_layer_and_keeps_X() -> None:
    adata = _toy()
    x_before = adata.X.copy()
    pp.pp.normalize(adata, target_sum=1e4)
    assert "lognorm" in adata.layers
    # X untouched (non-destructive)
    np.testing.assert_array_equal(adata.X, x_before)


def test_normalize_scales_to_target_before_log() -> None:
    adata = _toy()
    pp.pp.normalize(adata, target_sum=100.0, log1p=False)
    row_sums = adata.layers["lognorm"].sum(axis=1)
    np.testing.assert_allclose(np.asarray(row_sums).ravel(), [100.0, 100.0])


def test_normalize_log1p_applied() -> None:
    adata = _toy()
    pp.pp.normalize(adata, target_sum=100.0, log1p=True)
    # with log1p, values are log1p of the scaled counts; max should be < scaled max
    scaled_max = 100.0  # a single cell fully in one gene would be 100
    assert adata.layers["lognorm"].max() < scaled_max


def test_normalize_handles_zero_count_cell() -> None:
    adata = AnnData(X=np.array([[0, 0, 0], [1, 2, 3]], dtype=np.float32))
    pp.pp.normalize(adata, target_sum=10.0, log1p=False)
    # zero-count cell stays all zeros, no NaN/inf
    assert np.isfinite(adata.layers["lognorm"]).all()
    np.testing.assert_array_equal(adata.layers["lognorm"][0], [0, 0, 0])


def test_normalize_sparse_matches_dense() -> None:
    dense = _toy()
    sparse = _toy()
    sparse.X = sp.csr_matrix(sparse.X)
    pp.pp.normalize(dense, target_sum=1e4)
    pp.pp.normalize(sparse, target_sum=1e4)
    got = sparse.layers["lognorm"]
    got = got.toarray() if sp.issparse(got) else got
    np.testing.assert_allclose(got, dense.layers["lognorm"], rtol=1e-5)


def test_normalize_records_provenance() -> None:
    adata = _toy()
    pp.pp.normalize(adata)
    hist = get_history(adata)
    assert hist[-1]["function"] == "pp.normalize"
    assert hist[-1]["params"]["log1p"] is True


def test_target_sum_none_uses_median() -> None:
    adata = _toy()
    pp.pp.normalize(adata, target_sum=None, log1p=False)
    # medians of [4, 10] = 7
    row_sums = np.asarray(adata.layers["lognorm"].sum(axis=1)).ravel()
    np.testing.assert_allclose(row_sums, [7.0, 7.0])
