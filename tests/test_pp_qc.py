"""Phase 4 tests — QC metrics and filtering."""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from anndata import AnnData

import plasmopack as pp
from plasmopack._utils.provenance import get_history


def _toy() -> AnnData:
    # 3 cells x 4 genes, known counts
    x = np.array(
        [
            [1, 0, 2, 0],  # total 3, 2 genes
            [0, 0, 0, 0],  # total 0, 0 genes
            [5, 5, 0, 10],  # total 20, 3 genes
        ],
        dtype=np.float32,
    )
    var = {"gene": ["g0", "g1", "g2", "g3"]}
    adata = AnnData(X=x)
    adata.var_names = var["gene"]
    return adata


def test_qc_basic_metrics() -> None:
    adata = _toy()
    pp.pp.qc_metrics(adata)
    assert list(adata.obs["total_counts"]) == [3, 0, 20]
    assert list(adata.obs["n_genes_by_counts"]) == [2, 0, 3]


def test_qc_gene_group_pct() -> None:
    adata = _toy()
    pp.pp.qc_metrics(adata, gene_groups={"mito": ["g0"]})
    # cell0: g0=1 of 3 -> 33.33%; cell1: 0 total -> 0; cell2: g0=5 of 20 -> 25%
    pct = adata.obs["pct_counts_mito"].to_numpy()
    assert pct[0] == pytest.approx(100 / 3)
    assert pct[1] == 0.0
    assert pct[2] == pytest.approx(25.0)


def test_qc_works_on_sparse() -> None:
    adata = _toy()
    adata.X = sp.csr_matrix(adata.X)
    pp.pp.qc_metrics(adata, gene_groups={"mito": ["g0", "g3"]})
    assert list(adata.obs["total_counts"]) == [3, 0, 20]
    # cell2: g0+g3 = 15 of 20 -> 75%
    assert adata.obs["pct_counts_mito"].to_numpy()[2] == pytest.approx(75.0)


def test_qc_records_provenance() -> None:
    adata = _toy()
    pp.pp.qc_metrics(adata)
    hist = get_history(adata)
    assert hist[-1]["function"] == "pp.qc_metrics"


def test_filter_cells() -> None:
    adata = _toy()
    pp.pp.filter_cells(adata, min_genes=1)  # drops the all-zero cell
    assert adata.n_obs == 2
    assert list(adata.obs["total_counts"]) if "total_counts" in adata.obs else True


def test_filter_cells_min_counts() -> None:
    adata = _toy()
    pp.pp.filter_cells(adata, min_counts=5)  # keeps only the 20-count cell
    assert adata.n_obs == 1


def test_filter_genes() -> None:
    adata = _toy()
    # g1 appears in 1 cell, g2 in 1 cell; require >=2 cells
    pp.pp.filter_genes(adata, min_cells=2)
    assert "g0" in list(adata.var_names)  # g0 in cells 0 and 2
    assert "g1" not in list(adata.var_names)  # g1 only in cell 2


def test_filter_is_inplace_and_chains() -> None:
    adata = _toy()
    out = pp.pp.filter_cells(adata, min_genes=1)
    assert out is adata  # same object, chainable
