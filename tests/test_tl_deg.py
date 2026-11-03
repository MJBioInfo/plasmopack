"""Phase 5 tests — differential expression with selectable methods."""

from __future__ import annotations

import importlib.util

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import plasmopack as pp
from plasmopack._utils.provenance import get_history

_HAS_SKLEARN = importlib.util.find_spec("sklearn") is not None
_HAS_PYDESEQ2 = importlib.util.find_spec("pydeseq2") is not None


def _toy(seed: int = 0) -> AnnData:
    """Two groups (A/B) x 2 replicates each; gene0 up in A, gene1 up in B."""
    rng = np.random.default_rng(seed)
    n_per = 40
    genes = 6
    blocks = []
    obs_group = []
    obs_sample = []
    for grp, s_ids in [("A", ["A1", "A2"]), ("B", ["B1", "B2"])]:
        for s in s_ids:
            base = rng.poisson(2.0, size=(n_per, genes)).astype(float)
            if grp == "A":
                base[:, 0] += rng.poisson(15.0, size=n_per)  # gene0 marker of A
            else:
                base[:, 1] += rng.poisson(15.0, size=n_per)  # gene1 marker of B
            blocks.append(base)
            obs_group += [grp] * n_per
            obs_sample += [s] * n_per
    x = np.vstack(blocks)
    adata = AnnData(X=x)
    adata.var_names = [f"gene{i}" for i in range(genes)]
    adata.obs["group"] = pd.Categorical(obs_group)
    adata.obs["sample"] = pd.Categorical(obs_sample)
    adata.layers["counts"] = x.copy()
    # a log-normalized layer for cell-level methods
    lib = x.sum(axis=1, keepdims=True)
    adata.layers["lognorm"] = np.log1p(x / np.where(lib > 0, lib, 1) * 1e4)
    return adata


# -- registry ---------------------------------------------------------------
def test_list_methods_table() -> None:
    tbl = pp.tl.list_deg_methods()
    assert "wilcoxon" in tbl.index
    assert bool(tbl.loc["wilcoxon", "available"]) is True
    assert bool(tbl.loc["mast", "available"]) is False
    # core methods need no extra
    assert tbl.loc["wilcoxon", "extra"] == "(core)"


def test_unknown_method_raises() -> None:
    adata = _toy()
    with pytest.raises(ValueError, match="Unknown DEG method"):
        pp.tl.deg(adata, "group", method="nope")


# -- wilcoxon (core) --------------------------------------------------------
def test_wilcoxon_detects_markers() -> None:
    adata = _toy()
    pp.tl.deg(adata, "group", method="wilcoxon")
    res = adata.uns["deg"]["group"]["wilcoxon"]
    a = res[res["group"] == "A"].set_index("gene")
    # gene0 up in A: positive logFC and significant
    assert a.loc["gene0", "logFC"] > 1
    assert a.loc["gene0", "padj"] < 0.05
    # gene1 (B's marker) should be down in A
    assert a.loc["gene1", "logFC"] < -1


def test_ttest_runs_and_stores() -> None:
    adata = _toy()
    pp.tl.deg(adata, "group", method="t-test")
    assert "t-test" in adata.uns["deg"]["group"]
    cols = set(adata.uns["deg"]["group"]["t-test"].columns)
    assert {"group", "gene", "logFC", "pval", "padj", "method"} <= cols


def test_provenance_recorded() -> None:
    adata = _toy()
    pp.tl.deg(adata, "group", method="wilcoxon")
    assert get_history(adata)[-1]["function"] == "tl.deg"


# -- pseudobulk (core, replicate-aware) -------------------------------------
def test_pseudobulk_requires_sample_key() -> None:
    adata = _toy()
    with pytest.raises(ValueError, match="requires\\s+sample_key"):
        pp.tl.deg(adata, "group", method="pseudobulk")


def test_pseudobulk_detects_markers() -> None:
    adata = _toy()
    pp.tl.deg(adata, "group", method="pseudobulk", sample_key="sample")
    res = adata.uns["deg"]["group"]["pseudobulk"]
    a = res[res["group"] == "A"].set_index("gene")
    assert a.loc["gene0", "logFC"] > 1  # gene0 up in A at the replicate level


def test_pseudobulk_too_few_replicates_errors() -> None:
    adata = _toy()
    # collapse to a single sample per group -> not enough replicates
    adata.obs["sample"] = adata.obs["group"].astype(str)
    with pytest.raises(ValueError, match="at least 2 replicates"):
        pp.tl.deg(adata, "group", method="pseudobulk", sample_key="sample")


# -- optional methods: clear errors when deps missing -----------------------
def test_mast_not_implemented() -> None:
    adata = _toy()
    with pytest.raises(NotImplementedError, match="mast"):
        pp.tl.deg(adata, "group", method="mast")


@pytest.mark.skipif(_HAS_SKLEARN, reason="sklearn installed")
def test_logreg_missing_dep_message() -> None:
    adata = _toy()
    with pytest.raises(ImportError, match=r"plasmopack\[ml\]"):
        pp.tl.deg(adata, "group", method="logreg")


@pytest.mark.skipif(_HAS_PYDESEQ2, reason="pydeseq2 installed")
def test_deseq2_missing_dep_message() -> None:
    adata = _toy()
    with pytest.raises(ImportError, match=r"plasmopack\[deg\]"):
        pp.tl.deg(adata, "group", method="deseq2", sample_key="sample")
