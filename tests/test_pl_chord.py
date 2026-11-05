"""Phase 10 tests — pathway chord plot (smoke + error paths)."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless, no display needed

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from anndata import AnnData  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

import plasmopack as pp  # noqa: E402
from plasmopack.io import GeneSet, GeneSets  # noqa: E402


def _prepared() -> AnnData:
    rng = np.random.default_rng(0)
    n = 40
    x = rng.poisson(2.0, size=(60, n)).astype(float)
    x[:30, 0:8] += rng.poisson(12, size=(30, 8))
    adata = AnnData(X=x)
    adata.var_names = [f"g{i}" for i in range(n)]
    adata.obs["group"] = pd.Categorical(["A"] * 30 + ["B"] * 30)
    adata.layers["lognorm"] = np.log1p(x / x.sum(1, keepdims=True) * 1e4)
    pp.tl.deg(adata, "group", method="wilcoxon")
    sets = GeneSets(
        [
            GeneSet("P1", "pathway one", [f"g{i}" for i in range(8)]),
            GeneSet("P2", "pathway two", [f"g{i}" for i in range(5, 15)]),
        ],
        metadata={"version": "t"},
    )
    pp.tl.enrich(
        adata,
        sets,
        groupby="group",
        group="A",
        direction="both",
        min_set=3,
        key_added="A",
    )
    return adata


def test_chord_returns_figure() -> None:
    adata = _prepared()
    fig = pp.pl.pathway_chord(
        adata, "A", groupby="group", group="A", n_pathways=2, n_genes=8
    )
    assert isinstance(fig, Figure)


def test_chord_respects_n_genes() -> None:
    adata = _prepared()
    # should not raise with a small gene count
    fig = pp.pl.pathway_chord(
        adata, "A", groupby="group", group="A", n_pathways=2, n_genes=3
    )
    assert isinstance(fig, Figure)


def test_chord_missing_enrichment_errors() -> None:
    adata = _prepared()
    with pytest.raises(KeyError, match="No enrichment result"):
        pp.tl.pathway_chord_missing = None
        pp.pl.pathway_chord(adata, "NOPE", groupby="group", group="A")


def test_chord_missing_deg_errors() -> None:
    adata = _prepared()
    with pytest.raises(KeyError, match="No DEG result"):
        pp.pl.pathway_chord(adata, "A", groupby="group", group="A", deg_method="t-test")
