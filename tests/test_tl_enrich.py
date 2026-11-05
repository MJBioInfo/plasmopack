"""Phase 7 tests — over-representation enrichment (ORA)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import plasmopack as pp
from plasmopack._utils.provenance import get_history
from plasmopack.io import GeneSet, GeneSets


def _sets() -> GeneSets:
    return GeneSets(
        sets=[
            GeneSet("PATH:1", "translation", [f"g{i}" for i in range(10)], "BP"),
            GeneSet("PATH:2", "glycolysis", [f"g{i}" for i in range(10, 20)], "BP"),
            GeneSet("PATH:small", "tiny", ["g0", "g1"], "BP"),  # below min_set
        ],
        metadata={"version": "TEST-1"},
    )


def _adata_with_background(n_genes: int = 100) -> AnnData:
    x = np.ones((5, n_genes), dtype=np.float32)
    adata = AnnData(X=x)
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    return adata


def test_explicit_query_enriches_right_set() -> None:
    adata = _adata_with_background()
    # query is exactly PATH:1 members -> PATH:1 should be the top hit
    pp.tl.enrich(
        adata,
        _sets(),
        genes=[f"g{i}" for i in range(10)],
        min_set=5,
        key_added="q",
    )
    res = adata.uns["enrichment"]["q"]
    assert res.iloc[0]["set_id"] == "PATH:1"
    assert res.iloc[0]["n_overlap"] == 10
    assert res.iloc[0]["padj"] < 0.05


def test_small_set_filtered_out() -> None:
    adata = _adata_with_background()
    pp.tl.enrich(adata, _sets(), genes=["g0", "g1"], min_set=5, key_added="q")
    res = adata.uns["enrichment"]["q"]
    assert "PATH:small" not in set(res["set_id"])


def test_no_query_args_raises() -> None:
    adata = _adata_with_background()
    with pytest.raises(ValueError, match="genes=.*or both groupby"):
        pp.tl.enrich(adata, _sets())


def test_provenance_records_gmt_version() -> None:
    adata = _adata_with_background()
    pp.tl.enrich(adata, _sets(), genes=["g0", "g1", "g2"], min_set=2, key_added="q")
    h = get_history(adata)[-1]
    assert h["function"] == "tl.enrich"
    # custom GeneSets -> collection label "custom" -> source key "custom_version"
    assert h["sources"]["custom_version"] == "TEST-1"


def test_from_deg_path() -> None:
    adata = _adata_with_background()
    # stub a DEG result: g0..g9 significant & up in group "A"
    deg_df = pd.DataFrame(
        {
            "group": "A",
            "gene": [f"g{i}" for i in range(20)],
            "logFC": [3.0] * 10 + [0.0] * 10,
            "pval": [1e-5] * 10 + [0.9] * 10,
            "padj": [1e-4] * 10 + [0.9] * 10,
            "method": "wilcoxon",
        }
    )
    adata.uns["deg"] = {"group": {"wilcoxon": deg_df}}

    pp.tl.enrich(
        adata,
        _sets(),
        groupby="group",
        group="A",
        direction="up",
        min_set=5,
    )
    res = adata.uns["enrichment"]["group:A:up"]
    # the 10 up genes are PATH:1 -> top hit
    assert res.iloc[0]["set_id"] == "PATH:1"


def test_top_n_ranked_selection() -> None:
    adata = _adata_with_background()
    # 20 genes; gene0..g9 strong, g10..g19 weak. top_n=5 by combined should
    # pick 5 of the strong ones -> PATH:1 still enriched.
    deg_df = pd.DataFrame(
        {
            "group": "A",
            "gene": [f"g{i}" for i in range(20)],
            "logFC": [5.0 - 0.1 * i for i in range(10)] + [0.2] * 10,
            "pval": [1e-8] * 10 + [0.5] * 10,
            "padj": [1e-7] * 10 + [0.5] * 10,
            "method": "wilcoxon",
        }
    )
    adata.uns["deg"] = {"group": {"wilcoxon": deg_df}}
    pp.tl.enrich(
        adata,
        _sets(),
        groupby="group",
        group="A",
        direction="up",
        top_n=5,
        rank_by="combined",
        min_set=5,
        key_added="topn",
    )
    res = adata.uns["enrichment"]["topn"]
    # 5 selected genes all in PATH:1 -> it's the hit
    assert res.iloc[0]["set_id"] == "PATH:1"
    assert res.iloc[0]["n_overlap"] == 5


def test_rank_by_invalid_raises() -> None:
    adata = _adata_with_background()
    adata.uns["deg"] = {
        "group": {
            "wilcoxon": pd.DataFrame(
                {
                    "group": ["A"],
                    "gene": ["g0"],
                    "logFC": [2.0],
                    "pval": [1e-3],
                    "padj": [1e-2],
                    "method": ["wilcoxon"],
                }
            )
        }
    }
    with pytest.raises(ValueError, match="rank_by must be"):
        pp.tl.enrich(
            adata,
            _sets(),
            groupby="group",
            group="A",
            top_n=1,
            rank_by="bogus",
            min_set=1,
        )


def test_missing_deg_result_errors() -> None:
    adata = _adata_with_background()
    with pytest.raises(KeyError, match="No DEG result"):
        pp.tl.enrich(adata, _sets(), groupby="group", group="A")


def test_reads_gmt_from_path(tmp_path) -> None:
    from plasmopack.io import write_gmt

    path = write_gmt(_sets(), tmp_path / "s.gmt")
    adata = _adata_with_background()
    pp.tl.enrich(
        adata, path, genes=[f"g{i}" for i in range(10)], min_set=5, key_added="q"
    )
    assert adata.uns["enrichment"]["q"].iloc[0]["set_id"] == "PATH:1"
