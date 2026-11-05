"""Tests — default gene sets: detection, parsing, bundled loading, resolution."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import plasmopack as pp
from plasmopack.datasets import (
    detect_organism,
    gene_sets,
    list_bundled_gene_sets,
)
from plasmopack.datasets._genesets import _parse_uniprot_tsv, _terms_to_genesets


# -- organism detection -----------------------------------------------------
def test_detect_pberghei() -> None:
    a = AnnData(X=np.ones((3, 4), dtype=np.float32))
    a.var_names = ["PBANKA_010010", "PBANKA_020020", "PBANKA_030030", "PBANKA_040040"]
    assert detect_organism(a) == "pberghei"


def test_detect_pfalciparum() -> None:
    a = AnnData(X=np.ones((2, 3), dtype=np.float32))
    a.var_names = ["PF3D7_0100100", "PF3D7_0200200", "PF3D7_0300300"]
    assert detect_organism(a) == "pfalciparum"


def test_detect_returns_none_for_unknown_ids() -> None:
    a = AnnData(X=np.ones((2, 2), dtype=np.float32))
    a.var_names = ["FOO1", "BAR2"]
    assert detect_organism(a) is None


# -- UniProt TSV parsing (pure) ---------------------------------------------
def test_parse_uniprot_tsv() -> None:
    rows = [
        "A1\tPlasmoDB:PBANKA_010010;\ttranslation [GO:0006412]\tcytosol [GO:0005829]\t",
        "A2\tPlasmoDB:PBANKA_020020;\ttranslation [GO:0006412]\t\tbinding [GO:0005488]",
        "A3\tsomethingelse\ttranslation [GO:0006412]\t\t",  # no PBANKA -> skipped
    ]
    terms = _parse_uniprot_tsv(rows, "PBANKA_")
    assert set(terms["GO:0006412"][2]) == {"PBANKA_010010", "PBANKA_020020"}
    assert terms["GO:0006412"][0] == "translation"
    assert terms["GO:0006412"][1] == "BP"
    assert terms["GO:0005488"][1] == "MF"  # from go_f column


def test_terms_to_genesets_size_filter() -> None:
    terms = {
        "GO:1": ("big", "BP", {"g1", "g2", "g3"}),
        "GO:2": ("small", "MF", {"g1"}),
    }
    gs = _terms_to_genesets(
        terms, organism="pberghei", version="v", min_set=2, max_set=500
    )
    assert gs.ids == ["GO:1"]
    assert gs.metadata["organism"] == "pberghei"


# -- bundled snapshots ------------------------------------------------------
def test_four_organisms_bundled() -> None:
    bundled = set(list_bundled_gene_sets())
    assert {"pberghei", "pfalciparum", "pknowlesi", "tgondii"} <= bundled


def test_load_bundled_pberghei() -> None:
    gs = gene_sets("pberghei", "GO")
    assert len(gs) > 100
    assert gs.metadata["organism"] == "pberghei"
    # a well-known ribosome term should be present
    assert "GO:0006412" in gs  # translation


def test_aspect_filter() -> None:
    allgo = gene_sets("pberghei", "GO")
    bp = gene_sets("pberghei", "BP")
    assert 0 < len(bp) < len(allgo)
    assert all(s.description == "BP" for s in bp)


def test_bad_collection_raises() -> None:
    with pytest.raises(ValueError, match="collection must be"):
        gene_sets("pberghei", "NONSENSE")


def test_kegg_bundled_and_named() -> None:
    kegg = gene_sets("pberghei", "KEGG")
    assert len(kegg) > 20
    assert kegg.metadata["collection"] == "KEGG"
    # readable pathway names survive (not just KEGG ids)
    assert any("metabol" in s.name.lower() or "glycol" in s.name.lower() for s in kegg)


def test_unbundled_organism_raises() -> None:
    with pytest.raises(FileNotFoundError, match="No bundled"):
        gene_sets("pvivax", "GO")


# -- enrich resolves the string -> bundled sets -----------------------------
def test_enrich_with_go_string_autodetect() -> None:
    # tiny adata with real PBANKA_ ids (auto-detected), stubbed DEG
    genes = sorted(gene_sets("pberghei", "GO")["GO:0006412"].genes)[:30]
    a = AnnData(X=np.ones((4, len(genes)), dtype=np.float32))
    a.var_names = genes
    a.uns["deg"] = {
        "time": {
            "wilcoxon": pd.DataFrame(
                {
                    "group": "24h",
                    "gene": genes,
                    "logFC": [3.0] * len(genes),
                    "pval": [1e-6] * len(genes),
                    "padj": [1e-5] * len(genes),
                    "method": "wilcoxon",
                }
            )
        }
    }
    pp.tl.enrich(a, "GO", groupby="time", group="24h", direction="up", key_added="e")
    res = a.uns["enrichment"]["e"]
    # translation term should be a top hit (its genes are the query)
    assert "GO:0006412" in set(res["set_id"])


def test_enrich_multi_collection_one_table() -> None:
    # GO + KEGG in one call -> single table with a 'collection' column
    go = gene_sets("pberghei", "GO")
    genes = sorted(go["GO:0006412"].genes)[:30]
    a = AnnData(X=np.ones((4, len(genes)), dtype=np.float32))
    a.var_names = genes
    a.uns["deg"] = {
        "time": {
            "wilcoxon": pd.DataFrame(
                {
                    "group": "24h",
                    "gene": genes,
                    "logFC": [3.0] * len(genes),
                    "pval": [1e-6] * len(genes),
                    "padj": [1e-5] * len(genes),
                    "method": "wilcoxon",
                }
            )
        }
    }
    pp.tl.enrich(
        a,
        ["GO", "KEGG"],
        groupby="time",
        group="24h",
        direction="up",
        key_added="multi",
    )
    res = a.uns["enrichment"]["multi"]
    assert "collection" in res.columns
    assert set(res["collection"].unique()) <= {"GO", "KEGG"}
    assert "GO" in set(res["collection"])


def test_enrich_go_string_no_organism_detectable() -> None:
    a = AnnData(X=np.ones((2, 2), dtype=np.float32))
    a.var_names = ["FOO1", "BAR2"]
    a.uns["deg"] = {
        "g": {
            "wilcoxon": pd.DataFrame(
                {
                    "group": ["A"],
                    "gene": ["FOO1"],
                    "logFC": [2.0],
                    "pval": [1e-3],
                    "padj": [1e-2],
                    "method": ["wilcoxon"],
                }
            )
        }
    }
    with pytest.raises(ValueError, match="could not be detected"):
        pp.tl.enrich(a, "GO", groupby="g", group="A")
