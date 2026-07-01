"""Phase 2 tests — the GeneSets in-memory model."""

from __future__ import annotations

from plasmopack.io import GeneSet, GeneSets


def test_geneset_dedupes_preserving_order() -> None:
    s = GeneSet(id="GO:1", name="one", genes=["a", "b", "a", "c", "b"])
    assert s.genes == ["a", "b", "c"]
    assert len(s) == 3


def test_genesets_container_protocol() -> None:
    gs = GeneSets(
        sets=[
            GeneSet("GO:1", "one", ["a", "b"]),
            GeneSet("GO:2", "two", ["b", "c", "d"]),
        ]
    )
    assert len(gs) == 2
    assert "GO:1" in gs
    assert "GO:9" not in gs
    assert gs["GO:2"].genes == ["b", "c", "d"]
    assert gs.ids == ["GO:1", "GO:2"]
    assert [s.id for s in gs] == ["GO:1", "GO:2"]


def test_gene_universe_and_sizes() -> None:
    gs = GeneSets(
        sets=[
            GeneSet("GO:1", "one", ["a", "b"]),
            GeneSet("GO:2", "two", ["b", "c", "d"]),
        ]
    )
    assert gs.gene_universe() == {"a", "b", "c", "d"}
    assert gs.sizes() == {"GO:1": 2, "GO:2": 3}


def test_filter_by_size() -> None:
    gs = GeneSets(
        sets=[
            GeneSet("small", "s", ["a"]),
            GeneSet("mid", "m", ["a", "b", "c"]),
            GeneSet("big", "b", ["a", "b", "c", "d", "e"]),
        ]
    )
    filtered = gs.filter_by_size(min_genes=2, max_genes=4)
    assert filtered.ids == ["mid"]
    # original is untouched
    assert gs.ids == ["small", "mid", "big"]
    # filter parameters recorded in metadata
    assert filtered.metadata["filtered"] == {"min_genes": 2, "max_genes": 4}
