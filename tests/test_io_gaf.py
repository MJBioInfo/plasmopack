"""Phase 2 tests — GAF parsing into gene sets.

Uses tests/data/sample_pberghei.gaf, a synthetic 8-annotation fixture covering
three GO terms across BP/MF/CC, a duplicate annotation, and one NOT-qualified
row that must be excluded.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plasmopack.io import read_gaf

_GAF = Path(__file__).parent / "data" / "sample_pberghei.gaf"


def test_groups_by_go_term() -> None:
    gs = read_gaf(_GAF)
    # three GO terms present (0006412 BP, 0003735 MF, 0005737 CC)
    assert set(gs.ids) == {"GO:0006412", "GO:0003735", "GO:0005737"}


def test_not_qualified_annotation_excluded() -> None:
    gs = read_gaf(_GAF)
    # PBANKA_050050 is annotated to GO:0006412 with NOT -> must be absent
    assert "PBANKA_050050" not in gs["GO:0006412"].genes


def test_duplicate_annotation_deduped() -> None:
    gs = read_gaf(_GAF)
    # PBANKA_010010 -> GO:0006412 appears twice in the file; count once
    assert gs["GO:0006412"].genes.count("PBANKA_010010") == 1
    assert set(gs["GO:0006412"].genes) == {
        "PBANKA_010010",
        "PBANKA_020020",
        "PBANKA_030030",
    }


def test_aspect_filter() -> None:
    bp = read_gaf(_GAF, aspect="BP")
    assert set(bp.ids) == {"GO:0006412"}

    mf = read_gaf(_GAF, aspect="MF")
    assert set(mf.ids) == {"GO:0003735"}


def test_aspect_label_stored_as_description() -> None:
    gs = read_gaf(_GAF)
    assert gs["GO:0006412"].description == "BP"
    assert gs["GO:0005737"].description == "CC"


def test_version_and_organism_recorded() -> None:
    gs = read_gaf(_GAF, version="PlasmoDB-68", organism="pberghei")
    assert gs.metadata["version"] == "PlasmoDB-68"
    assert gs.metadata["organism"] == "pberghei"


def test_invalid_aspect_raises() -> None:
    with pytest.raises(ValueError, match="aspect must be"):
        read_gaf(_GAF, aspect="XX")


def test_gaf_to_gmt_pipeline(tmp_path: Path) -> None:
    """The flagship path: GAF -> GeneSets -> filtered -> GMT -> read back."""
    from plasmopack.io import read_gmt, write_gmt

    gs = read_gaf(_GAF, version="PlasmoDB-68", organism="pberghei")
    gs = gs.filter_by_size(min_genes=2)  # drops the single-gene CC term
    path = write_gmt(gs, tmp_path / "pberghei_GO.gmt")

    back = read_gmt(path)
    assert set(back.ids) == {"GO:0006412", "GO:0003735"}
    assert back.metadata["version"] == "PlasmoDB-68"
