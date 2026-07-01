"""Phase 1 tests — provenance stamping."""

from __future__ import annotations

from anndata import AnnData

from plasmopack._utils.provenance import (
    HISTORY_KEY,
    UNS_ROOT,
    get_history,
    record_provenance,
)


def test_record_creates_history(dummy_pberghei: AnnData) -> None:
    rec = record_provenance(
        dummy_pberghei,
        "tl.enrich",
        params={"groupby": "time", "method": "local_ora"},
        sources={"gmt_version": "PlasmoDB-68"},
    )
    assert rec["function"] == "tl.enrich"
    assert rec["params"]["groupby"] == "time"
    assert rec["sources"]["gmt_version"] == "PlasmoDB-68"
    assert "timestamp" in rec
    assert "plasmopack_version" in rec

    hist = get_history(dummy_pberghei)
    assert len(hist) == 1
    assert hist[0]["function"] == "tl.enrich"


def test_multiple_records_accumulate(dummy_pberghei: AnnData) -> None:
    record_provenance(dummy_pberghei, "pp.qc_metrics")
    record_provenance(dummy_pberghei, "tl.deg", params={"groupby": "time"})
    record_provenance(dummy_pberghei, "tl.enrich")

    hist = get_history(dummy_pberghei)
    assert [r["function"] for r in hist] == [
        "pp.qc_metrics",
        "tl.deg",
        "tl.enrich",
    ]


def test_history_survives_h5ad_roundtrip(dummy_pberghei: AnnData, tmp_path) -> None:
    import anndata as ad

    record_provenance(dummy_pberghei, "tl.deg", params={"method": "wilcoxon"})
    path = tmp_path / "roundtrip.h5ad"
    dummy_pberghei.write(path)

    reloaded = ad.read_h5ad(path)
    hist = get_history(reloaded)
    assert len(hist) == 1
    assert hist[0]["function"] == "tl.deg"
    assert hist[0]["params"]["method"] == "wilcoxon"


def test_uns_root_structure(dummy_pberghei: AnnData) -> None:
    record_provenance(dummy_pberghei, "pp.qc_metrics")
    assert UNS_ROOT in dummy_pberghei.uns
    assert HISTORY_KEY in dummy_pberghei.uns[UNS_ROOT]
