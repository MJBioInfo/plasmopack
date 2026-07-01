"""Phase 2 tests — AnnData read/write wrappers."""

from __future__ import annotations

from pathlib import Path

from anndata import AnnData

from plasmopack.io import read_h5ad, write_h5ad


def test_h5ad_roundtrip(dummy_pberghei: AnnData, tmp_path: Path) -> None:
    path = write_h5ad(dummy_pberghei, tmp_path / "x.h5ad")
    assert path.exists()

    back = read_h5ad(path)
    assert back.shape == dummy_pberghei.shape
    assert list(back.var["gene_id"]) == list(dummy_pberghei.var["gene_id"])
