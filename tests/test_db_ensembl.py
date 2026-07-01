"""Phase 3 tests — Ensembl adapter and public client (network-free).

Uses the canned tests/data/ensembl_sample.json (a real P. falciparum MSP2
lookup) with an injected transport, so no network is used.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import plasmopack as pp
from plasmopack._utils.cache import DiskCache
from plasmopack._utils.http import HTTPError
from plasmopack.adapters import ensembl
from plasmopack.adapters.ensembl import parse_gene

_SAMPLE = json.loads(
    (Path(__file__).parent / "data" / "ensembl_sample.json").read_text()
)


def _fake_transport(url: str) -> dict:
    return _SAMPLE


# -- parsing (pure) ---------------------------------------------------------
def test_parse_gene_extracts_fields() -> None:
    g = parse_gene(_SAMPLE)
    assert g.id == "PF3D7_0206800"
    assert g.symbol == "MSP2"
    assert g.biotype == "protein_coding"
    assert g.species == "plasmodium_falciparum"
    assert g.chromosome == "Pf3D7_02_v3"
    assert g.start == 271576
    assert g.end == 274917
    assert g.strand == -1
    assert g.assembly == "GCA000002765v3"
    assert "Merozoite surface antigen 2" in g.description


def test_parse_gene_tolerates_missing_fields() -> None:
    g = parse_gene({"id": "X1"})
    assert g.id == "X1"
    assert g.symbol == ""
    assert g.start is None
    assert g.strand is None


# -- lookup with injected transport ----------------------------------------
def test_lookup_uses_transport_not_network() -> None:
    genes = ensembl.lookup(["PF3D7_0206800"], transport=_fake_transport)
    assert len(genes) == 1
    assert genes[0].symbol == "MSP2"


def test_lookup_preserves_order() -> None:
    def transport(url: str) -> dict:
        gid = url.rsplit("/", 1)[-1]
        return {"id": gid}

    genes = ensembl.lookup(["A", "B", "C"], transport=transport)
    assert [g.id for g in genes] == ["A", "B", "C"]


def test_unknown_id_gives_helpful_error() -> None:
    def transport_404(url: str) -> dict:
        raise HTTPError("not found", status=404)

    with pytest.raises(ValueError, match="no gene with id"):
        ensembl.lookup_one("BOGUS", transport=transport_404)


# -- caching ----------------------------------------------------------------
def test_cache_avoids_second_call(tmp_path: Path) -> None:
    calls = {"n": 0}

    def counting(url: str) -> dict:
        calls["n"] += 1
        return _SAMPLE

    cache = DiskCache(tmp_path, namespace="ensembl")
    ensembl.lookup_one("PF3D7_0206800", transport=counting, cache=cache)
    ensembl.lookup_one("PF3D7_0206800", transport=counting, cache=cache)
    assert calls["n"] == 1


# -- public surface ---------------------------------------------------------
def test_public_lookup_accepts_single_string(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("PLASMOPACK_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr("plasmopack.adapters.ensembl.http_get_json", _fake_transport)
    genes = pp.db.ensembl.lookup("PF3D7_0206800")
    assert isinstance(genes, list)
    assert genes[0].symbol == "MSP2"
