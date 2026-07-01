"""Phase 3 tests — UniProt adapter and public client (network-free).

The transport is faked with the canned tests/data/uniprot_sample.json, so no
network is used. This is the durability pattern: real fetching is isolated
behind an injectable transport.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import plasmopack as pp
from plasmopack._utils.cache import DiskCache
from plasmopack._utils.http import HTTPError
from plasmopack.adapters import uniprot
from plasmopack.adapters.uniprot import parse_record

_SAMPLE = json.loads(
    (Path(__file__).parent / "data" / "uniprot_sample.json").read_text()
)


def _fake_transport(url: str) -> dict:
    return _SAMPLE


# -- parsing (pure) ---------------------------------------------------------
def test_parse_record_extracts_fields() -> None:
    rec = parse_record(_SAMPLE)
    assert rec.accession == "Q8I3H7"
    assert rec.entry_name == "SAMPLE_PLAF7"
    assert rec.protein_name == "Sample malaria protein"
    assert rec.gene_names == ["SAMP1"]
    assert rec.organism.startswith("Plasmodium falciparum")
    assert rec.taxon_id == 36329
    assert rec.length == 242
    assert rec.function == "Involved in a sample biological process."


def test_parse_record_falls_back_to_submission_name() -> None:
    # Typical TrEMBL (unreviewed) Plasmodium entry: no recommendedName, the
    # name lives under submissionNames (real UniProt schema, e.g. C6KSZ8).
    trembl = {
        "primaryAccession": "C6KSZ8",
        "proteinDescription": {
            "submissionNames": [
                {"fullName": {"value": "RNA and export factor binding protein"}}
            ]
        },
    }
    rec = parse_record(trembl)
    assert rec.protein_name == "RNA and export factor binding protein"


def test_parse_record_tolerates_missing_fields() -> None:
    rec = parse_record({"primaryAccession": "X1"})
    assert rec.accession == "X1"
    assert rec.protein_name == ""
    assert rec.gene_names == []
    assert rec.function is None
    assert rec.taxon_id is None


# -- fetch with injected transport -----------------------------------------
def test_fetch_uses_transport_not_network() -> None:
    recs = uniprot.fetch(["Q8I3H7"], transport=_fake_transport)
    assert len(recs) == 1
    assert recs[0].accession == "Q8I3H7"


def test_fetch_preserves_order() -> None:
    def transport(url: str) -> dict:
        acc = url.rsplit("/", 1)[-1].removesuffix(".json")
        return {"primaryAccession": acc}

    recs = uniprot.fetch(["A", "B", "C"], transport=transport)
    assert [r.accession for r in recs] == ["A", "B", "C"]


# -- caching behaviour ------------------------------------------------------
def test_cache_avoids_second_network_call(tmp_path: Path) -> None:
    calls = {"n": 0}

    def counting_transport(url: str) -> dict:
        calls["n"] += 1
        return _SAMPLE

    cache = DiskCache(tmp_path, namespace="uniprot")
    r1 = uniprot.fetch_one("Q8I3H7", transport=counting_transport, cache=cache)
    r2 = uniprot.fetch_one("Q8I3H7", transport=counting_transport, cache=cache)

    assert r1.accession == r2.accession == "Q8I3H7"
    assert calls["n"] == 1  # second call served from cache


# -- friendly error on gene-name-as-accession ------------------------------
def test_gene_name_as_accession_gives_helpful_error() -> None:
    def transport_400(url: str) -> dict:
        raise HTTPError("bad request", status=400)

    with pytest.raises(ValueError, match="not a valid UniProt accession"):
        uniprot.fetch_one("AMA1", transport=transport_400)


# -- search -----------------------------------------------------------------
def test_search_parses_results() -> None:
    def transport(url: str) -> dict:
        assert "search?" in url
        return {"results": [_SAMPLE, {"primaryAccession": "X2"}]}

    recs = uniprot.search("gene:SAMP1", transport=transport)
    assert [r.accession for r in recs] == ["Q8I3H7", "X2"]


def test_search_empty_results() -> None:
    recs = uniprot.search("gene:NOPE", transport=lambda url: {"results": []})
    assert recs == []


def test_public_search_builds_organism_query(monkeypatch, tmp_path) -> None:
    captured = {}

    def fake_adapter_search(query, *, limit, cache):  # type: ignore[no-untyped-def]
        captured["query"] = query
        captured["limit"] = limit
        return []

    monkeypatch.setattr("plasmopack.db.uniprot._adapter.search", fake_adapter_search)
    pp.db.uniprot.search("AMA1", organism="pfalciparum", limit=5)
    # UniProt uses the species taxon (5833) with subtree-matching taxonomy_id
    assert captured["query"] == "gene:AMA1 AND taxonomy_id:5833"
    assert captured["limit"] == 5


# -- public surface ---------------------------------------------------------
def test_public_fetch_accepts_single_string(tmp_path, monkeypatch) -> None:
    # point cache at tmp and stub the adapter's default transport
    monkeypatch.setenv("PLASMOPACK_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr("plasmopack.adapters.uniprot.http_get_json", _fake_transport)
    recs = pp.db.uniprot.fetch("Q8I3H7")
    assert isinstance(recs, list)
    assert recs[0].accession == "Q8I3H7"
