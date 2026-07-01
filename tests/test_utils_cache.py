"""Phase 3 tests — the on-disk cache and cache-key helper."""

from __future__ import annotations

from pathlib import Path

from plasmopack._utils.cache import DiskCache, make_key


def test_make_key_is_stable_and_distinct() -> None:
    assert make_key("uniprot", "Q1") == make_key("uniprot", "Q1")
    assert make_key("uniprot", "Q1") != make_key("uniprot", "Q2")


def test_set_get_roundtrip(tmp_path: Path) -> None:
    cache = DiskCache(tmp_path, namespace="uniprot")
    assert cache.get("missing") is None

    cache.set("k", {"a": 1, "b": [2, 3]})
    assert cache.get("k") == {"a": 1, "b": [2, 3]}


def test_namespace_isolation(tmp_path: Path) -> None:
    a = DiskCache(tmp_path, namespace="uniprot")
    b = DiskCache(tmp_path, namespace="ensembl")
    a.set("k", {"who": "uniprot"})
    assert b.get("k") is None


def test_clear(tmp_path: Path) -> None:
    cache = DiskCache(tmp_path, namespace="x")
    cache.set("k1", {"v": 1})
    cache.set("k2", {"v": 2})
    cache.clear()
    assert cache.get("k1") is None
    assert cache.get("k2") is None


def test_env_var_override(tmp_path: Path, monkeypatch) -> None:
    from plasmopack._utils import cache as cache_mod

    monkeypatch.setenv("PLASMOPACK_CACHE_DIR", str(tmp_path / "custom"))
    assert cache_mod.default_cache_dir() == tmp_path / "custom"
