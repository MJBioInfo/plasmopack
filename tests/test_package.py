"""Phase 0 smoke tests — the package imports and exposes its namespaces."""

from __future__ import annotations

import plasmopack as pp


def test_version_is_string() -> None:
    assert isinstance(pp.__version__, str)
    assert pp.__version__


def test_namespaces_present() -> None:
    for name in ("io", "db", "pp", "tl", "pl", "datasets"):
        assert hasattr(pp, name), f"missing namespace pp.{name}"


def test_cci_submodule_present() -> None:
    assert hasattr(pp.tl, "cci")
