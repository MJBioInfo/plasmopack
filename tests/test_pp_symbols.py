"""Tests — gene symbols with ID fallback."""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData

import plasmopack as pp


def _adata() -> AnnData:
    a = AnnData(X=np.ones((2, 4), dtype=np.float32))
    a.var_names = ["PBANKA_010010", "PBANKA_020020", "PBANKA_030030", "PBANKA_040040"]
    return a


def test_from_col_with_fallback() -> None:
    a = _adata()
    a.var["Gene_symbol"] = ["MSP2", "", None, "AMA1"]  # 2 missing
    pp.pp.add_gene_symbols(a, from_col="Gene_symbol")
    assert list(a.var["symbol"]) == [
        "MSP2",
        "PBANKA_020020",  # empty -> id
        "PBANKA_030030",  # None -> id
        "AMA1",
    ]


def test_autodetect_symbol_column() -> None:
    a = _adata()
    a.var["gene_name"] = ["RPS26", "RPL6", "nan", "-"]
    pp.pp.add_gene_symbols(a)  # auto-detects gene_name
    assert list(a.var["symbol"]) == ["RPS26", "RPL6", "PBANKA_030030", "PBANKA_040040"]


def test_no_column_falls_back_to_ids() -> None:
    a = _adata()
    pp.pp.add_gene_symbols(a)
    assert list(a.var["symbol"]) == list(a.var_names)


def test_missing_from_col_raises() -> None:
    a = _adata()
    with pytest.raises(KeyError, match="not found in adata.var"):
        pp.pp.add_gene_symbols(a, from_col="nope")


def test_symbol_column_has_no_missing() -> None:
    a = _adata()
    a.var["symbol_merged"] = [np.nan, "X", np.nan, "Y"]
    pp.pp.add_gene_symbols(a)
    assert a.var["symbol"].notna().all()
    assert (a.var["symbol"].astype(str).str.len() > 0).all()
