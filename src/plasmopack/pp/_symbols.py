"""Gene symbols with ID fallback.

Populates ``adata.var["symbol"]`` with a human-readable gene symbol (e.g.
``MSP2``) wherever one is available, and falls back to the gene ID (e.g.
``PBANKA_0206800``) where it is not — so plots and tables always show something
meaningful. Symbols come from an existing ``var`` column (auto-detected) or,
optionally, from Ensembl.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from plasmopack._utils.provenance import record_provenance

# var columns that commonly hold a gene symbol, in order of preference.
_SYMBOL_HINTS = (
    "symbol",
    "gene_symbol",
    "genesymbol",
    "gene_name",
    "genename",
    "display_name",
    "gene_label",
    "symbol_merged",
)


def _looks_like_symbol_col(name: str) -> bool:
    n = name.lower().replace(" ", "_")
    return any(h in n for h in _SYMBOL_HINTS)


def _detect_symbol_column(adata: AnnData) -> str | None:
    lowered = {str(c).lower().replace(" ", "_"): str(c) for c in adata.var.columns}
    for hint in _SYMBOL_HINTS:
        if hint in lowered:
            return lowered[hint]
    for c in adata.var.columns:
        if _looks_like_symbol_col(str(c)):
            return str(c)
    return None


def _clean(value: object, fallback: str) -> str:
    """Return a usable symbol, else the fallback ID."""
    if value is None:
        return fallback
    if isinstance(value, float) and np.isnan(value):
        return fallback
    s = str(value).strip()
    # treat empties and ID-looking values (same as the gene id) as 'no symbol'
    if not s or s.lower() in {"nan", "none", "na", "-"}:
        return fallback
    return s


def add_gene_symbols(
    adata: AnnData,
    *,
    from_col: str | None = None,
    fetch: bool = False,
    key_added: str = "symbol",
) -> AnnData:
    """Add ``adata.var[key_added]`` = gene symbol, falling back to the gene ID.

    Parameters
    ----------
    adata
        Annotated data; modified in place and returned.
    from_col
        A ``var`` column holding symbols. If ``None``, a likely column is
        auto-detected (``Gene_symbol``, ``gene_name``, ``display_name``, …).
    fetch
        If no suitable column is found (or to override), look symbols up from
        Ensembl via :func:`plasmopack.db.ensembl.lookup`. Network + slower; the
        results are cached. Genes without an Ensembl symbol keep their ID.
    key_added
        Destination ``var`` column (default ``"symbol"``).

    Returns
    -------
    AnnData
        The same object, for chaining. ``adata.var[key_added]`` has no missing
        values — every gene has a symbol or its ID.
    """
    ids = [str(v) for v in adata.var_names]
    symbols: list[str]

    col = from_col or (None if fetch else _detect_symbol_column(adata))
    source: str
    if col is not None:
        if col not in adata.var.columns:
            raise KeyError(f"from_col={col!r} not found in adata.var.")
        symbols = [_clean(v, gid) for v, gid in zip(adata.var[col], ids, strict=True)]
        source = f"var:{col}"
    elif fetch:
        from plasmopack.db import ensembl

        recs = ensembl.lookup(ids)
        mapping = {r.id: r.symbol for r in recs}
        symbols = [_clean(mapping.get(gid), gid) for gid in ids]
        source = "ensembl"
    else:
        # nothing to go on — everything falls back to the ID
        symbols = list(ids)
        source = "id_only"

    adata.var[key_added] = pd.Series(symbols, index=adata.var_names)
    record_provenance(
        adata,
        "pp.add_gene_symbols",
        params={
            "from_col": col,
            "fetch": fetch,
            "key_added": key_added,
            "source": source,
        },
    )
    return adata
