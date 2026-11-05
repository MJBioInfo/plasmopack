"""Preprocessing — QC, filtering, and normalization.

All functions operate in place on an :class:`~anndata.AnnData` (and also return
it, so calls can be chained) and stamp provenance into
``adata.uns["plasmopack"]["history"]``. See Phase 4 in ``publish/ROADMAP.md``.

Public functions
----------------
- :func:`qc_metrics` — per-cell QC into ``adata.obs``
- :func:`filter_cells`, :func:`filter_genes` — threshold filtering
- :func:`normalize` — library-size + log1p into a layer
- :func:`add_gene_symbols` — ``adata.var["symbol"]`` (symbol, else gene ID)
"""

from __future__ import annotations

from plasmopack.pp._normalize import normalize
from plasmopack.pp._qc import filter_cells, filter_genes, qc_metrics
from plasmopack.pp._symbols import add_gene_symbols

__all__ = [
    "add_gene_symbols",
    "filter_cells",
    "filter_genes",
    "normalize",
    "qc_metrics",
]
