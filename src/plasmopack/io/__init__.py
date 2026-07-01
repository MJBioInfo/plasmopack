"""Input/output — reading and writing files.

Covers AnnData (``.h5ad``), GAF parsing, and GMT gene-set files (the flagship
gene-set builder). The Seurat bridge is documented separately. See Phase 2 in
``publish/ROADMAP.md``.

Public functions
----------------
- :func:`read_h5ad`, :func:`write_h5ad` — AnnData I/O
- :func:`read_gaf` — PlasmoDB/VEuPathDB GO annotations -> gene sets
- :func:`read_gmt`, :func:`write_gmt` — the standard gene-set exchange format
- :class:`GeneSet`, :class:`GeneSets` — the in-memory gene-set model
"""

from __future__ import annotations

from plasmopack.io._anndata import read_h5ad, write_h5ad
from plasmopack.io._gaf import read_gaf
from plasmopack.io._genesets import GeneSet, GeneSets
from plasmopack.io._gmt import read_gmt, write_gmt

__all__ = [
    "GeneSet",
    "GeneSets",
    "read_gaf",
    "read_gmt",
    "read_h5ad",
    "write_gmt",
    "write_h5ad",
]
