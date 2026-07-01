"""plasmopack — AnnData-native toolkit for Plasmodium and VEuPathDB organisms.

The public API is organised into scanpy-style submodules::

    import plasmopack as pp

    pp.io    # read/write: AnnData, GFF, GAF, GMT, Seurat bridge
    pp.db    # database clients: VEuPathDB, UniProt, Ensembl, NCBI, ...
    pp.pp    # preprocessing: QC, normalization, ID normalization
    pp.tl    # tools: DEG, compositional, enrichment, TF, orthologs
    pp.pl    # plotting
    pp.datasets  # example data + organism registry

Every function that annotates an ``AnnData`` writes its results back into the
object using a documented slot convention and stamps provenance into
``adata.uns["plasmopack"]["history"]``.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from plasmopack import datasets, db, io, pl, pp, tl

try:
    __version__ = version("plasmopack")
except PackageNotFoundError:  # package not installed (e.g. running from source tree)
    __version__ = "0.0.0"

__all__ = [
    "__version__",
    "datasets",
    "db",
    "io",
    "pl",
    "pp",
    "tl",
]
