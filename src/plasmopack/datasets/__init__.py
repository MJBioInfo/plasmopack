"""Datasets — example data, the organism registry, and default gene sets.

See Phase 1 (organism registry) and Phase 11 (snapshots) in
``publish/ROADMAP.md``.
"""

from __future__ import annotations

from plasmopack.datasets._genesets import (
    build_gene_sets,
    detect_organism,
    gene_sets,
    list_bundled_gene_sets,
)
from plasmopack.datasets._registry import Organism, get_organism, list_organisms

__all__ = [
    "Organism",
    "build_gene_sets",
    "detect_organism",
    "gene_sets",
    "get_organism",
    "list_bundled_gene_sets",
    "list_organisms",
]
