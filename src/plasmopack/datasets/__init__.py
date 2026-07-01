"""Datasets — example data, the organism registry, and frozen DB snapshots.

See Phase 1 (organism registry) and Phase 11 (snapshots) in
``publish/ROADMAP.md``.
"""

from __future__ import annotations

from plasmopack.datasets._registry import Organism, get_organism, list_organisms

__all__ = ["Organism", "get_organism", "list_organisms"]
