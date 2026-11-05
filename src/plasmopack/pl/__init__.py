"""Plotting — matplotlib-based visualizations.

Every function returns a matplotlib ``Figure``. See Phase 10 in
``publish/ROADMAP.md``.

Public functions
----------------
- :func:`pathway_chord` — circos-style pathway <-> gene chord diagram,
  genes coloured by fold change
"""

from __future__ import annotations

from plasmopack.pl._chord import pathway_chord

__all__ = ["pathway_chord"]
