"""The ``GeneSets`` in-memory object.

A gene set is a named collection of genes (e.g. a GO term, a pathway, an
ortholog group). ``GeneSets`` is the working unit that plasmopack's I/O and
enrichment code passes around: it holds the sets *plus* provenance metadata
(where they came from, which database release, when), so a set built today is
still interpretable years later.

The object is deliberately plain — dataclasses over numpy/pandas only — so it
has no heavy dependencies and round-trips cleanly to GMT + JSON on disk.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Any


@dataclass
class GeneSet:
    """A single named set of genes.

    Attributes
    ----------
    id
        Stable identifier, e.g. ``"GO:0006412"``.
    name
        Human-readable label, e.g. ``"translation"``. Falls back to ``id``
        when no name is available (GAF files carry IDs, not term names).
    genes
        Member gene identifiers. Order is preserved but duplicates are removed
        on construction.
    description
        Optional free text (the GMT "description" column). Often the GO aspect
        (``"BP"``/``"CC"``/``"MF"``) or a source note.
    """

    id: str
    name: str
    genes: list[str]
    description: str = ""

    def __post_init__(self) -> None:
        # De-duplicate while preserving first-seen order.
        seen: dict[str, None] = {}
        for g in self.genes:
            seen.setdefault(g, None)
        self.genes = list(seen)

    def __len__(self) -> int:
        return len(self.genes)


@dataclass
class GeneSets:
    """An ordered collection of :class:`GeneSet` plus provenance metadata.

    Parameters
    ----------
    sets
        The gene sets.
    metadata
        Provenance dictionary. Conventional keys: ``source`` (e.g.
        ``"PlasmoDB GOA"``), ``version`` (e.g. ``"PlasmoDB-68"``),
        ``organism``, ``date``, ``id_namespace``. Kept JSON-serialisable.
    """

    sets: list[GeneSet] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self._index: dict[str, GeneSet] = {s.id: s for s in self.sets}

    # -- container protocol -------------------------------------------------
    def __len__(self) -> int:
        return len(self.sets)

    def __iter__(self) -> Iterator[GeneSet]:
        return iter(self.sets)

    def __contains__(self, set_id: str) -> bool:
        return set_id in self._index

    def __getitem__(self, set_id: str) -> GeneSet:
        return self._index[set_id]

    # -- convenience --------------------------------------------------------
    @property
    def ids(self) -> list[str]:
        """The set identifiers, in order."""
        return [s.id for s in self.sets]

    def gene_universe(self) -> set[str]:
        """All unique genes appearing across every set."""
        universe: set[str] = set()
        for s in self.sets:
            universe.update(s.genes)
        return universe

    def sizes(self) -> dict[str, int]:
        """Map each set id to its gene count."""
        return {s.id: len(s) for s in self.sets}

    def filter_by_size(
        self, min_genes: int = 1, max_genes: int | None = None
    ) -> GeneSets:
        """Return a new ``GeneSets`` keeping only sets within a size range.

        Standard GSEA practice drops very small and very large sets
        (defaults there are typically ``min_genes=5``, ``max_genes=500``).

        Parameters
        ----------
        min_genes
            Minimum number of genes (inclusive).
        max_genes
            Maximum number of genes (inclusive). ``None`` means no upper bound.
        """
        kept = [
            s
            for s in self.sets
            if len(s) >= min_genes and (max_genes is None or len(s) <= max_genes)
        ]
        new_meta = dict(self.metadata)
        new_meta["filtered"] = {"min_genes": min_genes, "max_genes": max_genes}
        return GeneSets(sets=kept, metadata=new_meta)
