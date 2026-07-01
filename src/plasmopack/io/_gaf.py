"""Read GAF (Gene Association Format) files into gene sets.

GAF is the standard Gene Ontology annotation format published by PlasmoDB /
VEuPathDB (and the GO Consortium) for every organism and release. It is
gene-centric: one row per (gene, GO term) annotation. plasmopack inverts it
into pathway-centric :class:`GeneSets` (one entry per GO term with its member
genes) — the form enrichment tools need.

GAF 2.x is tab-separated with 17 columns; lines beginning with ``!`` are
comments. The columns used here:

===  ======================  =========================================
Col  Name                    Use
===  ======================  =========================================
2    DB Object ID            gene identifier
4    Qualifier               skip rows containing ``NOT``
5    GO ID                   the gene-set id
9    Aspect                  P/F/C -> BP/MF/CC (optional filter)
===  ======================  =========================================
"""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

from plasmopack.io._genesets import GeneSet, GeneSets

# GAF aspect code -> GO namespace label.
_ASPECT_LABEL = {"P": "BP", "F": "MF", "C": "CC"}

# Zero-based column indices in a GAF 2.x row.
_COL_OBJECT_ID = 1
_COL_QUALIFIER = 3
_COL_GO_ID = 4
_COL_ASPECT = 8
_MIN_COLS = 9


def read_gaf(
    path: str | Path,
    *,
    aspect: str | None = None,
    version: str | None = None,
    organism: str | None = None,
) -> GeneSets:
    """Parse a GAF file into a :class:`GeneSets` grouped by GO term.

    Parameters
    ----------
    path
        Path to a ``.gaf`` (or ``.gaf.gz`` — not yet; plain text for now) file.
    aspect
        If given, keep only one namespace: ``"BP"``, ``"CC"``, or ``"MF"``.
        ``None`` keeps all three.
    version
        Optional database release label (e.g. ``"PlasmoDB-68"``) recorded in
        ``GeneSets.metadata["version"]``. This is the reproducibility anchor —
        supply it whenever you know it.
    organism
        Optional organism key recorded in metadata.

    Returns
    -------
    GeneSets
        One set per GO term. ``NOT``-qualified annotations are excluded.

    Raises
    ------
    ValueError
        If ``aspect`` is not one of BP/CC/MF/None.
    """
    if aspect is not None and aspect not in {"BP", "CC", "MF"}:
        raise ValueError(f"aspect must be BP, CC, MF, or None; got {aspect!r}")

    path = Path(path)
    # Preserve first-seen order of GO terms for deterministic output.
    term_genes: OrderedDict[str, list[str]] = OrderedDict()
    term_aspect: dict[str, str] = {}

    with path.open(encoding="utf-8") as fh:
        for raw in fh:
            if not raw or raw.startswith("!"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < _MIN_COLS:
                continue

            qualifier = cols[_COL_QUALIFIER]
            if "NOT" in qualifier.upper().split("|"):
                continue

            aspect_code = cols[_COL_ASPECT]
            label = _ASPECT_LABEL.get(aspect_code, aspect_code)
            if aspect is not None and label != aspect:
                continue

            go_id = cols[_COL_GO_ID]
            gene = cols[_COL_OBJECT_ID]
            if not go_id or not gene:
                continue

            term_genes.setdefault(go_id, []).append(gene)
            term_aspect.setdefault(go_id, label)

    sets = [
        GeneSet(
            id=go_id,
            name=go_id,  # GAF carries no term names; enrich later via .obo
            genes=genes,
            description=term_aspect.get(go_id, ""),
        )
        for go_id, genes in term_genes.items()
    ]

    metadata: dict[str, str | int] = {
        "source": "GAF",
        "source_file": path.name,
    }
    if version is not None:
        metadata["version"] = version
    if organism is not None:
        metadata["organism"] = organism
    if aspect is not None:
        metadata["aspect"] = aspect

    return GeneSets(sets=sets, metadata=metadata)
