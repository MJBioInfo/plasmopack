"""Read and write GMT (Gene Matrix Transposed) files.

GMT is the standard gene-set exchange format used by GSEA, gseapy, decoupler,
and most enrichment tools. Each line is tab-separated::

    <set_id>\\t<description>\\t<gene1>\\t<gene2>\\t...

plasmopack additionally writes a JSON *sidecar* (``<path>.meta.json``) carrying
provenance (source, database version, organism, date). This sidecar is what
makes a GMT reproducible: it records exactly which annotation release produced
the sets. The GMT itself stays 100% standard so any other tool can read it.
"""

from __future__ import annotations

import json
from pathlib import Path

from plasmopack.io._genesets import GeneSet, GeneSets

_SIDECAR_SUFFIX = ".meta.json"


def _sidecar_path(path: Path) -> Path:
    return path.with_name(path.name + _SIDECAR_SUFFIX)


def write_gmt(
    gene_sets: GeneSets,
    path: str | Path,
    *,
    write_sidecar: bool = True,
) -> Path:
    """Write a :class:`GeneSets` to a GMT file (+ provenance sidecar).

    Parameters
    ----------
    gene_sets
        The sets to write.
    path
        Destination ``.gmt`` path.
    write_sidecar
        If True (default), also write ``<path>.meta.json`` with the
        ``gene_sets.metadata`` plus a summary (set count, gene-universe size).

    Returns
    -------
    Path
        The path the GMT was written to.
    """
    path = Path(path)
    lines: list[str] = []
    for s in gene_sets:
        # description column must never be empty in strict GMT readers; use id.
        description = s.description or s.name or s.id
        fields = [s.id, description, *s.genes]
        lines.append("\t".join(fields))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    if write_sidecar:
        meta = dict(gene_sets.metadata)
        meta["_summary"] = {
            "n_sets": len(gene_sets),
            "n_genes_total": len(gene_sets.gene_universe()),
            "gmt_file": path.name,
        }
        _sidecar_path(path).write_text(
            json.dumps(meta, indent=2, sort_keys=True), encoding="utf-8"
        )
    return path


def read_gmt(path: str | Path) -> GeneSets:
    """Read a GMT file into a :class:`GeneSets`.

    If a sidecar ``<path>.meta.json`` exists next to the file, its contents are
    loaded into ``GeneSets.metadata``; otherwise metadata records just the
    source path.

    Parameters
    ----------
    path
        Path to a ``.gmt`` file.

    Raises
    ------
    ValueError
        If a non-empty line has fewer than the minimum 2 tab-separated fields.
    """
    path = Path(path)
    sets: list[GeneSet] = []
    with path.open(encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(
                    f"{path}:{lineno}: GMT line needs at least 2 tab-separated "
                    f"fields (id, description), got {len(fields)}"
                )
            set_id, description, *genes = fields
            # Drop empty trailing gene cells (files sometimes have trailing tabs).
            genes = [g for g in genes if g]
            sets.append(
                GeneSet(
                    id=set_id,
                    name=set_id,
                    genes=genes,
                    description=description,
                )
            )

    sidecar = _sidecar_path(path)
    if sidecar.exists():
        metadata = json.loads(sidecar.read_text(encoding="utf-8"))
    else:
        metadata = {"source": str(path)}

    return GeneSets(sets=sets, metadata=metadata)
