"""Default gene sets: organism detection, bundled snapshots, and live build.

Lets pathway enrichment work with zero setup::

    pp.tl.enrich(adata, gene_sets="GO")   # organism auto-detected from gene IDs

The organism is resolved from the gene-ID pattern in ``adata.var_names`` (via the
registry), the collection ("GO"/"BP"/"CC"/"MF") picks which sets, and a bundled,
version-stamped snapshot supplies the gene sets offline. :func:`build_gene_sets`
fetches the latest annotations from UniProt (GO terms + PlasmoDB/VEuPathDB gene
IDs) for when a user wants current data.
"""

from __future__ import annotations

import re
import urllib.parse
import urllib.request
from collections import defaultdict
from importlib.resources import files
from pathlib import Path

from anndata import AnnData

from plasmopack.datasets._registry import get_organism, list_organisms
from plasmopack.io import GeneSet, GeneSets, read_gmt, write_gmt

_ASPECT_FIELDS = {"BP": "go_p", "CC": "go_c", "MF": "go_f"}
# "name [GO:0006412]" -> ("GO:0006412", "name")
_GO_RE = re.compile(r"^(.*?)\s*\[(GO:\d+)\]$")
_UNIPROT = "https://rest.uniprot.org/uniprotkb/search"


# ---------------------------------------------------------------------------
# Organism detection from gene IDs
# ---------------------------------------------------------------------------
def detect_organism(adata: AnnData, *, min_fraction: float = 0.3) -> str | None:
    """Guess the organism from the gene-ID pattern in ``adata.var_names``.

    Returns the registry key (e.g. ``"pberghei"``) of the organism whose gene-ID
    regex matches the most var_names, provided at least ``min_fraction`` of genes
    match; otherwise ``None``.
    """
    names = [str(v) for v in adata.var_names]
    if not names:
        return None
    best: str | None = None
    best_frac = 0.0
    for key in list_organisms():
        org = get_organism(key)
        matched = sum(1 for n in names if org.matches_gene_id(n))
        frac = matched / len(names)
        if frac > best_frac:
            best, best_frac = key, frac
    return best if best_frac >= min_fraction else None


# ---------------------------------------------------------------------------
# Parsing UniProt TSV into GeneSets (pure, testable)
# ---------------------------------------------------------------------------
def _parse_uniprot_tsv(
    rows: list[str], gene_prefix: str
) -> dict[str, tuple[str, str, set[str]]]:
    """Parse UniProt TSV rows into ``{go_id: (name, aspect, {genes})}``.

    Columns expected: accession, xref_veupathdb, go_p, go_c, go_f.
    """
    pb_re = re.compile(rf"{re.escape(gene_prefix)}\w+")
    terms: dict[str, tuple[str, str, set[str]]] = {}
    gene_members: dict[str, set[str]] = defaultdict(set)
    term_meta: dict[str, tuple[str, str]] = {}
    for line in rows:
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        _, xref, go_p, go_c, go_f = parts[:5]
        m = pb_re.search(xref)
        if not m:
            continue
        gene = m.group(0)
        for aspect, field in zip(("BP", "CC", "MF"), (go_p, go_c, go_f), strict=True):
            for entry in field.split(";"):
                entry = entry.strip()
                if not entry:
                    continue
                gm = _GO_RE.match(entry)
                if not gm:
                    continue
                name, go_id = gm.group(1).strip(), gm.group(2)
                gene_members[go_id].add(gene)
                term_meta.setdefault(go_id, (name, aspect))
    for go_id, genes in gene_members.items():
        name, aspect = term_meta[go_id]
        terms[go_id] = (name, aspect, genes)
    return terms


def _terms_to_genesets(
    terms: dict[str, tuple[str, str, set[str]]],
    *,
    organism: str,
    version: str,
    min_set: int,
    max_set: int,
) -> GeneSets:
    sets = []
    for go_id, (name, aspect, genes) in terms.items():
        if min_set <= len(genes) <= max_set:
            sets.append(
                GeneSet(id=go_id, name=name, genes=sorted(genes), description=aspect)
            )
    sets.sort(key=lambda s: s.id)
    return GeneSets(
        sets=sets,
        metadata={
            "source": "UniProt",
            "collection": "GO",
            "organism": organism,
            "version": version,
            "min_set": min_set,
            "max_set": max_set,
        },
    )


# GMT has only one free-text column. To keep BOTH the readable GO name and the
# aspect (BP/CC/MF) through a GMT round-trip, encode them as "ASPECT|name" in
# that column and decode on load (aspect is 2 letters, never contains "|").
def _encode_for_gmt(gs: GeneSets) -> GeneSets:
    return GeneSets(
        sets=[
            GeneSet(s.id, s.name, list(s.genes), f"{s.description}|{s.name}")
            for s in gs
        ],
        metadata=dict(gs.metadata),
    )


def _decode_from_gmt(gs: GeneSets) -> GeneSets:
    out = []
    for s in gs:
        d = s.description
        if "|" in d and d.split("|", 1)[0] in {"BP", "CC", "MF"}:
            aspect, name = d.split("|", 1)
            out.append(GeneSet(s.id, name, list(s.genes), aspect))
        else:
            out.append(s)
    return GeneSets(sets=out, metadata=dict(gs.metadata))


# ---------------------------------------------------------------------------
# Live build from UniProt
# ---------------------------------------------------------------------------
def _uniprot_rows(taxon: int, *, max_pages: int = 40) -> list[str]:
    params = {
        "query": f"taxonomy_id:{taxon}",
        "format": "tsv",
        "fields": "accession,xref_veupathdb,go_p,go_c,go_f",
        "size": "500",
    }
    url: str | None = _UNIPROT + "?" + urllib.parse.urlencode(params)
    rows: list[str] = []
    pages = 0
    while url and pages < max_pages:
        req = urllib.request.Request(url, headers={"Accept": "text/plain"})
        with urllib.request.urlopen(req, timeout=90) as resp:  # noqa: S310
            text = resp.read().decode("utf-8")
            link = resp.headers.get("Link", "")
        lines = text.splitlines()
        rows.extend(lines[1:] if pages == 0 else lines[1:])
        m = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = m.group(1) if m else None
        pages += 1
    return rows


def build_gene_sets(
    organism: str,
    *,
    min_set: int = 5,
    max_set: int = 500,
    save_to: str | Path | None = None,
) -> GeneSets:
    """Build GO gene sets for an organism live from UniProt.

    Fetches, per protein, its PlasmoDB/VEuPathDB gene ID and GO terms (split by
    aspect), and returns a :class:`~plasmopack.io.GeneSets`. Pass ``save_to`` to
    also write a ``.gmt`` snapshot.

    Parameters
    ----------
    organism
        Registry key or alias (e.g. ``"pberghei"``).
    min_set, max_set
        Keep only GO terms with this many genes.
    save_to
        Optional path to write the resulting GMT (plus provenance sidecar).
    """
    org = get_organism(organism)
    rows = _uniprot_rows(org.species_taxonomy_id)
    terms = _parse_uniprot_tsv(rows, org.gene_id_prefix)
    import datetime as _dt

    version = f"UniProt-{_dt.date.today().isoformat()}"
    gs = _terms_to_genesets(
        terms, organism=org.key, version=version, min_set=min_set, max_set=max_set
    )
    if save_to is not None:
        write_gmt(_encode_for_gmt(gs), save_to)  # keep names + aspect in the GMT
    return gs


# ---------------------------------------------------------------------------
# Loading bundled / cached snapshots
# ---------------------------------------------------------------------------
_ASPECTS = ("GO", "BP", "CC", "MF")


def _bundled_path(organism: str) -> Path | None:
    try:
        p = files("plasmopack.datasets").joinpath("data", f"{organism}_GO.gmt")
        return Path(str(p)) if p.is_file() else None
    except (FileNotFoundError, ModuleNotFoundError):
        return None


def list_bundled_gene_sets() -> list[str]:
    """List organisms that ship with a bundled GO snapshot."""
    return [o for o in list_organisms() if _bundled_path(o) is not None]


def gene_sets(organism: str, collection: str = "GO") -> GeneSets:
    """Load default gene sets for an organism.

    Parameters
    ----------
    organism
        Registry key or alias (e.g. ``"pberghei"``).
    collection
        ``"GO"`` (all aspects) or a single aspect ``"BP"``/``"CC"``/``"MF"``.

    Returns
    -------
    GeneSets
        The bundled snapshot (filtered to the aspect if a single one is asked).

    Raises
    ------
    ValueError
        If ``collection`` is not recognised.
    FileNotFoundError
        If no bundled snapshot exists for the organism (suggests build_gene_sets).
    """
    if collection not in _ASPECTS:
        raise ValueError(f"collection must be one of {_ASPECTS}; got {collection!r}.")
    org = get_organism(organism)
    path = _bundled_path(org.key)
    if path is None:
        avail = ", ".join(list_bundled_gene_sets()) or "(none)"
        raise FileNotFoundError(
            f"No bundled gene sets for {org.key!r}. Bundled: {avail}. "
            f"Build them with pp.datasets.build_gene_sets({org.key!r}), or pass a "
            f"GeneSets object / .gmt path to enrich()."
        )
    gs = _decode_from_gmt(read_gmt(path))  # restore readable names + aspect
    if collection == "GO":
        return gs
    kept = [s for s in gs if s.description == collection]
    meta = dict(gs.metadata)
    meta["collection"] = collection
    return GeneSets(sets=kept, metadata=meta)
