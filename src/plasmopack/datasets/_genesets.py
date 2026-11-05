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


_TAGS = {"BP", "CC", "MF", "KEGG"}


def _decode_from_gmt(gs: GeneSets) -> GeneSets:
    out = []
    for s in gs:
        d = s.description
        if "|" in d and d.split("|", 1)[0] in _TAGS:
            tag, name = d.split("|", 1)
            out.append(GeneSet(s.id, name, list(s.genes), tag))
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


# KEGG organism codes for the bundled organisms.
_KEGG_CODE = {
    "pberghei": "pbe",
    "pfalciparum": "pfa",
    "pknowlesi": "pkn",
    "tgondii": "tgo",
}


def _http_text(url: str, *, timeout: int = 90) -> str:
    req = urllib.request.Request(url, headers={"Accept": "text/plain"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310
        text: str = resp.read().decode("utf-8")
    return text


def _uniprot_geneid_to_gene(
    taxon: int, prefix: str, *, max_pages: int = 40
) -> dict[str, str]:
    """Map NCBI GeneID -> VEuPathDB gene id (e.g. 55148354 -> PBANKA_...)."""
    params = {
        "query": f"taxonomy_id:{taxon}",
        "format": "tsv",
        "fields": "xref_geneid,xref_veupathdb",
        "size": "500",
    }
    url: str | None = _UNIPROT + "?" + urllib.parse.urlencode(params)
    pb = re.compile(rf"{re.escape(prefix)}\w+")
    out: dict[str, str] = {}
    pages = 0
    while url and pages < max_pages:
        req = urllib.request.Request(url, headers={"Accept": "text/plain"})
        with urllib.request.urlopen(req, timeout=90) as resp:  # noqa: S310
            text = resp.read().decode("utf-8")
            link = resp.headers.get("Link", "")
        for ln in text.splitlines()[1:]:
            parts = ln.split("\t")
            if len(parts) < 2:
                continue
            geneids, xref = parts[0], parts[1]
            m = pb.search(xref)
            if not m:
                continue
            for gid in geneids.split(";"):
                gid = gid.strip()
                if gid:
                    out[gid] = m.group(0)
        nx = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = nx.group(1) if nx else None
        pages += 1
    return out


def _build_kegg(org, *, min_set: int, max_set: int) -> GeneSets:  # type: ignore[no-untyped-def]
    import datetime as _dt

    code = _KEGG_CODE.get(org.key)
    if code is None:
        raise ValueError(
            f"No KEGG organism code for {org.key!r}. KEGG is available for: "
            f"{', '.join(_KEGG_CODE)}."
        )
    links = _http_text(f"https://rest.kegg.jp/link/{code}/pathway")
    names = {
        p.split("\t")[0]: p.split("\t")[1]
        for p in _http_text(f"https://rest.kegg.jp/list/pathway/{code}").splitlines()
        if "\t" in p
    }
    path_genes: dict[str, set[str]] = {}
    for line in links.splitlines():
        if "\t" not in line:
            continue
        p, g = line.split("\t")
        path_genes.setdefault(p.replace("path:", ""), set()).add(g.split(":", 1)[-1])

    # KEGG gene ids are either the VEuPathDB locus tag already (e.g. PF3D7_...)
    # or a numeric NCBI GeneID (e.g. pbe). Map numerics via UniProt; pass tags
    # through directly.
    all_genes = {g for s in path_genes.values() for g in s}
    tag_re = re.compile(rf"^{re.escape(org.gene_id_prefix)}\w+$")
    needs_map = any(not tag_re.match(g) for g in all_genes)
    geneid_map = (
        _uniprot_geneid_to_gene(org.species_taxonomy_id, org.gene_id_prefix)
        if needs_map
        else {}
    )

    def _to_gene(g: str) -> str | None:
        if tag_re.match(g):
            return g
        return geneid_map.get(g)

    sets = []
    suffix = re.compile(
        r"\s*-\s*" + re.escape(org.scientific_name.split(" (")[0]) + r".*$"
    )
    for p, genes in path_genes.items():
        mapped = sorted({v for g in genes if (v := _to_gene(g)) is not None})
        if min_set <= len(mapped) <= max_set:
            nm = suffix.sub("", names.get(p, p)).strip() or p
            sets.append(GeneSet(id=p, name=nm, genes=mapped, description="KEGG"))
    sets.sort(key=lambda s: s.id)
    return GeneSets(
        sets=sets,
        metadata={
            "source": "KEGG",
            "collection": "KEGG",
            "organism": org.key,
            "version": f"KEGG-{_dt.date.today().isoformat()}",
            "min_set": min_set,
            "max_set": max_set,
        },
    )


def build_gene_sets(
    organism: str,
    *,
    collection: str = "GO",
    min_set: int = 5,
    max_set: int = 500,
    save_to: str | Path | None = None,
) -> GeneSets:
    """Build gene sets for an organism live from UniProt (GO) or KEGG.

    For ``collection="GO"`` fetches, per protein, its PlasmoDB/VEuPathDB gene ID
    and GO terms by aspect. For ``collection="KEGG"`` fetches KEGG pathways and
    maps their genes to VEuPathDB gene IDs via UniProt. Pass ``save_to`` to also
    write a ``.gmt`` snapshot.

    Parameters
    ----------
    organism
        Registry key or alias (e.g. ``"pberghei"``).
    collection
        ``"GO"`` or ``"KEGG"``.
    min_set, max_set
        Keep only sets with this many genes.
    save_to
        Optional path to write the resulting GMT (plus provenance sidecar).
    """
    org = get_organism(organism)
    if collection == "KEGG":
        gs = _build_kegg(org, min_set=min_set, max_set=max_set)
    elif collection == "GO":
        import datetime as _dt

        rows = _uniprot_rows(org.species_taxonomy_id)
        terms = _parse_uniprot_tsv(rows, org.gene_id_prefix)
        version = f"UniProt-{_dt.date.today().isoformat()}"
        gs = _terms_to_genesets(
            terms, organism=org.key, version=version, min_set=min_set, max_set=max_set
        )
    else:
        raise ValueError(f"collection must be 'GO' or 'KEGG'; got {collection!r}.")
    if save_to is not None:
        write_gmt(_encode_for_gmt(gs), save_to)  # keep names (+aspect) in the GMT
    return gs


# ---------------------------------------------------------------------------
# Loading bundled / cached snapshots
# ---------------------------------------------------------------------------
# Collections: GO (BP/CC/MF aspects, one file) and KEGG (one file).
_COLLECTIONS = ("GO", "BP", "CC", "MF", "KEGG")
_FILE_FOR = {"GO": "GO", "BP": "GO", "CC": "GO", "MF": "GO", "KEGG": "KEGG"}


def _bundled_path(organism: str, collection: str = "GO") -> Path | None:
    kind = _FILE_FOR.get(collection, "GO")
    try:
        p = files("plasmopack.datasets").joinpath("data", f"{organism}_{kind}.gmt")
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
        ``"GO"`` (all GO aspects), a single aspect ``"BP"``/``"CC"``/``"MF"``,
        or ``"KEGG"`` (metabolic pathways).

    Returns
    -------
    GeneSets
        The bundled snapshot (filtered to the aspect for a single GO aspect).

    Raises
    ------
    ValueError
        If ``collection`` is not recognised.
    FileNotFoundError
        If no bundled snapshot exists for the organism (suggests build_gene_sets).
    """
    if collection not in _COLLECTIONS:
        raise ValueError(
            f"collection must be one of {_COLLECTIONS}; got {collection!r}."
        )
    org = get_organism(organism)
    path = _bundled_path(org.key, collection)
    if path is None:
        avail = ", ".join(list_bundled_gene_sets()) or "(none)"
        raise FileNotFoundError(
            f"No bundled {collection} gene sets for {org.key!r}. Bundled GO: {avail}. "
            f"Build them with pp.datasets.build_gene_sets({org.key!r}, "
            f"collection={collection!r}), or pass a GeneSets / .gmt path to enrich()."
        )
    gs = _decode_from_gmt(read_gmt(path))  # restore readable names (+aspect)
    if collection in ("GO", "KEGG"):
        return gs
    kept = [s for s in gs if s.description == collection]  # one GO aspect
    meta = dict(gs.metadata)
    meta["collection"] = collection
    return GeneSets(sets=kept, metadata=meta)
