"""UniProt adapter (private).

Fetches UniProtKB protein records via the stable EBI REST API and parses the
JSON into a small typed :class:`UniProtRecord`. Parsing is separated from
fetching so it can be unit-tested with canned JSON and no network.

Durability features:
- responses are cached on disk (offline after first fetch),
- the transport is injectable, so tests never touch the network,
- only a handful of fields are extracted, insulating callers from UniProt's
  large and evolving response schema.

API reference: https://www.uniprot.org/help/api
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any

from plasmopack._utils.cache import DiskCache, make_key
from plasmopack._utils.http import http_get_json

_BASE_URL = "https://rest.uniprot.org/uniprotkb"

# A transport is anything that takes a URL and returns parsed JSON.
Transport = Callable[[str], Any]


@dataclass
class UniProtRecord:
    """A parsed subset of a UniProtKB entry.

    Attributes
    ----------
    accession
        Primary accession, e.g. ``"Q8I3H7"``.
    entry_name
        UniProtKB ID, e.g. ``"K7NX48_PLAF7"``.
    protein_name
        Recommended full protein name (empty string if none).
    gene_names
        Gene name(s) associated with the entry.
    organism
        Scientific organism name.
    taxon_id
        NCBI taxonomy id, if present.
    length
        Sequence length in residues, if present.
    function
        Text of the first FUNCTION comment, if present.
    """

    accession: str
    entry_name: str = ""
    protein_name: str = ""
    gene_names: list[str] = field(default_factory=list)
    organism: str = ""
    taxon_id: int | None = None
    length: int | None = None
    function: str | None = None


def _dig(data: Any, *keys: str | int) -> Any:
    """Safely walk nested dict/list structures; return None if any step misses."""
    current = data
    for key in keys:
        if isinstance(key, int):
            if isinstance(current, list) and -len(current) <= key < len(current):
                current = current[key]
            else:
                return None
        else:
            if isinstance(current, dict) and key in current:
                current = current[key]
            else:
                return None
    return current


def parse_record(data: dict[str, Any]) -> UniProtRecord:
    """Parse a UniProtKB JSON entry into a :class:`UniProtRecord`.

    Pure function — no I/O. Missing fields degrade to sensible empties rather
    than raising, since UniProt entries vary widely in completeness.
    """
    accession = str(_dig(data, "primaryAccession") or "")

    gene_names: list[str] = []
    genes = data.get("genes")
    if isinstance(genes, list):
        for gene in genes:
            name = _dig(gene, "geneName", "value")
            if isinstance(name, str) and name:
                gene_names.append(name)

    function: str | None = None
    comments = data.get("comments")
    if isinstance(comments, list):
        for comment in comments:
            if _dig(comment, "commentType") == "FUNCTION":
                text = _dig(comment, "texts", 0, "value")
                if isinstance(text, str):
                    function = text
                    break

    # Protein name: reviewed (Swiss-Prot) entries use recommendedName; many
    # Plasmodium proteins are unreviewed (TrEMBL) and carry the name under
    # submissionNames instead (older exports used submittedNames — accept both).
    protein_name = _dig(
        data, "proteinDescription", "recommendedName", "fullName", "value"
    )
    for fallback_key in ("submissionNames", "submittedNames"):
        if protein_name:
            break
        protein_name = _dig(
            data, "proteinDescription", fallback_key, 0, "fullName", "value"
        )

    taxon = _dig(data, "organism", "taxonId")
    length = _dig(data, "sequence", "length")

    return UniProtRecord(
        accession=accession,
        entry_name=str(_dig(data, "uniProtkbId") or ""),
        protein_name=str(protein_name or ""),
        gene_names=gene_names,
        organism=str(_dig(data, "organism", "scientificName") or ""),
        taxon_id=int(taxon) if isinstance(taxon, int) else None,
        length=int(length) if isinstance(length, int) else None,
        function=function,
    )


def fetch_one(
    accession: str,
    *,
    transport: Transport | None = None,
    cache: DiskCache | None = None,
) -> UniProtRecord:
    """Fetch and parse a single UniProt accession.

    Parameters
    ----------
    accession
        UniProtKB accession, e.g. ``"Q8I3H7"``.
    transport
        Callable returning parsed JSON for a URL. Defaults to the stdlib
        HTTP getter (resolved at call time so tests can monkeypatch it);
        tests inject a fake to stay offline.
    cache
        Optional :class:`DiskCache`. When given, a hit avoids the network
        entirely and a miss is stored after fetching.
    """
    # Resolve at call time (not as a default arg) so monkeypatching works.
    if transport is None:
        transport = http_get_json

    key = make_key("uniprot", accession)
    if cache is not None:
        cached = cache.get(key)
        if cached is not None:
            return parse_record(cached)

    url = f"{_BASE_URL}/{accession}.json"
    data = transport(url)
    if cache is not None:
        cache.set(key, data)
    return parse_record(data)


def fetch(
    accessions: Sequence[str],
    *,
    transport: Transport | None = None,
    cache: DiskCache | None = None,
) -> list[UniProtRecord]:
    """Fetch and parse several UniProt accessions, preserving input order."""
    return [fetch_one(acc, transport=transport, cache=cache) for acc in accessions]
