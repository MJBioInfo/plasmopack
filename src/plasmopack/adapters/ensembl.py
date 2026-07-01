"""Ensembl adapter (private).

Looks up gene records via the Ensembl REST API (``rest.ensembl.org``), which
serves Ensembl Genomes species — including *Plasmodium* and other VEuPathDB
organisms — through the same clean interface. One GET per gene returns
location, symbol, description, and biotype.

Follows the same durability pattern as the UniProt adapter: pure ``parse_gene``
separated from fetching, injectable transport, on-disk caching.

API reference: https://rest.ensembl.org/documentation/info/lookup
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

from plasmopack._utils.cache import DiskCache, make_key
from plasmopack._utils.http import HTTPError, http_get_json

_BASE_URL = "https://rest.ensembl.org"

Transport = Callable[[str], Any]


@dataclass
class EnsemblGene:
    """A parsed subset of an Ensembl gene lookup.

    Attributes
    ----------
    id
        Ensembl / stable gene id, e.g. ``"PF3D7_0206800"``.
    symbol
        Display name / gene symbol, e.g. ``"MSP2"`` (empty if none).
    description
        Free-text description, often with source cross-reference.
    biotype
        e.g. ``"protein_coding"``.
    species
        Ensembl species key, e.g. ``"plasmodium_falciparum"``.
    chromosome
        Sequence region name (chromosome/contig).
    start, end
        1-based genomic coordinates (None if absent).
    strand
        ``1`` or ``-1`` (None if absent).
    assembly
        Assembly name, e.g. ``"GCA000002765v3"``.
    """

    id: str
    symbol: str = ""
    description: str = ""
    biotype: str = ""
    species: str = ""
    chromosome: str = ""
    start: int | None = None
    end: int | None = None
    strand: int | None = None
    assembly: str = ""


def parse_gene(data: dict[str, Any]) -> EnsemblGene:
    """Parse an Ensembl lookup JSON object into an :class:`EnsemblGene`.

    Pure function — no I/O. Missing fields degrade to empties.
    """

    def _int(value: Any) -> int | None:
        return int(value) if isinstance(value, int) else None

    return EnsemblGene(
        id=str(data.get("id") or ""),
        symbol=str(data.get("display_name") or ""),
        description=str(data.get("description") or ""),
        biotype=str(data.get("biotype") or ""),
        species=str(data.get("species") or ""),
        chromosome=str(data.get("seq_region_name") or ""),
        start=_int(data.get("start")),
        end=_int(data.get("end")),
        strand=_int(data.get("strand")),
        assembly=str(data.get("assembly_name") or ""),
    )


def lookup_one(
    gene_id: str,
    *,
    transport: Transport | None = None,
    cache: DiskCache | None = None,
) -> EnsemblGene:
    """Look up a single gene id via Ensembl REST.

    Raises
    ------
    ValueError
        If Ensembl reports the id is invalid/unknown (HTTP 400).
    """
    if transport is None:
        transport = http_get_json

    key = make_key("ensembl_lookup", gene_id)
    if cache is not None:
        cached = cache.get(key)
        if cached is not None:
            return parse_gene(cached)

    url = f"{_BASE_URL}/lookup/id/{gene_id}"
    try:
        data = transport(url)
    except HTTPError as exc:
        if exc.status in (400, 404):
            raise ValueError(
                f"Ensembl has no gene with id {gene_id!r} "
                f"(check the identifier and that the species is in Ensembl)."
            ) from exc
        raise
    if cache is not None:
        cache.set(key, data)
    return parse_gene(data)


def lookup(
    gene_ids: Sequence[str],
    *,
    transport: Transport | None = None,
    cache: DiskCache | None = None,
) -> list[EnsemblGene]:
    """Look up several gene ids, preserving input order."""
    return [lookup_one(gid, transport=transport, cache=cache) for gid in gene_ids]
