"""Public UniProt client — ``pp.db.uniprot``.

A thin, stable surface over the private :mod:`plasmopack.adapters.uniprot`.
Callers get on-disk caching by default (offline after the first fetch) without
having to wire up the cache themselves.
"""

from __future__ import annotations

from collections.abc import Sequence

from plasmopack._utils.cache import DiskCache
from plasmopack.adapters import uniprot as _adapter
from plasmopack.adapters.uniprot import UniProtRecord
from plasmopack.datasets import get_organism

__all__ = ["UniProtRecord", "fetch", "search"]


def fetch(
    accessions: str | Sequence[str],
    *,
    use_cache: bool = True,
    cache_dir: str | None = None,
) -> list[UniProtRecord]:
    """Fetch UniProtKB records for one or more accessions.

    Parameters
    ----------
    accessions
        A single accession string or a sequence of them.
    use_cache
        If True (default), responses are cached on disk so repeated calls are
        offline and reproducible.
    cache_dir
        Optional cache directory. Defaults to ``$PLASMOPACK_CACHE_DIR`` or
        ``~/.plasmopack_cache``.

    Returns
    -------
    list[UniProtRecord]
        One record per input accession, in the same order.

    Examples
    --------
    >>> import plasmopack as pp
    >>> recs = pp.db.uniprot.fetch("Q8I3H7")           # doctest: +SKIP
    >>> recs[0].organism                                # doctest: +SKIP
    'Plasmodium falciparum (isolate 3D7)'
    """
    if isinstance(accessions, str):
        accessions = [accessions]

    cache = DiskCache(cache_dir, namespace="uniprot") if use_cache else None
    return _adapter.fetch(accessions, cache=cache)


def search(
    gene: str,
    *,
    organism: str | None = None,
    limit: int = 25,
    use_cache: bool = True,
    cache_dir: str | None = None,
) -> list[UniProtRecord]:
    """Search UniProtKB by gene name/symbol, optionally scoped to an organism.

    Use this when you have a gene *name* (e.g. ``"AMA1"``) rather than a UniProt
    accession. For a known accession, use :func:`fetch` instead.

    Parameters
    ----------
    gene
        Gene name or symbol, e.g. ``"AMA1"``.
    organism
        Optional organism key or alias from the registry (e.g. ``"pfalciparum"``
        or ``"pf"``). When given, results are restricted to that organism's
        NCBI taxon.
    limit
        Maximum number of records to return.
    use_cache
        Cache the query result on disk (offline + reproducible thereafter).
    cache_dir
        Optional cache directory override.

    Returns
    -------
    list[UniProtRecord]
        Matching records (possibly empty).

    Examples
    --------
    >>> import plasmopack as pp
    >>> hits = pp.db.uniprot.search("AMA1", organism="pfalciparum")   # doctest: +SKIP
    >>> hits[0].gene_names                                             # doctest: +SKIP
    ['AMA1']
    """
    query = f"gene:{gene}"
    if organism is not None:
        # UniProt indexes proteins at the species level and ``taxonomy_id:``
        # matches a taxon plus all its descendant strains, so the species
        # taxon is the correct, most inclusive filter.
        taxon = get_organism(organism).species_taxonomy_id
        query += f" AND taxonomy_id:{taxon}"

    cache = DiskCache(cache_dir, namespace="uniprot") if use_cache else None
    return _adapter.search(query, limit=limit, cache=cache)
