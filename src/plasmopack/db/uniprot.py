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

__all__ = ["UniProtRecord", "fetch"]


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
