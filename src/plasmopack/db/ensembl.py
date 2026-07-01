"""Public Ensembl client — ``pp.db.ensembl``.

A thin, stable surface over the private :mod:`plasmopack.adapters.ensembl`.
Works for *Plasmodium* and other Ensembl Genomes species. Caching is on by
default so repeated lookups are offline and reproducible.
"""

from __future__ import annotations

from collections.abc import Sequence

from plasmopack._utils.cache import DiskCache
from plasmopack.adapters import ensembl as _adapter
from plasmopack.adapters.ensembl import EnsemblGene

__all__ = ["EnsemblGene", "lookup"]


def lookup(
    gene_ids: str | Sequence[str],
    *,
    use_cache: bool = True,
    cache_dir: str | None = None,
) -> list[EnsemblGene]:
    """Look up gene records from Ensembl by stable gene id.

    Parameters
    ----------
    gene_ids
        A single gene id or a sequence of them, e.g. ``"PF3D7_0206800"``.
    use_cache
        If True (default), cache responses on disk (offline + reproducible).
    cache_dir
        Optional cache directory override.

    Returns
    -------
    list[EnsemblGene]
        One record per input id, in the same order.

    Examples
    --------
    >>> import plasmopack as pp
    >>> genes = pp.db.ensembl.lookup("PF3D7_0206800")   # doctest: +SKIP
    >>> genes[0].symbol                                  # doctest: +SKIP
    'MSP2'
    """
    if isinstance(gene_ids, str):
        gene_ids = [gene_ids]

    cache = DiskCache(cache_dir, namespace="ensembl") if use_cache else None
    return _adapter.lookup(gene_ids, cache=cache)
