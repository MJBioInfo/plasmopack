"""Provenance stamping.

Every plasmopack function that annotates an ``AnnData`` calls
:func:`record_provenance` so the object carries a durable, self-documenting
history of how it was produced: which function ran, with which parameters,
against which database/GMT version, when, and with which package version.

The history lives at ``adata.uns["plasmopack"]["history"]`` as a list of dicts.
This is the single mechanism that makes a saved ``.h5ad`` reproducible years
later — it is the gap we observed in existing manually-produced datasets.
"""

from __future__ import annotations

import datetime as _dt
import json
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from anndata import AnnData

# Key under adata.uns that namespaces everything plasmopack writes.
UNS_ROOT = "plasmopack"
HISTORY_KEY = "history"

# Records are stored as JSON strings, not raw dicts. A list of dicts cannot be
# written to .h5ad, but a list of strings round-trips cleanly. We serialise on
# write and parse on read so callers always work with plain dicts.


def _package_version() -> str:
    from plasmopack import __version__

    return __version__


def record_provenance(
    adata: AnnData,
    function: str,
    *,
    params: dict[str, Any] | None = None,
    sources: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Append a provenance record to ``adata.uns["plasmopack"]["history"]``.

    Parameters
    ----------
    adata
        The annotated data object being modified. Mutated in place.
    function
        Dotted name of the calling function, e.g. ``"tl.enrich"``.
    params
        User-facing parameters that affect the result (groupby, method, ...).
        Keep this JSON-serialisable so the object round-trips through ``.h5ad``.
    sources
        External data provenance, e.g. ``{"gmt_version": "PlasmoDB-68"}`` or
        ``{"apicotfdb_snapshot": "2024-01"}``.

    Returns
    -------
    dict
        The record that was appended (useful for testing).
    """
    record: dict[str, Any] = {
        "function": function,
        "params": dict(params or {}),
        "sources": dict(sources or {}),
        "timestamp": _dt.datetime.now(_dt.UTC).isoformat(),
        "plasmopack_version": _package_version(),
    }

    root = adata.uns.setdefault(UNS_ROOT, {})
    # adata.uns values can be loaded back as non-dict containers; be defensive.
    if not isinstance(root, dict):
        root = {}
        adata.uns[UNS_ROOT] = root

    # Store as a list of JSON strings so it survives the .h5ad round-trip.
    history = list(root.get(HISTORY_KEY, []))
    history.append(json.dumps(record))
    root[HISTORY_KEY] = history
    return record


def get_history(adata: AnnData) -> list[dict[str, Any]]:
    """Return the provenance history as a list of dicts (empty if none)."""
    root = adata.uns.get(UNS_ROOT, {})
    if not isinstance(root, dict):
        return []
    history = root.get(HISTORY_KEY, [])
    if history is None:
        return []
    parsed: list[dict[str, Any]] = []
    for item in history:
        parsed.append(item if isinstance(item, dict) else json.loads(item))
    return parsed
