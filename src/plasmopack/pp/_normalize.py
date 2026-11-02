"""Normalization for AnnData.

Library-size normalization followed by an optional log1p transform — the
standard first step before most single-cell analysis. Non-destructive by
default: the result is written to a layer and ``adata.X`` is left untouched,
so raw counts remain available (DEG pseudo-bulk needs them).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from anndata import AnnData

from plasmopack._utils.provenance import record_provenance


def normalize(
    adata: AnnData,
    *,
    target_sum: float | None = 1e4,
    log1p: bool = True,
    layer: str | None = None,
    key_added: str = "lognorm",
) -> AnnData:
    """Library-size normalize (and optionally log1p) into a new layer.

    Each cell's counts are scaled so they sum to ``target_sum`` (or to the
    per-cell median total when ``target_sum`` is ``None``), then optionally
    ``log1p``-transformed. The result goes to ``adata.layers[key_added]``;
    ``adata.X`` is not modified.

    Parameters
    ----------
    adata
        Annotated data; modified in place (adds a layer) and returned.
    target_sum
        Target total counts per cell after scaling. ``None`` uses the median
        of the per-cell totals.
    log1p
        Apply ``log1p`` after scaling.
    layer
        Source layer to normalize (defaults to ``adata.X``).
    key_added
        Destination layer name.

    Returns
    -------
    AnnData
        The same object, for chaining.
    """
    x = adata.layers[layer] if layer is not None else adata.X
    x = x.astype(np.float64) if not sp.issparse(x) else x.astype(np.float64)

    totals = np.asarray(x.sum(axis=1)).ravel()
    totals = np.where(totals > 0, totals, 1.0)  # avoid divide-by-zero

    effective_target = (
        float(target_sum) if target_sum is not None else float(np.median(totals))
    )
    scale = effective_target / totals  # per-cell factor

    if sp.issparse(x):
        scaler = sp.diags(scale)
        normed = scaler @ x
        normed = normed.tocsr()
        if log1p:
            normed.data = np.log1p(normed.data)
    else:
        normed = x * scale[:, None]
        if log1p:
            normed = np.log1p(normed)

    adata.layers[key_added] = normed

    record_provenance(
        adata,
        "pp.normalize",
        params={
            "target_sum": target_sum,
            "effective_target": effective_target,
            "log1p": log1p,
            "source_layer": layer,
            "key_added": key_added,
        },
    )
    return adata
