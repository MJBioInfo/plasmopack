"""Thin AnnData read/write wrappers.

These exist so users have a single, discoverable entry point
(``pp.io.read_h5ad``) consistent with the rest of the API, rather than reaching
for ``anndata`` directly. They add nothing surprising — just a stable surface
we can extend later (e.g. attaching the organism registry on read).
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
from anndata import AnnData


def read_h5ad(path: str | Path) -> AnnData:
    """Read an ``.h5ad`` file into an :class:`~anndata.AnnData`.

    Parameters
    ----------
    path
        Path to the ``.h5ad`` file.
    """
    return ad.read_h5ad(Path(path))


def write_h5ad(adata: AnnData, path: str | Path) -> Path:
    """Write an :class:`~anndata.AnnData` to ``.h5ad``.

    Parameters
    ----------
    adata
        The object to write.
    path
        Destination path.

    Returns
    -------
    Path
        The path written to.
    """
    path = Path(path)
    adata.write_h5ad(path)
    return path
